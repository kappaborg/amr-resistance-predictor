"""Regulator-grade metrics table with calibration, lineage-clustered bootstrap CIs, and DeLong tests
— across ALL organisms.

For every organism × drug, under the phylogeny-aware split (pooled held-out predictions from a full
5-fold lineage-grouped CV; thresholds picked on each fold's TRAIN out-of-fold predictions), the model
is an **isotonic-calibrated XGBoost** and we report:

  Very-Major Error (VME, R->S)   FDA/CLSI target  <= 1.5%   (<= 3% tolerated with CI caveats)
  Major Error (ME, S->R)                          <= 3%
  Categorical Agreement (CA)                       >= 90%
  Brier score (calibration)       lower is better   (non-negotiable #7: always calibrate + report it)
  Essential Agreement (EA)                         >= 90%   (K. pneumoniae MIC models, where available)

Every rate carries a **lineage-clustered** bootstrap 95% CI (resample whole MLST lineages, not i.i.d.
genomes — the honest variance given clonal correlation, the project's own thesis). Each drug's ROC-AUC
is compared to the honest **organism-aware** known-gene rules baseline (`rules_for`) with a DeLong
paired test.

    python -m src.evaluation.clinical_rigor --config config/config.yaml            # all organisms
    python -m src.evaluation.clinical_rigor --organism paeruginosa
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
import pickle

from scipy import stats
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import matthews_corrcoef, roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold, cross_val_predict
from xgboost import XGBClassifier

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from src.app.registry import ORGANISMS, ORG_ORDER  # noqa: E402
from src.organism_pipeline import rules_for  # noqa: E402 (organism-aware baseline)
from src.models.train_panel import rules_predict, pick_threshold  # noqa: E402
from src.evaluation.cross_species_transfer import fetch_labels  # noqa: E402 (robust label fetch)

MIC_METRICS = REPO / "results/metrics/mic_metrics.json"     # K. pneumoniae EA (Summary #19)
OUT = REPO / "results/reports/summary_24_clinical_rigor.md"
FDA = {"vme": 0.015, "vme_tolerated": 0.03, "me": 0.03, "ca": 0.90, "ea": 0.90}


# ---------- DeLong test (fast algorithm, Sun & Xu 2014) for two correlated ROC AUCs ----------
def _compute_midrank(x):
    J = np.argsort(x)
    z = x[J]
    N = len(x)
    T = np.zeros(N, dtype=float)
    i = 0
    while i < N:
        j = i
        while j < N and z[j] == z[i]:
            j += 1
        T[i:j] = 0.5 * (i + j - 1) + 1
        i = j
    T2 = np.empty(N, dtype=float)
    T2[J] = T
    return T2


def _fast_delong(preds_sorted, m):
    n = preds_sorted.shape[1] - m
    k = preds_sorted.shape[0]
    tx = np.empty([k, m]); ty = np.empty([k, n]); tz = np.empty([k, m + n])
    for r in range(k):
        tx[r] = _compute_midrank(preds_sorted[r, :m])
        ty[r] = _compute_midrank(preds_sorted[r, m:])
        tz[r] = _compute_midrank(preds_sorted[r])
    aucs = tz[:, :m].sum(axis=1) / m / n - (m + 1.0) / 2.0 / n
    v01 = (tz[:, :m] - tx) / n
    v10 = 1.0 - (tz[:, m:] - ty) / m
    sx = np.cov(v01); sy = np.cov(v10)
    delongcov = sx / m + sy / n
    return aucs, delongcov


def delong_p(y_true, score1, score2):
    order = np.argsort(-np.asarray(y_true))          # positives (1) first
    yt = np.asarray(y_true)[order]
    m = int((yt == 1).sum())
    if m == 0 or m == len(yt):
        return np.array([float("nan"), float("nan")]), float("nan")
    preds = np.vstack([np.asarray(score1, float)[order], np.asarray(score2, float)[order]])
    aucs, cov = _fast_delong(preds, m)
    var = cov[0, 0] + cov[1, 1] - 2 * cov[0, 1]
    if var <= 0:
        return aucs, float("nan")
    z = (aucs[0] - aucs[1]) / np.sqrt(var)
    return aucs, float(2 * stats.norm.sf(abs(z)))


# ---------- metrics + lineage-clustered bootstrap ----------
def rates(y, pred):
    y, pred = np.asarray(y), np.asarray(pred)
    nR, nS = (y == 1).sum(), (y == 0).sum()
    vme = ((y == 1) & (pred == 0)).sum() / nR if nR else np.nan
    me = ((y == 0) & (pred == 1)).sum() / nS if nS else np.nan
    ca = (y == pred).mean()
    return vme, me, ca


def brier(y, score):
    return float(np.mean((np.asarray(score, float) - np.asarray(y, float)) ** 2))


PREVALENCES = (0.10, 0.30, 0.50)   # realistic local resistance prevalences for PPV/NPV


def ppv_npv(sens, spec, prev):
    """Prevalence-adjusted PPV/NPV from sensitivity (=1-VME) and specificity (=1-ME)."""
    ppv = sens * prev / (sens * prev + (1 - spec) * (1 - prev)) if (sens * prev + (1 - spec) * (1 - prev)) else np.nan
    npv = spec * (1 - prev) / (spec * (1 - prev) + (1 - sens) * prev) if (spec * (1 - prev) + (1 - sens) * prev) else np.nan
    return float(ppv), float(npv)


def boot_ci(y, pred, sraw, scal, groups, n_boot, seed):
    """LINEAGE-CLUSTERED bootstrap: resample whole MLST lineages (blocks), not i.i.d. genomes.
    Genomes within a lineage are clonally correlated, so an i.i.d. bootstrap would under-estimate
    variance — this gives honest CIs consistent with the project's population-structure argument.
    AUC uses RAW scores (rank-preserving); Brier uses the CALIBRATED scores."""
    rng = np.random.RandomState(seed)
    y, pred, sraw, scal, groups = (np.asarray(y), np.asarray(pred), np.asarray(sraw),
                                   np.asarray(scal), np.asarray(groups))
    uniq = np.unique(groups)
    idx_by = {g: np.where(groups == g)[0] for g in uniq}
    vs, ms, cs, aucs, brs = [], [], [], [], []
    for _ in range(n_boot):
        pick = rng.choice(uniq, len(uniq), replace=True)
        idx = np.concatenate([idx_by[g] for g in pick])
        if len(set(y[idx])) < 2:
            continue
        v, m, c = rates(y[idx], pred[idx])
        vs.append(v); ms.append(m); cs.append(c)
        aucs.append(roc_auc_score(y[idx], sraw[idx])); brs.append(brier(y[idx], scal[idx]))
    q = lambda a: (float(np.percentile(a, 2.5)), float(np.percentile(a, 97.5))) if a else (np.nan, np.nan)
    return {"vme": q(vs), "me": q(ms), "ca": q(cs), "auc": q(aucs), "brier": q(brs)}


def pooled_oof(X, y, groups, cols, rules_tokens, seed, vme_target):
    """Full 5-fold lineage-grouped CV; per-fold isotonic-calibrated XGBoost, threshold on TRAIN
    inner-OOF. Returns pooled held-out (y, RAW score, CALIBRATED score, R/S pred, rules pred, groups).
    RAW scores are for AUC/DeLong (isotonic ties would deflate AUC); CALIBRATED scores are for Brier
    and the R/S decision at the VME-capped threshold."""
    Y, Sr, Sc, P, R, G = [], [], [], [], [], []
    outer = StratifiedGroupKFold(5, shuffle=True, random_state=seed)
    inner = StratifiedGroupKFold(4, shuffle=True, random_state=seed)
    for tr, te in outer.split(X, y, groups):
        spw = (y[tr] == 0).sum() / max((y[tr] == 1).sum(), 1)
        mdl = XGBClassifier(n_estimators=300, max_depth=4, learning_rate=0.1, subsample=0.9,
                            colsample_bytree=0.8, eval_metric="logloss", scale_pos_weight=spw,
                            random_state=seed, n_jobs=4)
        oof = cross_val_predict(mdl, X[tr], y[tr], groups=groups[tr], cv=inner,
                                method="predict_proba")[:, 1]
        iso = IsotonicRegression(out_of_bounds="clip").fit(oof, y[tr])
        thr = pick_threshold(y[tr], iso.transform(oof), vme_target)
        mdl.fit(X[tr], y[tr])
        raw = mdl.predict_proba(X[te])[:, 1]
        cal = iso.transform(raw)
        Y.extend(y[te]); Sr.extend(raw); Sc.extend(cal); P.extend((cal >= thr).astype(int))
        R.extend(rules_predict(cols, X[te], rules_tokens)); G.extend(groups[te])
    return np.array(Y), np.array(Sr), np.array(Sc), np.array(P), np.array(R), np.array(G)


def load_labels(org_key):
    org = ORGANISMS[org_key]
    drugs = list(org["drugs"])
    if org["labels_file"]:
        rows = list(csv.DictReader((REPO / org["labels_file"]).open()))
        return {d: {r["genome_id"]: r[d] for r in rows if r.get(d) in ("Resistant", "Susceptible")}
                for d in drugs}
    return {d: fetch_labels(org["taxon"], org["drugs"][d]) for d in drugs}


def flag(ok):
    return "✅" if ok else ("❌" if ok is not None else "—")


def evaluate_organism(org_key, seed, vme_target, n_boot, ea_map, preds_store):
    org = ORGANISMS[org_key]
    feats = pd.read_csv(REPO / org["features"], dtype={"genome_id": str}).set_index("genome_id")
    lineage = {r["genome_id"]: r["lineage"]
               for r in csv.DictReader((REPO / org["lineages"]).open())}
    cols = list(feats.columns)
    labels = load_labels(org_key)
    out = {}
    print(f"\n[{org['display']}]  {'drug':30s}{'VME':>9}{'ME':>9}{'CA':>8}{'MCC':>7}{'Brier':>8}{'AUC vs rules':>16}{'DeLong p':>11}")
    for drug in org["drugs"]:
        lab = labels[drug]
        g = [x for x in feats.index if lab.get(x) in ("Resistant", "Susceptible")]
        y = np.array([1 if lab[x] == "Resistant" else 0 for x in g])
        if y.sum() < 20 or (y == 0).sum() < 20:
            continue
        X = feats.loc[g].values.astype(int)
        groups = np.array([lineage.get(x, x) for x in g])
        Y, Sr, Sc, P, R, G = pooled_oof(X, y, groups, cols, rules_for(org_key, drug), seed, vme_target)
        vme, me, ca = rates(Y, P)
        br = brier(Y, Sc)                                 # Brier on calibrated scores
        mcc = float(matthews_corrcoef(Y, P))              # robust to class imbalance
        sens, spec = 1 - vme, 1 - me
        ppvnpv = {f"{int(pr*100)}%": ppv_npv(sens, spec, pr) for pr in PREVALENCES}
        ci = boot_ci(Y, P, Sr, Sc, G, n_boot, seed)
        (auc_m, auc_r), dp = delong_p(Y, Sr, R.astype(float))   # AUC/DeLong on RAW scores
        ea_v = ea_map.get(drug, {}).get("essential_agreement") if org_key == "kpneu" else None
        out[drug] = {"n": int(len(Y)), "n_R": int(Y.sum()), "vme": float(vme), "vme_ci": ci["vme"],
                     "me": float(me), "me_ci": ci["me"], "ca": float(ca), "ca_ci": ci["ca"],
                     "mcc": mcc, "ppv_npv_at_prev": ppvnpv,
                     "brier": br, "brier_ci": ci["brier"], "auc_model": float(auc_m),
                     "auc_ci": ci["auc"], "auc_rules": float(auc_r), "delong_p": dp, "ea": ea_v,
                     "vme_pass_15": bool(vme <= FDA["vme"]), "vme_pass_3": bool(vme <= FDA["vme_tolerated"]),
                     "me_pass": bool(me <= FDA["me"]), "ca_pass": bool(ca >= FDA["ca"]),
                     "ea_pass": (bool(ea_v >= FDA["ea"]) if ea_v is not None else None)}
        # save pooled held-out predictions (calibrated score + rules + groups) for DCA / risk-coverage
        preds_store[(org_key, drug)] = {"y": Y, "score": Sc, "rules": R, "groups": G}
        dps = f"{dp:.1e}" if dp == dp else "n/a"
        print(f"{'':2s}{drug:30s}{vme:>8.1%}{me:>9.1%}{ca:>8.1%}{mcc:>7.2f}{br:>8.3f}"
              f"{auc_m:>8.3f} vs {auc_r:.3f}{dps:>11}")
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/config.yaml")
    ap.add_argument("--organism", choices=ORG_ORDER, help="default: all 8")
    ap.add_argument("--n-boot", type=int, default=2000)
    args = ap.parse_args()
    cfg = yaml.safe_load((REPO / args.config).read_text())
    seed, vme_target = cfg["project"]["seed"], cfg["evaluation"]["vme_target"]
    ea_map = json.loads(MIC_METRICS.read_text()) if MIC_METRICS.exists() else {}

    orgs = [args.organism] if args.organism else ORG_ORDER
    preds_store = {}
    allres = {ok: evaluate_organism(ok, seed, vme_target, args.n_boot, ea_map, preds_store)
              for ok in orgs}
    (REPO / "results/metrics/pooled_predictions.pkl").write_bytes(pickle.dumps(preds_store))

    # ---- report ----
    md = ["# Summary #24 — Regulator-Grade Metrics (all organisms): Calibration, CIs & DeLong\n",
          f"**Date:** 2026-07-17 · **all {len(orgs)} organisms** · isotonic-calibrated XGBoost · pooled "
          "held-out predictions from full 5-fold **lineage-grouped** CV (thresholds on each fold's TRAIN "
          f"out-of-fold predictions, VME target {vme_target:.0%}) · {args.n_boot} **lineage-clustered** "
          "bootstrap resamples.\n",
          "FDA/CLSI vocabulary + bars: **VME ≤ 1.5%** (≤3% tolerated), **ME ≤ 3%**, **CA ≥ 90%**, "
          "**EA ≥ 90%**; **Brier** = calibration (lower better). ROC-AUC vs the **organism-aware** "
          "known-gene rules baseline by a **DeLong** paired test. CIs resample whole MLST lineages "
          "(clonal blocks), the honest variance under population structure.\n"]
    for ok in orgs:
        res = allres[ok]
        if not res:
            continue
        md.append(f"\n## {ORGANISMS[ok]['display']}\n")
        md += ["| Drug | n (R) | VME % [CI] | ME % [CI] | CA % | MCC | Brier | ROC-AUC [CI] | Rules AUC | "
               "DeLong p | PPV/NPV @30% | VME≤1.5 | ME≤3 |",
               "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
        for d, x in res.items():
            dps = f"{x['delong_p']:.1e}" if x["delong_p"] == x["delong_p"] else "n/a"
            ppv, npv = x["ppv_npv_at_prev"]["30%"]
            md.append(
                f"| {d} | {x['n']} ({x['n_R']}) | {x['vme']:.1%} [{x['vme_ci'][0]:.1%},{x['vme_ci'][1]:.1%}] "
                f"| {x['me']:.1%} [{x['me_ci'][0]:.1%},{x['me_ci'][1]:.1%}] | {x['ca']:.1%} | {x['mcc']:.2f} "
                f"| {x['brier']:.3f} | {x['auc_model']:.3f} [{x['auc_ci'][0]:.3f},{x['auc_ci'][1]:.3f}] "
                f"| {x['auc_rules']:.3f} | {dps} | {ppv:.2f}/{npv:.2f} | {flag(x['vme_pass_15'])} "
                f"| {flag(x['me_pass'])} |")
    md.append(
        "\n**Reading.** Every organism × drug now carries a **calibrated** model (isotonic; Brier "
        "reported), **lineage-clustered** 95% CIs, and a **DeLong** significance test versus the "
        "organism-aware gene-lookup baseline — the project's central claim (\"what ML adds over a "
        "lookup\"), significance-tested across all organisms. Thresholds target the clinical **VME ≤ 3%** "
        "operating point; where the strict FDA 1.5% bar or CA ≥ 90% is not met it is shown honestly "
        "(the high ME / borderline VME on low-prevalence-R drugs is the safety-first over-calling "
        "trade-off, reported not hidden). **MCC** (Matthews correlation) is robust to the class "
        "imbalance that inflates accuracy/CA; **PPV/NPV** are shown at a realistic 30% local resistance "
        "prevalence (full 10/30/50% set in `clinical_rigor.json`) since predictive value depends on "
        "prevalence, not just the test's sensitivity/specificity.")
    OUT.write_text("\n".join(md) + "\n")
    (REPO / "results/metrics/clinical_rigor.json").write_text(json.dumps(allres, indent=2))
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
