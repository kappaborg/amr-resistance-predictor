"""Top-tier addition — MIC (minimum inhibitory concentration) regression.

Predicts the *level* of resistance (continuous MIC, the clinical gold standard) instead of only
binary R/S, per drug, from the determinant feature matrix. Evaluated with the field-standard metric:
ESSENTIAL AGREEMENT (EA) = fraction of predictions within ±1 two-fold dilution of the true MIC
(i.e. |log2(pred) − log2(true)| ≤ 1) — the CLSI/FDA criterion for genotype-based MIC prediction.
Also reports RMSE (in doubling dilutions) and Pearson r, all under the phylogeny-aware split.

    python -m src.models.mic_regression --config config/config.yaml
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import yaml
from scipy.stats import pearsonr
from sklearn.model_selection import StratifiedGroupKFold
from xgboost import XGBRegressor

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from src.models.train_panel import RULES  # noqa: E402 (drug list)

AMR = "https://www.bv-brc.org/api/genome_amr/"
LAB = 'eq(evidence,%22Laboratory%20Method%22)'
DRUG_ENC = {"meropenem": "meropenem", "gentamicin": "gentamicin", "ciprofloxacin": "ciprofloxacin",
            "trimethoprim_sulfamethoxazole": "trimethoprim%2Fsulfamethoxazole", "cefoxitin": "cefoxitin"}
FEATURES = REPO / "data/processed/thin_slice_cipro_features.csv"
LINEAGES = REPO / "data/processed/thin_slice_cipro_lineages.csv"
OUT = REPO / "results/reports/summary_19_mic_regression.md"


def _parse_mic(m) -> tuple[float, str] | None:
    """(log2 MIC, censor) from a BV-BRC measurement string. censor: 'right' (>, MIC is a lower bound),
    'left' (<, MIC is an upper bound), or 'exact'."""
    s = str(m).strip()
    # STRICT vs NON-STRICT matters and must not be collapsed (audit, decision #60):
    #   ">=V"  -> true >= V      -> agreement needs pred >= V-1
    #   ">V"   -> true >= V+1    -> agreement needs pred >= V      (one dilution stricter)
    # Treating ">" with the ">=" rule scored 40.5% of ciprofloxacin measurements one full
    # doubling dilution too leniently, inflating Essential Agreement.
    if s.startswith(">="):
        censor = "right"
    elif s.startswith(">"):
        censor = "right_strict"
    elif s.startswith("<="):
        censor = "left"
    elif s.startswith("<"):
        censor = "left_strict"
    else:
        censor = "exact"
    try:
        v = float(re.sub(r"[<>=]", "", s).strip())
    except ValueError:
        return None
    return (np.log2(v), censor) if v > 0 else None


def fetch_mic_with_censor(drug_enc: str) -> dict[str, tuple[float, str]]:
    """genome_id -> (median log2(MIC mg/L), censor). Censor = the sign shared by a majority of the
    genome's measurements ('exact' if mixed), so right/left-censored (off-scale) MICs are flagged and
    can be scored honestly instead of collapsed onto the boundary dilution."""
    q = (f"and(eq(taxon_id,573),{LAB},eq(antibiotic,%22{drug_enc}%22))"
         f"&select(genome_id,measurement,measurement_unit)&limit(25000)")
    r = None
    for _ in range(4):
        r = requests.get(f"{AMR}?{q}", headers={"Accept": "application/json"}, timeout=180)
        if r.status_code == 200 and isinstance(r.json(), list):
            break
    byg: dict[str, list] = defaultdict(list)
    for x in r.json():
        m, unit = x.get("measurement"), x.get("measurement_unit")
        if not m or unit != "mg/L":
            continue
        p = _parse_mic(m)
        if p:
            byg[x["genome_id"]].append(p)
    out = {}
    for g, vs in byg.items():
        vals = [v for v, _ in vs]
        cens = [c for _, c in vs]
        maj = max(set(cens), key=cens.count)
        censor = maj if cens.count(maj) > len(cens) / 2 else "exact"
        out[g] = (float(np.median(vals)), censor)
    return out


def fetch_mic(drug_enc: str) -> dict[str, float]:
    """genome_id -> median log2(MIC mg/L) (value only; used for training and by esm2_mic)."""
    return {g: v for g, (v, _) in fetch_mic_with_censor(drug_enc).items()}


def essential_agreement(pred, true, censor) -> float:
    """CLSI/FDA Essential Agreement, CENSOR-AWARE. exact: within ±1 doubling dilution. right-censored
    (true ≥ V): a prediction ≥ V-1 is consistent with some true value within ±1, so it agrees.
    left-censored (true ≤ V): agrees if pred ≤ V+1."""
    pred, true = np.asarray(pred), np.asarray(true)
    agree = np.abs(pred - true) <= 1.0
    c = np.asarray([str(x) for x in censor])
    #   >=V : true >= V     -> pred >= V-1
    #   >V  : true >= V+1   -> pred >= V
    #   <=V : true <= V     -> pred <= V+1
    #   <V  : true <= V-1   -> pred <= V
    agree = np.where(c == "right",        pred >= true - 1.0, agree)
    agree = np.where(c == "right_strict", pred >= true,       agree)
    agree = np.where(c == "left",         pred <= true + 1.0, agree)
    agree = np.where(c == "left_strict",  pred <= true,       agree)
    return float(np.mean(agree))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/config.yaml")
    args = ap.parse_args()
    cfg = yaml.safe_load((REPO / args.config).read_text())
    seed, test_frac = cfg["project"]["seed"], cfg["split"]["test_fraction"]

    feats = pd.read_csv(FEATURES, dtype={"genome_id": str}).set_index("genome_id")
    lineage = {r["genome_id"]: r["lineage"] for r in csv.DictReader(LINEAGES.open())}

    md = ["# Summary #19 — MIC Regression (predicting the resistance *level*)\n",
          "Per-drug regression of continuous **log2(MIC, mg/L)** from determinant features, under a full "
          "5-fold **phylogeny-aware** CV (pooled held-out predictions). **Essential Agreement (EA)** = "
          "predictions within ±1 two-fold dilution (CLSI/FDA), scored **censor-aware**: right-/left-"
          "censored (off-scale `>`/`<`) MICs are treated as bounds, not collapsed onto the boundary "
          "dilution, so EA is not inflated. RMSE in doubling dilutions; Pearson r.\n",
          "| Drug | n (censored) | EA (±1 dilution) | RMSE (dilutions) | Pearson r |",
          "|---|---|---|---|---|"]
    results = {}
    print(f"{'drug':30s}{'n':>6}{'EA':>8}{'RMSE':>8}{'r':>7}")
    for drug, enc in DRUG_ENC.items():
        mic = fetch_mic_with_censor(enc)
        g = [x for x in feats.index if x in mic]
        if len(g) < 100:
            print(f"{drug:30s} too few MIC values ({len(g)})")
            continue
        X = feats.loc[g].values.astype(int)
        y = np.array([mic[x][0] for x in g])
        censor = np.array([mic[x][1] for x in g])
        groups = np.array([lineage[x] for x in g])
        n_cens = int((censor != "exact").sum())
        # full 5-fold lineage-grouped CV; pool held-out predictions (stable EA, not a single fold)
        ybin = (y >= np.median(y)).astype(int)
        P = np.full(len(y), np.nan)
        for tr, te in StratifiedGroupKFold(5, shuffle=True, random_state=seed).split(X, ybin, groups):
            model = XGBRegressor(n_estimators=400, max_depth=4, learning_rate=0.08, subsample=0.9,
                                 colsample_bytree=0.8, random_state=seed, n_jobs=4).fit(X[tr], y[tr])
            P[te] = model.predict(X[te])
        ea = essential_agreement(P, y, censor)                    # censor-aware EA
        rmse = float(np.sqrt(np.mean((P - y) ** 2)))
        r = float(pearsonr(P, y)[0]) if len(set(y)) > 1 else float("nan")
        results[drug] = {"n": len(g), "n_censored": n_cens, "essential_agreement": ea,
                         "rmse_dilutions": rmse, "pearson_r": r}
        md.append(f"| {drug} | {len(g)} ({n_cens} cens) | {ea:.1%} | {rmse:.2f} | {r:.3f} |")
        print(f"{drug:30s}{len(g):>6}{ea:>8.1%}{rmse:>8.2f}{r:>7.3f}")

    md.append("\n**Reading:** EA is the clinical gold-standard for MIC prediction — an EA around or "
              "above the ~90% CLSI expectation means the predicted MIC lands within one doubling "
              "dilution of the lab value, i.e. the model predicts not just *whether* a strain resists "
              "but *how strongly*. This is a richer, quantitative output than binary R/S, evaluated by "
              "the metric regulators actually use.")
    OUT.write_text("\n".join(md) + "\n")
    (REPO / "results/metrics/mic_metrics.json").write_text(json.dumps(results, indent=2))
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
