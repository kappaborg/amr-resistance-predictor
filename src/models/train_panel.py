"""Phase 7+8 — per-drug models for the full panel (one model per antibiotic).

For each drug: build a lineage-aware split (StratifiedGroupKFold on MLST lineage, using only the
genomes labeled for that drug), then compare the expert known-gene RULES baseline, logistic
regression, and CALIBRATED XGBoost. Reports ROC-AUC, PR-AUC, and — first — VME/ME, plus Brier
score for calibration. Everything reuses the already-built feature matrix + MLST lineages (no new
annotation).

Calibration is lineage-clean: XGBoost is fit on a subset of TRAIN lineages and calibrated
(isotonic) on a disjoint subset of TRAIN lineages, then evaluated on the held-out TEST lineages.

Usage:
    python -m src.models.train_panel --config config/config.yaml
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
from sklearn.calibration import CalibratedClassifierCV
from sklearn.frozen import FrozenEstimator
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import brier_score_loss
from sklearn.model_selection import GroupShuffleSplit, StratifiedGroupKFold

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from src.evaluation.metrics import full_report  # noqa: E402

FEATURES = REPO / "data/processed/thin_slice_cipro_features.csv"
LINEAGES = REPO / "data/processed/thin_slice_cipro_lineages.csv"
PANEL = REPO / "data/processed/panel_labels.csv"
OUT_JSON = REPO / "results/metrics/panel_metrics.json"
OUT_MD = REPO / "results/reports/summary_05_week2_panel.md"

# Expert known-gene rules per drug (substring match on AMRFinderPlus element symbols).
# Grounded in the determinants actually present in the matrix; excludes intrinsic/irrelevant genes.
RULES = {
    "meropenem": ["blaKPC", "blaNDM", "blaOXA-48", "blaOXA-181", "blaOXA-232", "blaVIM", "blaIMP", "blaGES"],
    "gentamicin": ["aac(3)", "ant(2'')-Ia", "armA", "rmt"],  # aac(6') excluded (amikacin/tobramycin)
    "ciprofloxacin": ["gyrA_", "parC_", "gyrB_", "parE_", "qnr", "aac(6')-Ib-cr", "qepA"],
    "trimethoprim_sulfamethoxazole": ["sul1", "sul2", "sul3", "dfrA", "dfrB"],
    "cefoxitin": ["blaCMY", "blaDHA", "blaACT", "blaFOX", "blaMOX"],  # plasmid AmpC (ESBLs excluded)
}


def rules_predict(cols, X, tokens):
    idx = [i for i, c in enumerate(cols) if any(t in c for t in tokens)]
    if not idx:
        return np.zeros(X.shape[0], dtype=int)
    return (X[:, idx].sum(axis=1) > 0).astype(int)


def run_drug(drug, feats, lineage, labels, seed, test_frac):
    sub = labels[labels[drug].isin(["Resistant", "Susceptible"])].copy()
    sub = sub[sub["genome_id"].isin(feats.index)]
    g = sub["genome_id"].values
    X = feats.loc[g].values.astype(int)
    y = (sub[drug].values == "Resistant").astype(int)
    groups = np.array([lineage[x] for x in g])
    cols = list(feats.columns)

    n_splits = max(2, round(1 / test_frac))
    sgkf = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    tr_idx, te_idx = next(sgkf.split(X, y, groups))
    Xtr, Xte, ytr, yte, gtr = X[tr_idx], X[te_idx], y[tr_idx], y[te_idx], groups[tr_idx]

    res = {"n_train": int(len(tr_idx)), "n_test": int(len(te_idx)),
           "test_R": int(yte.sum()), "test_S": int((yte == 0).sum()),
           "train_lineages": int(len(set(gtr))), "test_lineages": int(len(set(groups[te_idx])))}

    # 1. rules baseline
    rp = rules_predict(cols, Xte, RULES[drug])
    res["rules"] = full_report(yte.tolist(), rp.tolist(), rp.astype(float).tolist())

    # 2. logistic regression
    lr = LogisticRegression(max_iter=2000, class_weight="balanced", random_state=seed).fit(Xtr, ytr)
    res["logreg"] = full_report(yte.tolist(), lr.predict(Xte).tolist(), lr.predict_proba(Xte)[:, 1].tolist())

    # 3. calibrated XGBoost — fit on a subset of TRAIN lineages, calibrate on a disjoint subset
    from xgboost import XGBClassifier
    spw = (ytr == 0).sum() / max((ytr == 1).sum(), 1)
    fit_idx, cal_idx = next(GroupShuffleSplit(n_splits=1, test_size=0.25,
                                              random_state=seed).split(Xtr, ytr, gtr))
    base = XGBClassifier(n_estimators=300, max_depth=4, learning_rate=0.1, subsample=0.9,
                         colsample_bytree=0.8, eval_metric="logloss", scale_pos_weight=spw,
                         random_state=seed, n_jobs=4)
    base.fit(Xtr[fit_idx], ytr[fit_idx])
    cal = CalibratedClassifierCV(FrozenEstimator(base), method="isotonic").fit(Xtr[cal_idx], ytr[cal_idx])
    xs = cal.predict_proba(Xte)[:, 1]
    xp = (xs >= 0.5).astype(int)
    rep = full_report(yte.tolist(), xp.tolist(), xs.tolist())
    rep["brier"] = float(brier_score_loss(yte, xs))
    res["xgboost_calibrated"] = rep
    return res


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/config.yaml")
    args = ap.parse_args()
    cfg = yaml.safe_load((REPO / args.config).read_text())
    seed, test_frac = cfg["project"]["seed"], cfg["split"]["test_fraction"]

    feats = pd.read_csv(FEATURES, dtype={"genome_id": str}).set_index("genome_id")
    lineage = {r["genome_id"]: r["lineage"] for r in csv.DictReader(LINEAGES.open())}
    labels = pd.read_csv(PANEL, dtype=str)

    results = {}
    for drug in RULES:
        print(f"\n=== {drug} ===")
        results[drug] = run_drug(drug, feats, lineage, labels, seed, test_frac)
        r = results[drug]
        print(f"test {r['n_test']} genomes ({r['test_R']}R/{r['test_S']}S), "
              f"{r['test_lineages']} lineages")
        print(f"{'model':22s} {'ROC':>6} {'PR':>6} {'VME':>6} {'ME':>6}")
        for m in ("rules", "logreg", "xgboost_calibrated"):
            d = r[m]
            print(f"{m:22s} {d['roc_auc']:6.3f} {d['pr_auc']:6.3f} "
                  f"{d['very_major_error_rate']:6.3f} {d['major_error_rate']:6.3f}")

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(results, indent=2))
    write_report(results)
    print(f"\nwrote {OUT_JSON}\nwrote {OUT_MD}")
    return 0


def write_report(results: dict) -> None:
    lines = ["# Summary #5 — Week 2: Full-Panel Per-Drug Models\n",
             "One model per antibiotic, reusing the already-annotated 1472-genome set (no new "
             "download/annotation). Each drug: lineage-aware split (StratifiedGroupKFold on MLST), "
             "expert rules baseline vs logistic regression vs calibrated XGBoost. **VME/ME lead.**\n",
             "VME = resistant called susceptible (dangerous miss); ME = susceptible called resistant.\n"]
    for drug, r in results.items():
        lines.append(f"\n## {drug}  (test {r['n_test']}: {r['test_R']}R/{r['test_S']}S, "
                     f"{r['test_lineages']} unseen lineages)\n")
        lines.append("| Model | ROC-AUC | PR-AUC | VME | ME | Brier |")
        lines.append("|---|---|---|---|---|---|")
        for m, name in [("rules", "Rules baseline"), ("logreg", "Logistic regression"),
                        ("xgboost_calibrated", "XGBoost (calibrated)")]:
            d = r[m]
            brier = f"{d['brier']:.3f}" if "brier" in d else "—"
            lines.append(f"| {name} | {d['roc_auc']:.3f} | {d['pr_auc']:.3f} | "
                         f"{d['very_major_error_rate']:.3f} | {d['major_error_rate']:.3f} | {brier} |")
    OUT_MD.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
