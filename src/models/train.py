"""Phase 6 — thin-slice model for ONE drug (ciprofloxacin), the Week-1 end-to-end proof.

Trains under the phylogeny-aware split and compares three things honestly:
  1. Known-gene RULES baseline  — predict Resistant iff any known cipro determinant is present
     (gyrA/parC point mutations, qnr*, aac(6')-Ib-cr, oqxAB, qepA). The transparent floor.
  2. Logistic-regression baseline — interpretable linear model on the determinant matrix.
  3. Gradient boosting (XGBoost) — the primary model.
All report ROC-AUC, PR-AUC, and — first — very-major / major error rates.

Usage:
    python -m src.models.train --config config/config.yaml
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
from sklearn.linear_model import LogisticRegression

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.evaluation.metrics import full_report  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
FEATURES = REPO / "data/processed/thin_slice_cipro_features.csv"
SPLIT = REPO / "data/processed/thin_slice_cipro_split.csv"
OUT = REPO / "results/metrics/thin_slice_cipro_metrics.json"

# Known ciprofloxacin resistance determinants (substring match on AMRFinderPlus element symbols).
CIPRO_RULES = ["gyrA_", "parC_", "gyrB_", "parE_", "qnr", "aac(6')-Ib-cr", "oqxA", "oqxB", "qepA"]


def rules_predict(feature_cols: list[str], X: np.ndarray) -> np.ndarray:
    """Baseline: Resistant (1) iff any known cipro determinant column is present."""
    idx = [i for i, c in enumerate(feature_cols) if any(tok in c for tok in CIPRO_RULES)]
    if not idx:
        return np.zeros(X.shape[0], dtype=int)
    return (X[:, idx].sum(axis=1) > 0).astype(int)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/config.yaml")
    args = ap.parse_args()
    cfg = yaml.safe_load((REPO / args.config).read_text())
    seed = cfg["project"]["seed"]

    feats = pd.read_csv(FEATURES).set_index("genome_id")
    split = pd.read_csv(SPLIT).set_index("genome_id")
    df = split.join(feats, how="inner")
    feature_cols = list(feats.columns)
    X = df[feature_cols].values.astype(int)
    y = (df["label"] == "Resistant").astype(int).values
    is_train = (df["fold"] == "train").values
    Xtr, Xte, ytr, yte = X[is_train], X[~is_train], y[is_train], y[~is_train]
    print(f"train {Xtr.shape} | test {Xte.shape} | {len(feature_cols)} determinants")

    results = {}

    # 1. rules baseline
    rp = rules_predict(feature_cols, Xte)
    results["rules_baseline"] = full_report(yte.tolist(), rp.tolist(), rp.astype(float).tolist())

    # 2. logistic regression
    lr = LogisticRegression(max_iter=2000, class_weight="balanced", random_state=seed)
    lr.fit(Xtr, ytr)
    lp, ls = lr.predict(Xte), lr.predict_proba(Xte)[:, 1]
    results["logreg"] = full_report(yte.tolist(), lp.tolist(), ls.tolist())

    # 3. gradient boosting (XGBoost)
    try:
        from xgboost import XGBClassifier
        spw = (ytr == 0).sum() / max((ytr == 1).sum(), 1)
        xgb = XGBClassifier(n_estimators=300, max_depth=4, learning_rate=0.1,
                            subsample=0.9, colsample_bytree=0.8, eval_metric="logloss",
                            scale_pos_weight=spw, random_state=seed, n_jobs=4)
        xgb.fit(Xtr, ytr)
        xp, xs = xgb.predict(Xte), xgb.predict_proba(Xte)[:, 1]
        results["xgboost"] = full_report(yte.tolist(), xp.tolist(), xs.tolist())
    except Exception as e:  # noqa: BLE001
        results["xgboost"] = {"error": str(e)}

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(results, indent=2))

    print(f"\n{'model':16s} {'ROC-AUC':>8} {'PR-AUC':>8} {'VME':>7} {'ME':>7}")
    for name, m in results.items():
        if "error" in m:
            print(f"{name:16s}  ERROR: {m['error']}")
            continue
        print(f"{name:16s} {m['roc_auc']:8.3f} {m['pr_auc']:8.3f} "
              f"{m['very_major_error_rate']:7.3f} {m['major_error_rate']:7.3f}")
    print(f"\nwrote {OUT}")
    print("VME = resistant called susceptible (dangerous); ME = susceptible called resistant.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
