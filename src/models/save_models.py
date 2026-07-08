"""Phase 11 prep — fit and persist deployable per-drug models for the demo.

For each drug: fit XGBoost on ALL labeled genomes (deployment uses all data), fit an isotonic
calibration curve on group-out-of-fold predictions, and pick the VME-capped operating threshold on
those OOF predictions. Saves one bundle per drug: model + calibrator + threshold + feature columns.

Honest-performance numbers come from the lineage-held-out evaluation (summaries #4–6); these
deployment models simply use all available data, as is standard once a method is validated.

Usage:
    python -m src.models.save_models --config config/config.yaml
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import yaml
from sklearn.isotonic import IsotonicRegression
from sklearn.model_selection import StratifiedGroupKFold, cross_val_predict
from xgboost import XGBClassifier

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from src.models.train_panel import RULES, pick_threshold  # noqa: E402

FEATURES = REPO / "data/processed/thin_slice_cipro_features.csv"
LINEAGES = REPO / "data/processed/thin_slice_cipro_lineages.csv"
PANEL = REPO / "data/processed/panel_labels.csv"
MODELS = REPO / "results/models"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/config.yaml")
    args = ap.parse_args()
    cfg = yaml.safe_load((REPO / args.config).read_text())
    seed, vme_target = cfg["project"]["seed"], cfg["evaluation"]["vme_target"]

    feats = pd.read_csv(FEATURES, dtype={"genome_id": str}).set_index("genome_id")
    cols = list(feats.columns)
    lineage = {r["genome_id"]: r["lineage"] for r in csv.DictReader(LINEAGES.open())}
    labels = pd.read_csv(PANEL, dtype=str)
    MODELS.mkdir(parents=True, exist_ok=True)

    for drug in RULES:
        sub = labels[labels[drug].isin(["Resistant", "Susceptible"])]
        sub = sub[sub["genome_id"].isin(feats.index)]
        g = sub["genome_id"].values
        X = feats.loc[g].values.astype(int)
        y = (sub[drug].values == "Resistant").astype(int)
        groups = np.array([lineage[x] for x in g])

        spw = (y == 0).sum() / max((y == 1).sum(), 1)
        base = XGBClassifier(n_estimators=300, max_depth=4, learning_rate=0.1, subsample=0.9,
                             colsample_bytree=0.8, eval_metric="logloss", scale_pos_weight=spw,
                             random_state=seed, n_jobs=4)
        cv = StratifiedGroupKFold(n_splits=4, shuffle=True, random_state=seed)
        oof = cross_val_predict(base, X, y, groups=groups, cv=cv, method="predict_proba")[:, 1]
        iso = IsotonicRegression(out_of_bounds="clip").fit(oof, y)
        thr = pick_threshold(y, iso.transform(oof), vme_target)
        base.fit(X, y)  # deployment model on all data

        joblib.dump({"model": base, "calibrator": iso, "threshold": float(thr),
                     "feature_cols": cols, "drug": drug, "rules": RULES[drug],
                     "n_train": int(len(y)), "vme_target": vme_target},
                    MODELS / f"{drug}.joblib")
        print(f"saved {drug}: n={len(y)} R={int(y.sum())} thr={thr:.2f} -> results/models/{drug}.joblib")
    return 0


if __name__ == "__main__":
    sys.exit(main())
