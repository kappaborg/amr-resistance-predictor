"""Fit and persist deployable per-drug models for the demo — all eight organisms.

For each organism × drug: fit XGBoost on ALL labelled genomes (deployment uses all data), fit an
isotonic calibrator on group-out-of-fold predictions, and pick the VME-capped operating threshold on
those OOF predictions. Saves one bundle per drug under results/models/<organism>/<drug>.joblib:
model + calibrator + threshold + feature columns + AMRFinderPlus organism string + rules.

Honest-performance numbers come from the lineage-held-out evaluation (summaries #4-6, #16-18); these
deployment models simply use all available data, as is standard once a method is validated.

    python -m src.models.save_models                 # all organisms
    python -m src.models.save_models --organism kpneu
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
from src.app.registry import ORGANISMS, ORG_ORDER  # noqa: E402
from src.evaluation.cross_species_transfer import fetch_labels  # noqa: E402 (robust label fetch)
from src.organism_pipeline import rules_for  # noqa: E402 (organism-aware known-gene baseline)
from src.models.train_panel import pick_threshold  # noqa: E402

MODELS = REPO / "results/models"
MIN_PER_CLASS = 20   # need enough of each class to calibrate + pick a threshold


def load_labels(org_key: str, org: dict, drugs: list[str]) -> dict[str, dict[str, str]]:
    """drug -> {genome_id: 'Resistant'|'Susceptible'}. Saved panel for K. pneumoniae, else live."""
    if org["labels_file"]:
        rows = list(csv.DictReader((REPO / org["labels_file"]).open()))
        return {d: {r["genome_id"]: r[d] for r in rows if r.get(d) in ("Resistant", "Susceptible")}
                for d in drugs}
    out = {}
    for d in drugs:
        out[d] = fetch_labels(org["taxon"], org["drugs"][d])
    return out


def save_organism(org_key: str, seed: int, vme_target: float) -> list[str]:
    org = ORGANISMS[org_key]
    feats = pd.read_csv(REPO / org["features"], dtype={"genome_id": str}).set_index("genome_id")
    cols = list(feats.columns)
    lineage = {r["genome_id"]: r["lineage"]
               for r in csv.DictReader((REPO / org["lineages"]).open())}
    drugs = list(org["drugs"])
    labels = load_labels(org_key, org, drugs)
    outdir = MODELS / org_key
    outdir.mkdir(parents=True, exist_ok=True)

    saved = []
    print(f"\n=== {org['display']} ({org_key}) ===")
    for drug in drugs:
        lab = labels[drug]
        g = [x for x in feats.index if lab.get(x) in ("Resistant", "Susceptible")]
        y = np.array([1 if lab[x] == "Resistant" else 0 for x in g])
        if y.sum() < MIN_PER_CLASS or (y == 0).sum() < MIN_PER_CLASS:
            print(f"  {drug:32s} SKIP — too few per class (R={int(y.sum())} S={int((y==0).sum())})")
            continue
        X = feats.loc[g].values.astype(int)
        groups = np.array([lineage.get(x, x) for x in g])

        spw = (y == 0).sum() / max((y == 1).sum(), 1)
        base = XGBClassifier(n_estimators=300, max_depth=4, learning_rate=0.1, subsample=0.9,
                             colsample_bytree=0.8, eval_metric="logloss", scale_pos_weight=spw,
                             random_state=seed, n_jobs=4)
        cv = StratifiedGroupKFold(n_splits=4, shuffle=True, random_state=seed)
        oof = cross_val_predict(base, X, y, groups=groups, cv=cv, method="predict_proba")[:, 1]
        iso = IsotonicRegression(out_of_bounds="clip").fit(oof, y)
        thr = pick_threshold(y, iso.transform(oof), vme_target)
        base.fit(X, y)   # deployment model on all data

        rules = rules_for(org_key, drug)
        joblib.dump({"model": base, "calibrator": iso, "threshold": float(thr),
                     "feature_cols": cols, "drug": drug, "rules": rules,
                     "organism": org_key, "amrfinder": org["amrfinder"],
                     "display": org["display"], "n_train": int(len(y)),
                     "n_resistant": int(y.sum()), "vme_target": vme_target},
                    outdir / f"{drug}.joblib")
        saved.append(drug)
        print(f"  {drug:32s} n={len(y):>4} R={int(y.sum()):>4} thr={thr:.2f} -> {org_key}/{drug}.joblib")
    return saved


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/config.yaml")
    ap.add_argument("--organism", choices=ORG_ORDER, help="default: all eight")
    args = ap.parse_args()
    cfg = yaml.safe_load((REPO / args.config).read_text())
    seed, vme_target = cfg["project"]["seed"], cfg["evaluation"]["vme_target"]

    orgs = [args.organism] if args.organism else ORG_ORDER
    total = 0
    for ok in orgs:
        total += len(save_organism(ok, seed, vme_target))
    print(f"\nsaved {total} deployable models across {len(orgs)} organism(s) -> results/models/<org>/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
