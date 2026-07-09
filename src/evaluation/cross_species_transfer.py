"""Round-3 — cross-species transfer (zero-shot generalization across organisms).

Trains a per-drug model on ONE organism and tests it on ANOTHER, using the shared determinant
vocabulary, for the drugs both panels share (ciprofloxacin, gentamicin, TMP-SMX). Answers: does the
determinant→phenotype relationship transfer across species without any target-organism training?
This is the multi-species-transfer question the field flags as the frontier — and we can do it
directly from the features we already have.

Usage:
    python -m src.evaluation.cross_species_transfer --config config/config.yaml
"""
from __future__ import annotations

import argparse
import csv
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from src.organism_pipeline import ORGANISMS  # noqa: E402

AMR_API = "https://www.bv-brc.org/api/genome_amr/"
LAB = 'eq(evidence,%22Laboratory%20Method%22)'


def fetch_labels(taxon: int, drug_enc: str) -> dict[str, str]:
    """Robust: retry, verify the response is a list, keep consistent-label genomes only."""
    q = (f"and(eq(taxon_id,{taxon}),{LAB},eq(antibiotic,%22{drug_enc}%22))"
         f"&select(genome_id,resistant_phenotype)&limit(25000)")
    for attempt in range(4):
        try:
            r = requests.get(f"{AMR_API}?{q}", headers={"Accept": "application/json"}, timeout=180)
            data = r.json()
            if r.status_code == 200 and isinstance(data, list):
                byg: dict[str, set] = defaultdict(set)
                for x in data:
                    p = x.get("resistant_phenotype")
                    if p in ("Resistant", "Susceptible"):
                        byg[x["genome_id"]].add(p)
                return {g: next(iter(v)) for g, v in byg.items() if len(v) == 1}
        except Exception:  # noqa: BLE001
            pass
        time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"label fetch failed for taxon {taxon} / {drug_enc}")

SHARED = ["ciprofloxacin", "gentamicin", "trimethoprim_sulfamethoxazole"]
KP_FEATURES = REPO / "data/processed/thin_slice_cipro_features.csv"
KP_PANEL = REPO / "data/processed/panel_labels.csv"
EC_FEATURES = REPO / "data/processed/ecoli_features.csv"
OUT = REPO / "results/reports/summary_16_cross_species.md"


def kp_labels_from_file(path: Path) -> dict[str, dict[str, str]]:
    """Reuse the saved K. pneumoniae panel labels (no re-fetch)."""
    rows = list(csv.DictReader(path.open()))
    return {d: {r["genome_id"]: r[d] for r in rows if r.get(d) in ("Resistant", "Susceptible")}
            for d in SHARED}


def load_org(feat_path, taxon, drugs_enc):
    feats = pd.read_csv(feat_path, dtype={"genome_id": str}).set_index("genome_id")
    labs = {col: fetch_labels(taxon, enc) for col, enc in drugs_enc.items()}
    return feats, labs


def xy(feats, labs, drug, cols):
    g = [x for x in feats.index if labs[drug].get(x) in ("Resistant", "Susceptible")]
    X = feats.loc[g, cols].values.astype(int)
    y = np.array([1 if labs[drug][x] == "Resistant" else 0 for x in g])
    return X, y


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/config.yaml")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    kp_feats = pd.read_csv(KP_FEATURES, dtype={"genome_id": str}).set_index("genome_id")
    kp_labs = kp_labels_from_file(KP_PANEL)
    ec_feats, ec_labs = load_org(EC_FEATURES, 562, {d: ORGANISMS["ecoli"]["drugs"][d] for d in SHARED})
    common = [c for c in kp_feats.columns if c in set(ec_feats.columns)]
    print(f"shared determinant features: {len(common)} "
          f"(K.pneu {kp_feats.shape[1]}, E.coli {ec_feats.shape[1]})")

    def fit_eval(Xtr, ytr, Xte, yte):
        m = LogisticRegression(max_iter=2000, class_weight="balanced",
                               random_state=args.seed).fit(Xtr, ytr)
        return roc_auc_score(yte, m.predict_proba(Xte)[:, 1])

    md = ["# Summary #16 — Cross-Species Transfer (zero-shot across organisms)\n",
          f"Train a per-drug model on one organism, test on the other, over the **{len(common)} shared "
          "determinants**. Zero target-organism training. ROC-AUC.\n",
          "| Drug | K.pneu→E.coli | E.coli→K.pneu | (within K.pneu) | (within E.coli) |",
          "|---|---|---|---|---|"]
    print(f"\n{'drug':30s}{'KP->EC':>9}{'EC->KP':>9}{'KP(in)':>9}{'EC(in)':>9}")
    for drug in SHARED:
        Xk, yk = xy(kp_feats, kp_labs, drug, common)
        Xe, ye = xy(ec_feats, ec_labs, drug, common)
        # 70/30 within-organism reference splits (random; transfer is the headline)
        rng = np.random.RandomState(args.seed)
        def split(X, y):
            idx = rng.permutation(len(y)); c = int(0.7 * len(y))
            return X[idx[:c]], y[idx[:c]], X[idx[c:]], y[idx[c:]]
        kp_tr_X, kp_tr_y, kp_te_X, kp_te_y = split(Xk, yk)
        ec_tr_X, ec_tr_y, ec_te_X, ec_te_y = split(Xe, ye)
        kp2ec = fit_eval(Xk, yk, Xe, ye)             # train ALL K.pneu, test ALL E.coli
        ec2kp = fit_eval(Xe, ye, Xk, yk)
        kp_in = fit_eval(kp_tr_X, kp_tr_y, kp_te_X, kp_te_y)
        ec_in = fit_eval(ec_tr_X, ec_tr_y, ec_te_X, ec_te_y)
        md.append(f"| {drug} | **{kp2ec:.3f}** | **{ec2kp:.3f}** | {kp_in:.3f} | {ec_in:.3f} |")
        print(f"{drug:30s}{kp2ec:>9.3f}{ec2kp:>9.3f}{kp_in:>9.3f}{ec_in:>9.3f}")

    md.append("\n**Reading:** a cross-species ROC-AUC well above 0.5 means the determinant→phenotype "
              "mapping the model learned on one organism **transfers to a different species without "
              "retraining** — strong evidence it captures real mechanism, not organism-specific "
              "lineage structure. Gaps vs the within-organism reference quantify the species shift.")
    OUT.write_text("\n".join(md) + "\n")
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
