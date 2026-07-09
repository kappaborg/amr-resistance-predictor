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


def fetch_mic(drug_enc: str) -> dict[str, float]:
    """genome_id -> log2(MIC in mg/L). Parses censored values (<=, >, =); keeps median per genome."""
    q = (f"and(eq(taxon_id,573),{LAB},eq(antibiotic,%22{drug_enc}%22))"
         f"&select(genome_id,measurement,measurement_unit)&limit(25000)")
    for attempt in range(4):
        r = requests.get(f"{AMR}?{q}", headers={"Accept": "application/json"}, timeout=180)
        if r.status_code == 200 and isinstance(r.json(), list):
            break
    byg: dict[str, list] = defaultdict(list)
    for x in r.json():
        m, unit = x.get("measurement"), x.get("measurement_unit")
        if not m or unit != "mg/L":
            continue
        num = re.sub(r"[<>=]", "", str(m)).strip()   # strip censor signs
        try:
            v = float(num)
        except ValueError:
            continue
        if v > 0:
            byg[x["genome_id"]].append(np.log2(v))
    return {g: float(np.median(vs)) for g, vs in byg.items() if vs}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/config.yaml")
    args = ap.parse_args()
    cfg = yaml.safe_load((REPO / args.config).read_text())
    seed, test_frac = cfg["project"]["seed"], cfg["split"]["test_fraction"]

    feats = pd.read_csv(FEATURES, dtype={"genome_id": str}).set_index("genome_id")
    lineage = {r["genome_id"]: r["lineage"] for r in csv.DictReader(LINEAGES.open())}

    md = ["# Summary #19 — MIC Regression (predicting the resistance *level*)\n",
          "Per-drug regression of continuous **log2(MIC, mg/L)** from determinant features, under the "
          "phylogeny-aware split. **Essential Agreement (EA)** = predictions within ±1 two-fold "
          "dilution (the CLSI/FDA metric for genotype-based MIC prediction); higher is better. RMSE "
          "in doubling dilutions; Pearson r.\n",
          "| Drug | n (MIC) | test | EA (±1 dilution) | RMSE (dilutions) | Pearson r |",
          "|---|---|---|---|---|---|"]
    results = {}
    print(f"{'drug':30s}{'n':>6}{'EA':>8}{'RMSE':>8}{'r':>7}")
    for drug, enc in DRUG_ENC.items():
        mic = fetch_mic(enc)
        g = [x for x in feats.index if x in mic]
        if len(g) < 100:
            print(f"{drug:30s} too few MIC values ({len(g)})")
            continue
        X = feats.loc[g].values.astype(int)
        y = np.array([mic[x] for x in g])
        groups = np.array([lineage[x] for x in g])
        # stratify the group split on binarised MIC (median) just to balance folds
        ybin = (y >= np.median(y)).astype(int)
        tr, te = next(StratifiedGroupKFold(max(2, round(1 / test_frac)), shuffle=True,
                                           random_state=seed).split(X, ybin, groups))
        model = XGBRegressor(n_estimators=400, max_depth=4, learning_rate=0.08, subsample=0.9,
                             colsample_bytree=0.8, random_state=seed, n_jobs=4).fit(X[tr], y[tr])
        pred = model.predict(X[te])
        ea = float(np.mean(np.abs(pred - y[te]) <= 1.0))          # within ±1 doubling dilution
        rmse = float(np.sqrt(np.mean((pred - y[te]) ** 2)))
        r = float(pearsonr(pred, y[te])[0]) if len(set(y[te])) > 1 else float("nan")
        results[drug] = {"n": len(g), "n_test": int(len(te)), "essential_agreement": ea,
                         "rmse_dilutions": rmse, "pearson_r": r}
        md.append(f"| {drug} | {len(g)} | {len(te)} | **{ea:.1%}** | {rmse:.2f} | {r:.3f} |")
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
