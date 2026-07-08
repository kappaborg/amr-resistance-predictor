"""Phase 10 (Week 3) — global SHAP per drug + biological validation.

For each drug: train XGBoost under the lineage-aware split, compute SHAP values on the held-out
(unseen-lineage) test set, rank determinants by mean |SHAP|, and check the top features against
KNOWN resistance mechanisms. Agreement with established biology is evidence the model learned real
signal (the proposal's biological-validation step, to be reviewed by the microbiologist).

Outputs: results/figures/shap_<drug>.png (beeswarm), results/metrics/shap_top.json,
results/reports/summary_06_shap.md.

Usage:
    python -m src.interpret.shap_analysis --config config/config.yaml
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import shap  # noqa: E402
import yaml  # noqa: E402
from sklearn.model_selection import StratifiedGroupKFold  # noqa: E402
from xgboost import XGBClassifier  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
FEATURES = REPO / "data/processed/thin_slice_cipro_features.csv"
LINEAGES = REPO / "data/processed/thin_slice_cipro_lineages.csv"
PANEL = REPO / "data/processed/panel_labels.csv"
FIGS = REPO / "results/figures"
DRUGS = ["meropenem", "gentamicin", "ciprofloxacin", "trimethoprim_sulfamethoxazole", "cefoxitin"]

# Known resistance mechanisms per drug: token -> human-readable mechanism (for biological validation).
KNOWN = {
    "meropenem": {"blaKPC": "KPC carbapenemase", "blaNDM": "NDM metallo-carbapenemase",
                  "blaOXA-48": "OXA-48 carbapenemase", "blaOXA-181": "OXA-48-like",
                  "blaOXA-232": "OXA-48-like", "blaVIM": "VIM carbapenemase",
                  "blaIMP": "IMP carbapenemase", "ompK35": "porin loss", "ompK36": "porin loss"},
    "gentamicin": {"aac(3)": "aminoglycoside acetyltransferase", "ant(2'')": "nucleotidyltransferase",
                   "armA": "16S rRNA methyltransferase", "rmt": "16S rRNA methyltransferase"},
    "ciprofloxacin": {"gyrA": "gyrase QRDR mutation", "parC": "topoisomerase QRDR mutation",
                      "gyrB": "gyrase QRDR", "parE": "topoisomerase QRDR", "qnr": "PMQR (Qnr)",
                      "aac(6')-Ib-cr": "PMQR (aac-cr)", "oqxA": "OqxAB efflux",
                      "oqxB": "OqxAB efflux", "qepA": "QepA efflux"},
    "trimethoprim_sulfamethoxazole": {"sul1": "sulfonamide resistance", "sul2": "sulfonamide",
                                      "sul3": "sulfonamide", "dfrA": "dihydrofolate reductase",
                                      "dfrB": "dihydrofolate reductase"},
    "cefoxitin": {"blaCMY": "CMY AmpC", "blaDHA": "DHA AmpC", "blaACT": "ACT AmpC",
                  "blaFOX": "FOX AmpC", "blaMOX": "MOX AmpC", "ompK35": "porin loss",
                  "ompK36": "porin loss"},
}


def known_mechanism(feature: str, drug: str) -> str | None:
    for token, mech in KNOWN[drug].items():
        if token in feature:
            return mech
    return None


def run_drug(drug, feats, lineage, labels, seed, test_frac):
    sub = labels[labels[drug].isin(["Resistant", "Susceptible"])]
    sub = sub[sub["genome_id"].isin(feats.index)]
    g = sub["genome_id"].values
    X = feats.loc[g].values.astype(int)
    y = (sub[drug].values == "Resistant").astype(int)
    groups = np.array([lineage[x] for x in g])
    cols = list(feats.columns)

    sgkf = StratifiedGroupKFold(n_splits=max(2, round(1 / test_frac)), shuffle=True, random_state=seed)
    tr, te = next(sgkf.split(X, y, groups))
    spw = (y[tr] == 0).sum() / max((y[tr] == 1).sum(), 1)
    model = XGBClassifier(n_estimators=300, max_depth=4, learning_rate=0.1, subsample=0.9,
                          colsample_bytree=0.8, eval_metric="logloss", scale_pos_weight=spw,
                          random_state=seed, n_jobs=4).fit(X[tr], y[tr])

    expl = shap.TreeExplainer(model)
    sv = expl.shap_values(X[te])  # (n_test, n_features)
    mean_abs = np.abs(sv).mean(0)
    # signed direction: mean SHAP among genomes where the feature is present (does presence push to R?)
    order = np.argsort(mean_abs)[::-1]
    top = []
    for i in order[:15]:
        present = X[te][:, i] == 1
        direction = float(sv[present, i].mean()) if present.any() else 0.0
        top.append({"feature": cols[i], "mean_abs_shap": float(mean_abs[i]),
                    "pushes_toward": "Resistant" if direction > 0 else "Susceptible",
                    "known_mechanism": known_mechanism(cols[i], drug)})

    # beeswarm plot
    FIGS.mkdir(parents=True, exist_ok=True)
    plt.figure()
    shap.summary_plot(sv, X[te], feature_names=cols, max_display=12, show=False)
    plt.title(f"SHAP — {drug} (test: unseen lineages)")
    plt.tight_layout()
    plt.savefig(FIGS / f"shap_{drug}.png", dpi=120)
    plt.close()
    return {"n_test": int(len(te)), "top": top}


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
    md = ["# Summary #6 — Week 3: SHAP Interpretation + Biological Validation\n",
          "Global SHAP (XGBoost, TreeExplainer) per drug on the held-out **unseen-lineage** test set. "
          "Top determinants checked against known mechanisms — ✓ = matches established biology "
          "(⚕ microbiology to confirm). Beeswarm plots in `results/figures/shap_<drug>.png`.\n"]
    for drug in DRUGS:
        print(f"\n=== {drug} ===")
        r = run_drug(drug, feats, lineage, labels, seed, test_frac)
        results[drug] = r
        n_known = sum(1 for t in r["top"][:10] if t["known_mechanism"])
        md.append(f"\n## {drug}  ({n_known}/10 top features match known mechanisms)\n")
        md.append("| Rank | Determinant | mean\\|SHAP\\| | pushes toward | known mechanism |")
        md.append("|---|---|---|---|---|")
        for rank, t in enumerate(r["top"][:10], 1):
            mech = f"✓ {t['known_mechanism']}" if t["known_mechanism"] else "— (investigate)"
            md.append(f"| {rank} | {t['feature']} | {t['mean_abs_shap']:.3f} | "
                      f"{t['pushes_toward']} | {mech} |")
        print(f"  top-10 known-mechanism agreement: {n_known}/10")
        for t in r["top"][:5]:
            print(f"    {t['feature']:22s} |SHAP|={t['mean_abs_shap']:.3f} -> {t['pushes_toward']}"
                  f"  {'✓'+t['known_mechanism'] if t['known_mechanism'] else ''}")

    (REPO / "results/metrics/shap_top.json").write_text(json.dumps(results, indent=2))
    (REPO / "results/reports/summary_06_shap.md").write_text("\n".join(md) + "\n")
    print("\nwrote results/metrics/shap_top.json, results/reports/summary_06_shap.md, figures/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
