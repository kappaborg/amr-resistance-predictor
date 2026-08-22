"""Phase 3/10 visuals — generate the headline result figures for the report/poster.

Produces:
  results/figures/benchmark_summary.png  — ROC-AUC per drug (ML vs rules baseline) + published band,
                                            and the VME/ME clinical trade-off.
  results/figures/roc_curves.png          — per-drug ROC curves (rules / logreg / XGBoost).
  results/figures/calibration.png         — per-drug reliability curves (calibrated XGBoost).

Recomputes test-set scores under the same lineage-aware split as the evaluation.

Usage:
    python -m src.evaluation.figures --config config/config.yaml
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import yaml  # noqa: E402
from sklearn.calibration import calibration_curve  # noqa: E402
from sklearn.frozen import FrozenEstimator  # noqa: E402
from sklearn.calibration import CalibratedClassifierCV  # noqa: E402
from sklearn.isotonic import IsotonicRegression  # noqa: E402
from sklearn.linear_model import LogisticRegression  # noqa: E402
from sklearn.metrics import roc_auc_score, roc_curve  # noqa: E402
from sklearn.model_selection import StratifiedGroupKFold, cross_val_predict  # noqa: E402
from xgboost import XGBClassifier  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from src.models.train_panel import RULES, rules_predict, pick_threshold  # noqa: E402

FEATURES = REPO / "data/processed/thin_slice_cipro_features.csv"
LINEAGES = REPO / "data/processed/thin_slice_cipro_lineages.csv"
PANEL = REPO / "data/processed/panel_labels.csv"
FIGS = REPO / "results/figures"
DRUGS = list(RULES)
# NOTE (decision #60): a "published band" of 0.85-0.96 used to be shaded on the ROC-AUC axis here.
# It was withdrawn: it carried no citation, and the published numbers it was loosely based on
# (Nguyen et al. 2018) are MIC essential-agreement rates, not ROC-AUC. Shading one metric's range
# on another metric's axis invited a false cross-study comparison. See REPORT.md section 4.1.


def drug_scores(drug, feats, lineage, labels, seed, test_frac):
    sub = labels[labels[drug].isin(["Resistant", "Susceptible"])]
    sub = sub[sub["genome_id"].isin(feats.index)]
    g = sub["genome_id"].values
    X = feats.loc[g].values.astype(int)
    y = (sub[drug].values == "Resistant").astype(int)
    groups = np.array([lineage[x] for x in g])
    cols = list(feats.columns)
    sgkf = StratifiedGroupKFold(n_splits=max(2, round(1 / test_frac)), shuffle=True, random_state=seed)
    tr, te = next(sgkf.split(X, y, groups))

    lr = LogisticRegression(max_iter=2000, class_weight="balanced", random_state=seed).fit(X[tr], y[tr])
    lr_s = lr.predict_proba(X[te])[:, 1]
    spw = (y[tr] == 0).sum() / max((y[tr] == 1).sum(), 1)
    base = XGBClassifier(n_estimators=300, max_depth=4, learning_rate=0.1, subsample=0.9,
                         colsample_bytree=0.8, eval_metric="logloss", scale_pos_weight=spw,
                         random_state=seed, n_jobs=4).fit(X[tr], y[tr])
    oof = cross_val_predict(base, X[tr], y[tr], groups=groups[tr],
                            cv=StratifiedGroupKFold(4, shuffle=True, random_state=seed),
                            method="predict_proba")[:, 1]
    iso = IsotonicRegression(out_of_bounds="clip").fit(oof, y[tr])
    xgb_s = iso.transform(base.predict_proba(X[te])[:, 1])
    rule_s = rules_predict(cols, X[te], RULES[drug]).astype(float)
    # Operating threshold from TRAIN out-of-fold predictions only. Previously the threshold was
    # picked on the test fold and then scored on that same fold, which made the plotted VME bars
    # satisfy the 3% target by construction (audit, decision #60).
    lr_oof = cross_val_predict(
        LogisticRegression(max_iter=2000, class_weight="balanced", random_state=seed),
        X[tr], y[tr], groups=groups[tr],
        cv=StratifiedGroupKFold(4, shuffle=True, random_state=seed),
        method="predict_proba")[:, 1]
    return {"y": y[te], "rules": rule_s, "logreg": lr_s, "xgb": xgb_s, "lr_oof": lr_oof, "y_tr": y[tr]}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/config.yaml")
    args = ap.parse_args()
    cfg = yaml.safe_load((REPO / args.config).read_text())
    seed, test_frac = cfg["project"]["seed"], cfg["split"]["test_fraction"]
    feats = pd.read_csv(FEATURES, dtype={"genome_id": str}).set_index("genome_id")
    lineage = {r["genome_id"]: r["lineage"] for r in csv.DictReader(LINEAGES.open())}
    labels = pd.read_csv(PANEL, dtype=str)
    FIGS.mkdir(parents=True, exist_ok=True)

    S = {d: drug_scores(d, feats, lineage, labels, seed, test_frac) for d in DRUGS}
    names = [d.replace("trimethoprim_sulfamethoxazole", "TMP-SMX") for d in DRUGS]

    # --- Figure 1: benchmark summary (ROC-AUC bars + clinical error) ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    x = np.arange(len(DRUGS))
    ml_auc = [roc_auc_score(S[d]["y"], S[d]["logreg"]) for d in DRUGS]
    rule_auc = [roc_auc_score(S[d]["y"], S[d]["rules"]) for d in DRUGS]
    ax1.bar(x - 0.2, ml_auc, 0.4, label="ML (logreg, lineage split)", color="#2c7fb8")
    ax1.bar(x + 0.2, rule_auc, 0.4, label="known-gene rules baseline", color="#c0c0c0")
    ax1.set_xticks(x); ax1.set_xticklabels(names, rotation=25, ha="right")
    ax1.set_ylabel("ROC-AUC"); ax1.set_ylim(0.4, 1.0)
    ax1.set_title("Discrimination on unseen lineages vs baseline"); ax1.legend(fontsize=8)
    for xi, a in zip(x - 0.2, ml_auc):
        ax1.text(xi, a + 0.01, f"{a:.2f}", ha="center", fontsize=8)

    # honest: threshold from held-out TRAIN predictions, never from the test fold
    thr = {d: pick_threshold(S[d]["y_tr"], S[d]["lr_oof"], cfg["evaluation"]["vme_target"])
           for d in DRUGS}
    vme = [((S[d]["y"] == 1) & (S[d]["logreg"] < thr[d])).sum() / max((S[d]["y"] == 1).sum(), 1)
           for d in DRUGS]
    me = [((S[d]["y"] == 0) & (S[d]["logreg"] >= thr[d])).sum() / max((S[d]["y"] == 0).sum(), 1)
          for d in DRUGS]
    ax2.bar(x - 0.2, vme, 0.4, label="very-major error (miss R)", color="#d7301f")
    ax2.bar(x + 0.2, me, 0.4, label="major error (false R)", color="#fdae61")
    ax2.axhline(cfg["evaluation"]["vme_target"], ls="--", color="k", lw=1, label="VME target 3%")
    ax2.set_xticks(x); ax2.set_xticklabels(names, rotation=25, ha="right")
    ax2.set_ylabel("error rate"); ax2.set_title("Clinical errors — threshold set on TRAIN, measured on unseen lineages")
    ax2.legend(fontsize=8)
    fig.suptitle("Reading Resistance — honest benchmark (K. pneumoniae, 3,850 genomes)", fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIGS / "benchmark_summary.png", dpi=130)
    plt.close(fig)

    # --- Figure 2: per-drug ROC curves ---
    fig, axes = plt.subplots(1, 5, figsize=(20, 4))
    for ax, d, nm in zip(axes, DRUGS, names):
        for key, lab, col in [("rules", "rules", "#999"), ("logreg", "logreg", "#2c7fb8"),
                              ("xgb", "XGBoost", "#31a354")]:
            fpr, tpr, _ = roc_curve(S[d]["y"], S[d][key])
            ax.plot(fpr, tpr, col, label=f"{lab} ({roc_auc_score(S[d]['y'], S[d][key]):.2f})")
        ax.plot([0, 1], [0, 1], "k:", lw=0.8)
        ax.set_title(nm); ax.set_xlabel("FPR"); ax.legend(fontsize=7, loc="lower right")
    axes[0].set_ylabel("TPR")
    fig.suptitle("Per-drug ROC — unseen-lineage test", fontweight="bold")
    fig.tight_layout(); fig.savefig(FIGS / "roc_curves.png", dpi=130); plt.close(fig)

    # --- Figure 3: calibration (reliability) curves ---
    fig, axes = plt.subplots(1, 5, figsize=(20, 4))
    for ax, d, nm in zip(axes, DRUGS, names):
        frac, mean = calibration_curve(S[d]["y"], S[d]["xgb"], n_bins=8, strategy="quantile")
        ax.plot(mean, frac, "o-", color="#31a354")
        ax.plot([0, 1], [0, 1], "k:", lw=0.8)
        ax.set_title(nm); ax.set_xlabel("predicted P(R)")
    axes[0].set_ylabel("observed fraction R")
    fig.suptitle("Calibration (isotonic XGBoost) — unseen-lineage test", fontweight="bold")
    fig.tight_layout(); fig.savefig(FIGS / "calibration.png", dpi=130); plt.close(fig)

    print("wrote results/figures/{benchmark_summary,roc_curves,calibration}.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
