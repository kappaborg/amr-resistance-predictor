"""Conformal prediction — statistically-guaranteed uncertainty + a clinical abstain ("defer to
phenotypic testing") option, per drug.

Method: class-conditional (Mondrian) inductive conformal prediction. For each class c we take the
(1-alpha) quantile of the nonconformity scores (1 - calibrated P(true class)) on a lineage-disjoint
CALIBRATION set, then a test point's prediction set includes label c iff its nonconformity <= q_c.
This gives a per-class coverage guarantee: P(Y in set | Y = c) >= 1 - alpha. Setting a small alpha
for the RESISTANT class bounds the very-major error rate WITH A STATISTICAL GUARANTEE.

Prediction sets: {R} or {S} = a confident call; {R,S} or {} = uncertain -> defer to lab testing.

Everything is lineage-disjoint: train-proper / calibration / test share no MLST sequence type.

Usage:
    python -m src.models.conformal --config config/config.yaml --alpha 0.05
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
import yaml  # noqa: E402
from sklearn.frozen import FrozenEstimator  # noqa: E402
from sklearn.calibration import CalibratedClassifierCV  # noqa: E402
from sklearn.model_selection import GroupShuffleSplit, StratifiedGroupKFold  # noqa: E402
from xgboost import XGBClassifier  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from src.models.train_panel import RULES  # noqa: E402

FEATURES = REPO / "data/processed/thin_slice_cipro_features.csv"
LINEAGES = REPO / "data/processed/thin_slice_cipro_lineages.csv"
PANEL = REPO / "data/processed/panel_labels.csv"
OUT_JSON = REPO / "results/metrics/conformal_metrics.json"
FIG = REPO / "results/figures/conformal.png"


def conformal_drug(drug, feats, lineage, labels, seed, alpha):
    sub = labels[labels[drug].isin(["Resistant", "Susceptible"])]
    sub = sub[sub["genome_id"].isin(feats.index)]
    g = sub["genome_id"].values
    X = feats.loc[g].values.astype(int)
    y = (sub[drug].values == "Resistant").astype(int)
    groups = np.array([lineage[x] for x in g])

    sgkf = StratifiedGroupKFold(n_splits=4, shuffle=True, random_state=seed)
    tr, te = next(sgkf.split(X, y, groups))
    # split train into proper-train / calibration, lineage-disjoint
    fit, cal = next(GroupShuffleSplit(n_splits=1, test_size=0.35, random_state=seed)
                    .split(X[tr], y[tr], groups[tr]))
    Xfit, yfit = X[tr][fit], y[tr][fit]
    Xcal, ycal = X[tr][cal], y[tr][cal]
    Xte, yte = X[te], y[te]

    spw = (yfit == 0).sum() / max((yfit == 1).sum(), 1)
    base = XGBClassifier(n_estimators=300, max_depth=4, learning_rate=0.1, subsample=0.9,
                         colsample_bytree=0.8, eval_metric="logloss", scale_pos_weight=spw,
                         random_state=seed, n_jobs=4).fit(Xfit, yfit)
    cal_model = CalibratedClassifierCV(FrozenEstimator(base), method="isotonic").fit(Xcal, ycal)

    def pR(X_):
        return cal_model.predict_proba(X_)[:, 1]

    # class-conditional nonconformity quantiles on the calibration set
    p_cal = pR(Xcal)
    # nonconformity for the TRUE label: 1 - P(true class)
    nc_R = 1 - p_cal[ycal == 1]        # scores for truly-Resistant calibration points
    nc_S = 1 - (1 - p_cal[ycal == 0])  # = p_cal for truly-Susceptible points
    # conformal quantile with finite-sample correction
    def qhat(scores):
        n = len(scores)
        k = min(n, int(np.ceil((n + 1) * (1 - alpha))))
        return np.sort(scores)[k - 1] if n else 1.0
    qR, qS = qhat(nc_R), qhat(nc_S)

    # test prediction sets
    p_te = pR(Xte)
    inc_R = (1 - p_te) <= qR          # include "Resistant" in the set?
    inc_S = p_te <= qS                # include "Susceptible"?
    setsize = inc_R.astype(int) + inc_S.astype(int)
    singleton = setsize == 1
    confident_R = inc_R & ~inc_S
    confident_S = inc_S & ~inc_R

    # empirical per-class coverage (should be >= 1-alpha)
    cov_R = float((inc_R[yte == 1]).mean()) if (yte == 1).any() else float("nan")
    cov_S = float((inc_S[yte == 0]).mean()) if (yte == 0).any() else float("nan")
    # on CONFIDENT (singleton) calls: error rates
    conf_mask = singleton
    n_conf = int(conf_mask.sum())
    pred_conf = confident_R[conf_mask].astype(int)  # 1 if called R
    yconf = yte[conf_mask]
    vme_conf = float(((yconf == 1) & (pred_conf == 0)).sum() / max((yconf == 1).sum(), 1))
    me_conf = float(((yconf == 0) & (pred_conf == 1)).sum() / max((yconf == 0).sum(), 1))
    return {
        "alpha": alpha, "n_test": int(len(yte)),
        "coverage_R": cov_R, "coverage_S": cov_S,          # guarantee check: >= 1-alpha
        "confident_rate": float(singleton.mean()),          # fraction with a decisive call
        "defer_rate": float((~singleton).mean()),           # {R,S} or {} -> send to lab
        "n_confident": n_conf,
        "vme_confident": vme_conf, "me_confident": me_conf,  # errors among confident calls
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/config.yaml")
    ap.add_argument("--alpha", type=float, default=0.05, help="target per-class error (VME bound)")
    args = ap.parse_args()
    cfg = yaml.safe_load((REPO / args.config).read_text())
    seed = cfg["project"]["seed"]

    feats = pd.read_csv(FEATURES, dtype={"genome_id": str}).set_index("genome_id")
    lineage = {r["genome_id"]: r["lineage"] for r in csv.DictReader(LINEAGES.open())}
    labels = pd.read_csv(PANEL, dtype=str)

    results = {}
    print(f"Conformal prediction (class-conditional, alpha={args.alpha} → guaranteed ≥{1-args.alpha:.0%} "
          f"per-class coverage)\n")
    print(f"{'drug':30s}{'covR':>6}{'covS':>6}{'confident%':>11}{'defer%':>8}"
          f"{'VME_conf':>9}{'ME_conf':>8}")
    for drug in RULES:
        r = conformal_drug(drug, feats, lineage, labels, seed, args.alpha)
        results[drug] = r
        print(f"{drug:30s}{r['coverage_R']:6.2f}{r['coverage_S']:6.2f}"
              f"{r['confident_rate']*100:10.1f}%{r['defer_rate']*100:7.1f}%"
              f"{r['vme_confident']:9.3f}{r['me_confident']:8.3f}")

    OUT_JSON.write_text(json.dumps(results, indent=2))

    # figure: confident vs deferred, and VME on confident calls
    drugs = list(results)
    names = [d.replace("trimethoprim_sulfamethoxazole", "TMP-SMX") for d in drugs]
    conf = [results[d]["confident_rate"] for d in drugs]
    defer = [results[d]["defer_rate"] for d in drugs]
    vme = [results[d]["vme_confident"] for d in drugs]
    x = np.arange(len(drugs))
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5))
    a1.bar(x, conf, color="#2c7fb8", label="confident call")
    a1.bar(x, defer, bottom=conf, color="#fdae61", label="defer to lab testing")
    a1.set_xticks(x); a1.set_xticklabels(names, rotation=25, ha="right"); a1.set_ylabel("fraction of strains")
    a1.set_title(f"Selective prediction (α={args.alpha})"); a1.legend(fontsize=8)
    a2.bar(x, vme, color="#d7301f"); a2.axhline(args.alpha, ls="--", color="k", lw=1, label=f"α={args.alpha}")
    a2.set_xticks(x); a2.set_xticklabels(names, rotation=25, ha="right")
    a2.set_ylabel("VME among confident calls"); a2.set_title("Guaranteed low miss-rate on confident calls")
    a2.legend(fontsize=8)
    fig.suptitle("Conformal prediction — guaranteed uncertainty + abstention (K. pneumoniae)", fontweight="bold")
    fig.tight_layout(); fig.savefig(FIG, dpi=130); plt.close(fig)
    print(f"\nwrote {OUT_JSON}, {FIG}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
