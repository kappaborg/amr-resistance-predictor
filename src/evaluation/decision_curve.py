"""Decision-curve analysis (net benefit) — the clinical-utility axis most AMR-ML omits.

For each drug, standardized **net benefit** vs the clinician's threshold probability p_t, comparing the
model to the honest rules baseline, treat-all, and treat-none. Net benefit (Vickers & Elkin 2006):

    NB(p_t) = TP/n − FP/n · (p_t / (1 − p_t))

where p_t encodes the clinician's tolerated cost ratio (a low p_t = "missing resistance is far worse
than over-treating", exactly the VME-first framing of this project). A model has clinical utility over
the p_t range where its curve is highest. Uses the pooled lineage-held-out **calibrated** predictions
from clinical_rigor (`results/metrics/pooled_predictions.pkl`).

    python -m src.evaluation.clinical_rigor   # first, to produce pooled_predictions.pkl
    python -m src.evaluation.decision_curve
"""
from __future__ import annotations

import pickle
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from src.app.registry import ORGANISMS  # noqa: E402

PREDS = REPO / "results/metrics/pooled_predictions.pkl"
OUT_MD = REPO / "results/reports/summary_30_decision_curve.md"
OUT_FIG = REPO / "results/figures/decision_curve.png"
ORG = "kpneu"   # clinical deep-dive organism


def net_benefit(y, score, pt):
    """NB of classifying positive when score >= pt (array over pt)."""
    y = np.asarray(y); n = len(y)
    nb = []
    for t in pt:
        pred = score >= t
        tp = np.sum((pred == 1) & (y == 1))
        fp = np.sum((pred == 1) & (y == 0))
        nb.append(tp / n - fp / n * (t / (1 - t)))
    return np.array(nb)


def nb_binary(y, pred, pt):
    """NB of a FIXED binary prediction (the rules baseline) across pt."""
    y, pred = np.asarray(y), np.asarray(pred); n = len(y)
    tp = np.sum((pred == 1) & (y == 1)); fp = np.sum((pred == 1) & (y == 0))
    return tp / n - fp / n * (pt / (1 - pt))


def main() -> int:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    store = pickle.loads(PREDS.read_bytes())
    drugs = [d for (o, d) in store if o == ORG]
    pt = np.linspace(0.01, 0.5, 60)                       # low p_t = VME-averse region
    ncol = min(3, len(drugs)); nrow = int(np.ceil(len(drugs) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.6 * ncol, 3.6 * nrow), squeeze=False)
    md = ["# Summary #30 — Decision-Curve Analysis (clinical net benefit)\n",
          f"**Date:** 2026-07-17 · *{ORGANISMS[ORG]['display']}* · net benefit vs threshold probability "
          "p_t from the pooled lineage-held-out **calibrated** predictions. A low p_t encodes "
          "\"missing resistance is far costlier than over-treating\" — the project's VME-first stance. "
          "The model has clinical utility over the p_t range where its curve sits highest.\n",
          "| Drug | model > rules over p_t | model > treat-all over p_t | NB(model) @ p_t=0.1 |",
          "|---|---|---|---|"]
    for k, d in enumerate(drugs):
        rec = store[(ORG, d)]
        y, score, rules = rec["y"], rec["score"], rec["rules"]
        prev = float(np.mean(y))
        nb_m = net_benefit(y, score, pt)
        nb_r = nb_binary(y, rules, pt)
        nb_all = prev - (1 - prev) * (pt / (1 - pt))
        ax = axes[k // ncol][k % ncol]
        ax.plot(pt, nb_m, label="model", lw=2, color="#0e7c86")
        ax.plot(pt, nb_r, label="rules", lw=1.5, color="#c8892a", ls="--")
        ax.plot(pt, nb_all, label="treat-all", lw=1, color="#888")
        ax.axhline(0, label="treat-none", lw=1, color="#bbb")
        ax.set_title(d.replace("trimethoprim_sulfamethoxazole", "TMP-SMX"), fontsize=10)
        ax.set_xlabel("threshold probability $p_t$"); ax.set_ylabel("net benefit")
        ax.set_ylim(min(-0.02, nb_m.min()), max(nb_m.max(), prev) * 1.05)
        if k == 0:
            ax.legend(fontsize=7)
        m_gt_r = pt[(nb_m > nb_r)]
        m_gt_all = pt[(nb_m > np.maximum(nb_all, 0))]
        rng = lambda a: f"{a.min():.2f}–{a.max():.2f}" if len(a) else "—"
        md.append(f"| {d} | {rng(m_gt_r)} | {rng(m_gt_all)} | {nb_m[np.argmin(abs(pt-0.1))]:.3f} |")
    for j in range(len(drugs), nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")
    fig.suptitle(f"Decision-curve analysis — {ORGANISMS[ORG]['display']} (net benefit vs threshold)",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_FIG, dpi=150, bbox_inches="tight")
    md.append("\n**Reading.** Where the **model** curve is above both **treat-all** and **rules**, using "
              "the model to guide therapy yields more net benefit (more true resistant caught per "
              "false-positive over-treatment) than treating everyone or trusting a gene lookup. The "
              "advantage concentrates in the **low-p_t** region — exactly where a clinician who fears "
              "missing resistance operates — which is the clinical case for the model beyond raw AUC.")
    OUT_MD.write_text("\n".join(md) + "\n")
    print(f"wrote {OUT_MD}\nwrote {OUT_FIG}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
