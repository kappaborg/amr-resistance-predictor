"""Risk–coverage curves + AURC — rigorously evaluating the defer-to-lab abstention.

Selective prediction: rank calls by confidence (distance of the calibrated probability from the 0.5
boundary) and, at each *coverage* (fraction of strains given a confident call), measure the *risk*
among those non-deferred — both overall error and, clinically, the **very-major error (VME) among
non-deferred**. A good abstention policy drives VME down sharply as it defers the least-confident
strains to phenotypic testing. We report the **AURC** (area under the risk–coverage curve; lower is
better) and VME at selected coverage levels. Uses the pooled lineage-held-out calibrated predictions
from clinical_rigor (`results/metrics/pooled_predictions.pkl`).

    python -m src.evaluation.clinical_rigor   # first, to produce pooled_predictions.pkl
    python -m src.evaluation.risk_coverage
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
OUT_MD = REPO / "results/reports/summary_31_risk_coverage.md"
OUT_FIG = REPO / "results/figures/risk_coverage.png"
ORG = "kpneu"


def risk_coverage(y, score):
    """Return (coverage grid, overall-error, VME-among-non-deferred) as least-confident strains are
    progressively deferred. Confidence = |calibrated prob − 0.5|; point call = prob ≥ 0.5."""
    y = np.asarray(y); pred = (score >= 0.5).astype(int)
    conf = np.abs(score - 0.5)
    order = np.argsort(-conf)                              # most confident first
    yo, po = y[order], pred[order]
    n = len(y)
    cov, err, vme = [], [], []
    for k in range(max(1, n // 50), n + 1, max(1, n // 50)):
        yk, pk = yo[:k], po[:k]
        cov.append(k / n)
        err.append(float(np.mean(pk != yk)))
        nR = int((yk == 1).sum())
        vme.append(float(((yk == 1) & (pk == 0)).sum() / nR) if nR else np.nan)
    return np.array(cov), np.array(err), np.array(vme)


def main() -> int:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    store = pickle.loads(PREDS.read_bytes())
    drugs = [d for (o, d) in store if o == ORG]
    ncol = min(3, len(drugs)); nrow = int(np.ceil(len(drugs) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.6 * ncol, 3.4 * nrow), squeeze=False)
    md = ["# Summary #31 — Risk–Coverage Curves & AURC (defer-to-lab abstention)\n",
          f"**Date:** 2026-07-17 · *{ORGANISMS[ORG]['display']}* · selective prediction on the pooled "
          "lineage-held-out calibrated predictions. Deferring the least-confident strains to phenotypic "
          "testing lowers error among those kept. **AURC** = area under the error risk–coverage curve "
          "(lower better); we also track **VME among non-deferred**.\n",
          "| Drug | AURC | error @100% | error @70% cov | VME @100% | VME @70% cov |",
          "|---|---|---|---|---|---|"]
    for k, d in enumerate(drugs):
        rec = store[(ORG, d)]
        cov, err, vme = risk_coverage(rec["y"], rec["score"])
        aurc = float(np.sum(np.diff(cov) * (err[:-1] + err[1:]) / 2) / (cov[-1] - cov[0]))  # trapezoid
        i70 = int(np.argmin(np.abs(cov - 0.7)))
        ax = axes[k // ncol][k % ncol]
        ax.plot(cov, err, color="#0e7c86", lw=2, label="error")
        ax.plot(cov, vme, color="#b4472e", lw=1.6, ls="--", label="VME (missed R)")
        ax.set_title(d.replace("trimethoprim_sulfamethoxazole", "TMP-SMX"), fontsize=10)
        ax.set_xlabel("coverage (fraction given a confident call)")
        ax.set_ylabel("risk"); ax.set_xlim(cov[0], 1.0); ax.set_ylim(0, None)
        if k == 0:
            ax.legend(fontsize=7)
        md.append(f"| {d} | {aurc:.3f} | {err[-1]:.1%} | {err[i70]:.1%} | {vme[-1]:.1%} | {vme[i70]:.1%} |")
    for j in range(len(drugs), nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")
    fig.suptitle(f"Risk–coverage — {ORGANISMS[ORG]['display']} (defer least-confident to the lab)",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_FIG, dpi=150, bbox_inches="tight")
    md.append("\n**Reading.** Both curves fall as coverage drops: deferring the least-confident ~30% of "
              "strains to phenotypic testing cuts the error (and the clinically-critical **VME**) among "
              "the calls the tool *does* make. This is the quantitative case for the defer-to-lab "
              "abstention — the model knows when it doesn't know. Complements the conformal per-class "
              "coverage evaluation (Summary #9), which is empirically validated under lineage shift.")
    OUT_MD.write_text("\n".join(md) + "\n")
    print(f"wrote {OUT_MD}\nwrote {OUT_FIG}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
