"""The 8-organism generalization figure — per-organism × drug ROC-AUC heatmap.

Reads the authoritative calibrated evaluation (`results/metrics/clinical_rigor.json`, produced by
`clinical_rigor.py`: isotonic-calibrated XGBoost, lineage-held-out, raw-score AUC) and renders a
heatmap of unseen-lineage ROC-AUC for every organism × drug — the visual support for the project's
"8 WHO-priority pathogens" generalization claim. Cells annotated; drug not modelled in an organism is
left blank. Colour scale 0.5 (chance) → 1.0.

    python -m src.evaluation.figure_generalization
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from src.app.registry import ORGANISMS, ORG_ORDER  # noqa: E402

RIGOR = REPO / "results/metrics/clinical_rigor.json"
OUT = REPO / "results/figures/generalization_heatmap.png"

# column order: group drugs roughly by class for readability
DRUG_ORDER = ["meropenem", "imipenem", "ceftazidime", "ceftriaxone", "cefoxitin", "ampicillin",
              "penicillin", "oxacillin", "ciprofloxacin", "gentamicin", "tobramycin", "amikacin",
              "vancomycin", "erythromycin", "clindamycin", "tetracycline", "chloramphenicol",
              "trimethoprim_sulfamethoxazole"]
DRUG_SHORT = {"trimethoprim_sulfamethoxazole": "TMP-SMX"}


def main() -> int:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    data = json.loads(RIGOR.read_text())
    orgs = [o for o in ORG_ORDER if data.get(o)]
    drugs = [d for d in DRUG_ORDER if any(d in data[o] for o in orgs)]
    M = np.full((len(orgs), len(drugs)), np.nan)
    for i, o in enumerate(orgs):
        for j, d in enumerate(drugs):
            if d in data[o]:
                M[i, j] = data[o][d]["auc_model"]

    fig, ax = plt.subplots(figsize=(0.62 * len(drugs) + 3, 0.55 * len(orgs) + 2))
    cmap = plt.cm.RdYlGn.copy()
    cmap.set_bad("#e8e8e8")
    im = ax.imshow(np.ma.masked_invalid(M), cmap=cmap, vmin=0.5, vmax=1.0, aspect="auto")
    ax.set_xticks(range(len(drugs)))
    ax.set_xticklabels([DRUG_SHORT.get(d, d) for d in drugs], rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(orgs)))
    ax.set_yticklabels([ORGANISMS[o]["display"] for o in orgs], fontsize=9, fontstyle="italic")
    for i in range(len(orgs)):
        for j in range(len(drugs)):
            if not np.isnan(M[i, j]):
                ax.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center", fontsize=7,
                        color="black" if M[i, j] > 0.72 else "white")
    ax.set_title("Unseen-lineage ROC-AUC across 8 WHO-priority pathogens\n"
                 "(calibrated XGBoost, lineage-held-out; grey = drug not in that organism's panel)",
                 fontsize=10)
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cb.set_label("ROC-AUC (0.5 = chance)", fontsize=8)
    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=150, bbox_inches="tight")
    print(f"wrote {OUT}  ({len(orgs)} organisms × {len(drugs)} drugs, "
          f"{int(np.isfinite(M).sum())} models)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
