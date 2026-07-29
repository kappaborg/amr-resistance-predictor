# Poster

`poster.html` is a self-contained, print-ready academic poster (A0 landscape, 1189 × 841 mm).
It pulls the real figures from `../figures/` and the exact numbers from the evaluation.

## View it
```bash
open results/poster/poster.html          # macOS — opens in your default browser
```

## Export to PDF (for printing / submission)
1. Open `poster.html` in **Google Chrome** (best print fidelity).
2. **File → Print** (⌘P).
3. Set **Destination = Save as PDF**, **Layout = Landscape**, **Paper size = A0** (or the largest
   available; choose "Fit to page" if A0 isn't listed).
4. Enable **More settings → Background graphics** (so the colored headers/cards print).
5. Save.

To resize for a different board, change the `.poster{width; height}` values (and `@page size`) at the
top of `poster.html` — everything scales in millimetres.

## What's on it
Problem (two pitfalls) → approach → data/methods → per-drug results with **DeLong significance vs the
rules baseline** → the cefoxitin headline (0.52→0.90) → statistical rigor (leakage quantified, FDA/CLSI
CIs) → foundation models (ESM-2 adopted, GNN rejected) → cross-species transfer → SHAP interpretability
→ MIC + conformal abstention → 4-organism generalization + demo → honest limits & conclusion.

Figures referenced (must exist in `../figures/`): `benchmark_summary.png`, `cipro_transfer_matrix.png`,
`shap_meropenem.png`.
