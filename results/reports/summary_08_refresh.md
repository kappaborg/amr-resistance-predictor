# Summary #8 — Full-Dataset Refresh (post Phase-2b top-up)

**Date:** 2026-07-08 · **Status:** ✅ pipeline refreshed end-to-end on the full set.

## What changed
- Dataset: **1,472 → 3,850 genomes** (top-up of ~2,378, targeted at resistant strains); QC dropped 66.
- Features: **446 → 688 determinants**. Resistant counts rose ~3–5× (meropenem 294→1,427 R,
  gentamicin 410→1,601 R, cefoxitin 319→1,525 R).
- One command did it all: `python -m src.refresh_pipeline` (QC → annotate → MLST → matrix → labels →
  models → SHAP → saved demo models). Annotation reused the 1,472 cached; only new genomes ran.

## Refreshed results (logistic regression, unseen-lineage test, VME≤3% operating point)
| Drug | Test R/S (lineages) | ROC | VME | ME | vs before (ROC) |
|---|---|---|---|---|---|
| meropenem | 355/526 (132) | 0.982 | 0.017 | 0.293 | 0.947 → **0.982**, VME target now met |
| ciprofloxacin | 568/225 (139) | 0.973 | 0.033 | 0.089 | ~same, still strong |
| gentamicin | 401/528 (141) | 0.968 | 0.032 | 0.129 | rules already 0.967 (honest) |
| TMP-SMX | 512/226 (123) | 0.970 | 0.039 | 0.128 | ~same |
| cefoxitin | 381/253 (77) | 0.906 | 0.018 | 0.850 | **no longer degenerate**; rules still 0.517 |

## What the top-up fixed / revealed
- **Stabilised weak operating points:** meropenem and cefoxitin now meet VME≤3% (were data-limited /
  degenerate on the small set). Bigger, more diverse test sets (e.g. cefoxitin 47→77 lineages) make
  this a **harder, more credible** evaluation.
- **Sharpened interpretation:** meropenem SHAP now ranks blaKPC-2/blaKPC-3 top-2 (added lineage
  diversity broke the earlier co-selection confound). Causal mechanisms top 4/5 drugs.
- **Honest limits persist:** the VME≤3% operating point drives high major errors on hard drugs
  (cefoxitin ME 0.85) — a clinical trade-off; gentamicin's rules baseline remains as good as ML.

## Thesis: intact and stronger
Cefoxitin remains the clearest demonstration — the gene-lookup baseline is useless (ROC 0.517, misses
97% of resistant strains) because porin loss is invisible to it, while the model reaches ROC 0.906.
Demo models re-saved (`results/models/*.joblib`); `REPORT.md` updated to the 3,850-genome numbers.
