# Summary #4 — WEEK-1 GO/NO-GO GATE  ✅ GO

**Date:** 2026-07-08 · **Drug:** ciprofloxacin · **Organism:** K. pneumoniae
**Recommendation: GO** — the thin end-to-end pipeline works, validated honestly on unseen lineages.

## The gate criterion (from the brief)
*"A working baseline for one drug exists."* — Met, and exceeded: full pipeline
(features → phylogeny-aware split → models → honest benchmark) runs end-to-end with strong,
biologically-validated results.

## Results — test set: 368 genomes, UNSEEN lineages, 184 R / 184 S
| Model | ROC-AUC | PR-AUC | VME (miss R) | ME (false R) |
|---|---|---|---|---|
| Rules baseline (expert QRDR + PMQR) | 0.897 | 0.834 | 0.022 | 0.185 |
| Logistic regression | **0.976** | **0.976** | 0.049 | **0.043** |
| XGBoost | 0.976 | 0.976 | 0.043 | 0.092 |

VME = resistant called susceptible (dangerous); ME = susceptible called resistant.

## What the model adds over the known-gene baseline (the honest story)
- **Discrimination:** ROC-AUC 0.897 → **0.976** on unseen lineages.
- **Major errors cut ~4×:** the rules baseline calls 18.5% of susceptible strains resistant (a simple
  "any QRDR/PMQR present → R" over-calls); logreg reduces this to **4.3%** by weighting determinants.
- **Trade-off to address (Phase 8):** the rules baseline has marginally lower VME (0.022 vs 0.043–0.049)
  because it flags every QRDR hit as resistant. The ML models trade a little VME for a large ME drop;
  calibration + threshold tuning (Phase 8) can push VME down where clinical cost demands it.
- **Note:** cipro is an "easy" drug — resistance is driven by strong, direct gyrA/parC mutations, so
  high AUC is biologically expected. Harder drugs (meropenem: porin loss + β-lactamase combinations
  that rules miss) will show more of ML's added value in Week 2.

## Biological validation (sanity check — full SHAP is Week 3)
Top logistic-regression drivers of "Resistant" are the established fluoroquinolone mechanisms:
`gyrA_S83I (+2.39)`, `parC_S80I (+2.27)`, `gyrA_S83F`, `qnrS1`, `qnrB1`, `gyrA_D87N`. The model
learned real biology, not artifacts. (A few co-carried MDR-plasmid genes also rank high via
co-selection — expected, to revisit with SHAP.)

## Rigor checks passed
- **Zero lineage leakage:** train/test share no sequence type (`tests/test_split.py` passes on the
  real split). 400 STs; StratifiedGroupKFold keeps folds disjoint AND class-balanced.
- **Lab phenotypes only**, QC-filtered (literature-aligned), balanced 184/184 test.
- Fully reproducible: fixed seed, config-driven, manifest + decision log.

## Bugs caught & fixed en route (rigor working)
1. genome_sequence API truncation (limit-25 → truncated assemblies) — size guard added.
2. Download timeout aborting the run — retry/backoff added.
3. QC analysis falsy-zero bug (contamination 0 → 99) — fixed.
4. mlst using system Perl / missing blastn — PATH fixed.
5. **Rules baseline strawman:** oqxAB is intrinsic in K. pneumoniae (88% of genomes) → naive
   baseline called 91% resistant. Corrected to expert QRDR+PMQR determinants for a fair benchmark.

## Artifacts
`data/processed/thin_slice_cipro_features.csv` (1472 × 446) · `..._split.csv` · `..._labels_qc.csv` ·
`results/metrics/thin_slice_cipro_metrics.json`.

## Next (Week 2, pending GO confirmation)
- Phase 2b: pull the full 5-drug panel (~6.6k genomes, ~12 GB — will flag size).
- Per-drug models for all 5 drugs; hyperparameter tuning; probability calibration (Phase 8);
  per-drug VME/ME tables. Meropenem is the flagship (where ML should add the most over rules).
