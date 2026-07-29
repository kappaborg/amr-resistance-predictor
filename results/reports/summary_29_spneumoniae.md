# Summary #29 — Streptococcus pneumoniae added (organism #8)

**Date:** 2026-07-15 · Fourth Round-4 organism; third Gram-positive. WHO priority; heavily sequenced.

## Data
- **3,877 genomes** (2.43 GB, 23 truncated pulls correctly rejected) → **3,857 QC-passed** (length
  1.9–2.4 Mbp). AMRFinderPlus `--organism Streptococcus_pneumoniae` + MLST `spneumoniae`.
- **Panel: 4 drugs** — penicillin, erythromycin, tetracycline, TMP-SMX. Fluoroquinolones dropped
  (pneumococcal FQ-R is genuinely rare: cipro 3 R, levofloxacin 2 R).
- **Truncation-floor bug fixed en route:** the fixed 2.5 MB download floor (for K. pneumoniae's ~5.5 Mbp
  genomes) was falsely rejecting valid ~2.1 Mbp pneumococcal genomes; the floor now scales per organism
  (S. pneumoniae → 1.04 MB). Would also have entirely blocked C. jejuni.

## Results (unseen-lineage test, VME ≤ 3% operating point)
| Drug | ROC-AUC | PR-AUC | VME | ME | Rules ROC | Reading |
|---|---|---|---|---|---|---|
| **penicillin** | 0.961 | 0.923 | 0.066 | 0.058 | **0.860** | pbp-mosaic baseline (pre-registered `ORG_RULES`) works; ML +0.10 |
| erythromycin | 0.983 | 0.983 | 0.043 | 0.000 | 0.979 | erm/mef = mechanism; ML ≈ rules, well-balanced |
| tetracycline | 0.969 | 0.962 | 0.000 | 0.887 | 0.948 | tet genes; VME-safety over-calls (ME 0.89) |
| TMP-SMX | 0.858 | 0.843 | 0.067 | 0.403 | 0.500 | **co-selection-driven — see caveat** |

## Honesty checks
- **Penicillin — the pre-registered override paid off.** `(spneumoniae, penicillin) → PBP1a/2b/2x`
  (added in #45 before this organism existed) gives a real baseline of **0.860** (not the 0.500 a
  β-lactamase list would give); ML at 0.961 refines the specific mosaic patterns. Prediction drivers
  confirm it (`pbp2b`, `pbp2x`). Clean, honest ML-adds result.
- **TMP-SMX is a genuine scope limitation, NOT a fixable baseline.** Pneumococcal co-trimoxazole
  resistance is via **folA/folP chromosomal mutations**, which **AMRFinderPlus does not catalog**
  (folA/folP/dfr all 0% in the feature matrix). So neither the rules baseline (0.500 — the acquired
  sul/dfr genes are absent) **nor ML** sees the causal mechanism. ML's 0.858 comes from **co-selection**
  with correlated MDR determinants (erm/tet in resistant lineages), not the mechanism. Reported as a
  co-selection result, *not* as "ML adds +0.36" — an honest determinant-catalog boundary (same class as
  the k-mer/novel-gene scope limit).
- **No reliable species marker** (like E. coli/Salmonella) — species relies on the distinct
  `spneumoniae` MLST scheme; no collision with other organisms' markers.

## Status
- ✅ In the demo (8 organisms: + **Streptococcus pneumoniae**) + species check; 4 models saved; 13 tests pass.
- Remaining Round-4 queue: **C. jejuni** (last one).
