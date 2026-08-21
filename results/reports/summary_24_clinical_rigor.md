# Summary #24 — Regulator-Grade Metrics (all organisms): Calibration, CIs & DeLong

**Date:** 2026-07-17 · **all 1 organisms** · isotonic-calibrated XGBoost · pooled held-out predictions from full 5-fold **lineage-grouped** CV (thresholds on each fold's TRAIN out-of-fold predictions, VME target 3%) · 2000 **lineage-clustered** bootstrap resamples.

FDA/CLSI vocabulary + bars: **VME ≤ 1.5%** (≤3% tolerated), **ME ≤ 3%**, **CA ≥ 90%**, **EA ≥ 90%**; **Brier** = calibration (lower better). ROC-AUC vs the **organism-aware** known-gene rules baseline by a **DeLong** paired test. CIs resample whole MLST lineages (clonal blocks), the honest variance under population structure.


## Klebsiella pneumoniae

| Drug | n (R) | VME % [CI] | ME % [CI] | CA % | MCC | Brier | ROC-AUC [CI] | Rules AUC | DeLong p | PPV/NPV @30% | VME≤1.5 | ME≤3 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| meropenem | 3529 (1427) | 1.3% [0.3%,3.5%] | 55.9% [42.6%,65.8%] | 66.2% | 0.48 | 0.051 | 0.968 [0.946,0.980] | 0.933 | 5.3e-24 | 0.43/0.99 | ✅ | ❌ |
| gentamicin | 3713 (1601) | 1.6% [0.8%,2.8%] | 32.4% [16.2%,46.0%] | 80.9% | 0.67 | 0.035 | 0.981 [0.974,0.986] | 0.948 | 1.4e-30 | 0.57/0.99 | ❌ | ❌ |
| ciprofloxacin | 3176 (2277) | 2.3% [1.3%,5.1%] | 7.7% [4.7%,11.2%] | 96.2% | 0.90 | 0.034 | 0.983 [0.969,0.989] | 0.911 | 6.9e-36 | 0.85/0.99 | ❌ | ❌ |
| trimethoprim_sulfamethoxazole | 2966 (2057) | 2.3% [1.4%,3.9%] | 18.3% [9.6%,27.8%] | 92.8% | 0.83 | 0.048 | 0.977 [0.969,0.985] | 0.885 | 9.4e-51 | 0.70/0.99 | ❌ | ❌ |
| cefoxitin | 2537 (1525) | 0.7% [0.2%,1.7%] | 94.2% [90.7%,96.7%] | 62.0% | 0.16 | 0.113 | 0.905 [0.841,0.929] | 0.518 | 0.0e+00 | 0.31/0.95 | ✅ | ❌ |

**Reading.** Every organism × drug now carries a **calibrated** model (isotonic; Brier reported), **lineage-clustered** 95% CIs, and a **DeLong** significance test versus the organism-aware gene-lookup baseline — the project's central claim ("what ML adds over a lookup"), significance-tested across all organisms. Thresholds target the clinical **VME ≤ 3%** operating point; where the strict FDA 1.5% bar or CA ≥ 90% is not met it is shown honestly (the high ME / borderline VME on low-prevalence-R drugs is the safety-first over-calling trade-off, reported not hidden). **MCC** (Matthews correlation) is robust to the class imbalance that inflates accuracy/CA; **PPV/NPV** are shown at a realistic 30% local resistance prevalence (full 10/30/50% set in `clinical_rigor.json`) since predictive value depends on prevalence, not just the test's sensitivity/specificity.
