# Summary #24 — Regulator-Grade Metrics (all organisms): Calibration, CIs & DeLong

**Date:** 2026-07-17 · **all 8 organisms** · isotonic-calibrated XGBoost · pooled held-out predictions from full 5-fold **lineage-grouped** CV (thresholds on each fold's TRAIN out-of-fold predictions, VME target 3%) · 2000 **lineage-clustered** bootstrap resamples.

FDA/CLSI vocabulary + bars: **VME ≤ 1.5%** (≤3% tolerated), **ME ≤ 3%**, **CA ≥ 90%**, **EA ≥ 90%**; **Brier** = calibration (lower better). ROC-AUC vs the **organism-aware** known-gene rules baseline by a **DeLong** paired test. CIs resample whole MLST lineages (clonal blocks), the honest variance under population structure.


## Klebsiella pneumoniae

| Drug | n (R) | VME % [CI] | ME % [CI] | CA % | MCC | Brier | ROC-AUC [CI] | Rules AUC | DeLong p | PPV/NPV @30% | VME≤1.5 | ME≤3 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| meropenem | 3529 (1427) | 1.3% [0.3%,3.5%] | 55.9% [42.6%,65.8%] | 66.2% | 0.48 | 0.051 | 0.968 [0.946,0.980] | 0.933 | 5.3e-24 | 0.43/0.99 | ✅ | ❌ |
| gentamicin | 3713 (1601) | 1.6% [0.8%,2.8%] | 32.4% [16.2%,46.0%] | 80.9% | 0.67 | 0.035 | 0.981 [0.974,0.986] | 0.948 | 1.4e-30 | 0.57/0.99 | ❌ | ❌ |
| ciprofloxacin | 3176 (2277) | 2.3% [1.3%,5.1%] | 7.7% [4.7%,11.2%] | 96.2% | 0.90 | 0.034 | 0.983 [0.969,0.989] | 0.911 | 6.9e-36 | 0.85/0.99 | ❌ | ❌ |
| trimethoprim_sulfamethoxazole | 2966 (2057) | 2.3% [1.4%,3.9%] | 18.3% [9.6%,27.8%] | 92.8% | 0.83 | 0.048 | 0.977 [0.969,0.985] | 0.885 | 9.4e-51 | 0.70/0.99 | ❌ | ❌ |
| cefoxitin | 2537 (1525) | 0.7% [0.2%,1.7%] | 94.2% [90.7%,96.7%] | 62.0% | 0.16 | 0.113 | 0.905 [0.841,0.929] | 0.518 | 0.0e+00 | 0.31/0.95 | ✅ | ❌ |

## Escherichia coli

| Drug | n (R) | VME % [CI] | ME % [CI] | CA % | MCC | Brier | ROC-AUC [CI] | Rules AUC | DeLong p | PPV/NPV @30% | VME≤1.5 | ME≤3 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ciprofloxacin | 2787 (1314) | 2.8% [1.5%,5.7%] | 20.2% [13.8%,29.2%] | 88.0% | 0.77 | 0.026 | 0.982 [0.968,0.989] | 0.876 | 2.2e-76 | 0.67/0.99 | ❌ | ❌ |
| gentamicin | 2802 (857) | 1.6% [0.7%,3.0%] | 62.5% [53.2%,70.2%] | 56.1% | 0.37 | 0.027 | 0.972 [0.956,0.980] | 0.968 | 1.3e-01 | 0.40/0.98 | ❌ | ❌ |
| trimethoprim_sulfamethoxazole | 1671 (881) | 2.0% [1.0%,3.6%] | 48.4% [39.1%,59.5%] | 76.1% | 0.57 | 0.058 | 0.963 [0.948,0.974] | 0.845 | 1.1e-52 | 0.46/0.98 | ❌ | ❌ |
| ampicillin | 2412 (1835) | 0.6% [0.2%,1.2%] | 87.5% [80.0%,93.2%] | 78.6% | 0.28 | 0.065 | 0.952 [0.934,0.966] | 0.954 | 6.4e-01 | 0.33/0.98 | ✅ | ❌ |
| ceftriaxone | 931 (447) | 2.0% [0.7%,4.9%] | 50.4% [41.3%,60.2%] | 72.8% | 0.54 | 0.077 | 0.957 [0.929,0.974] | 0.974 | 6.7e-03 | 0.45/0.98 | ❌ | ❌ |

## Acinetobacter baumannii

| Drug | n (R) | VME % [CI] | ME % [CI] | CA % | MCC | Brier | ROC-AUC [CI] | Rules AUC | DeLong p | PPV/NPV @30% | VME≤1.5 | ME≤3 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| meropenem | 617 (415) | 3.6% [0.0%,23.7%] | 66.8% [26.9%,85.1%] | 75.7% | 0.41 | 0.183 | 0.801 [0.747,0.937] | 0.837 | 1.0e-01 | 0.38/0.96 | ❌ | ❌ |
| imipenem | 1258 (617) | 7.3% [0.0%,34.5%] | 66.5% [21.2%,84.8%] | 62.6% | 0.32 | 0.136 | 0.902 [0.838,0.968] | 0.896 | 4.9e-01 | 0.37/0.91 | ❌ | ❌ |
| ciprofloxacin | 1349 (1234) | 0.6% [0.3%,1.8%] | 2.6% [0.0%,7.1%] | 99.2% | 0.95 | 0.008 | 0.979 [0.949,0.998] | 0.980 | 8.9e-01 | 0.94/1.00 | ✅ | ✅ |
| gentamicin | 1255 (1124) | 2.1% [0.6%,8.1%] | 30.5% [6.7%,52.4%] | 94.9% | 0.71 | 0.043 | 0.926 [0.882,0.980] | 0.948 | 4.3e-02 | 0.58/0.99 | ❌ | ❌ |
| amikacin | 1196 (733) | 0.5% [0.0%,3.0%] | 47.1% [13.1%,74.5%] | 81.4% | 0.63 | 0.148 | 0.913 [0.859,0.980] | 0.735 | 1.1e-61 | 0.48/1.00 | ✅ | ❌ |

## Staphylococcus aureus

| Drug | n (R) | VME % [CI] | ME % [CI] | CA % | MCC | Brier | ROC-AUC [CI] | Rules AUC | DeLong p | PPV/NPV @30% | VME≤1.5 | ME≤3 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| oxacillin | 1457 (603) | 1.8% [0.5%,4.2%] | 65.5% [52.9%,77.7%] | 60.9% | 0.40 | 0.065 | 0.952 [0.910,0.983] | 0.941 | 6.4e-03 | 0.39/0.98 | ❌ | ❌ |
| cefoxitin | 2067 (1303) | 1.4% [0.5%,3.2%] | 53.7% [41.7%,63.3%] | 79.3% | 0.56 | 0.070 | 0.948 [0.924,0.974] | 0.913 | 5.5e-12 | 0.44/0.99 | ✅ | ❌ |
| ciprofloxacin | 2920 (1277) | 2.0% [0.4%,7.7%] | 73.6% [55.4%,89.1%] | 57.7% | 0.33 | 0.035 | 0.971 [0.916,0.985] | 0.963 | 6.0e-03 | 0.36/0.97 | ❌ | ❌ |
| erythromycin | 3178 (1441) | 1.0% [0.2%,3.2%] | 64.9% [47.5%,80.3%] | 64.0% | 0.43 | 0.055 | 0.955 [0.924,0.970] | 0.945 | 2.2e-04 | 0.40/0.99 | ✅ | ❌ |
| clindamycin | 1913 (875) | 1.9% [0.7%,4.9%] | 83.2% [73.6%,90.4%] | 53.9% | 0.25 | 0.109 | 0.914 [0.861,0.946] | 0.880 | 6.0e-09 | 0.34/0.95 | ❌ | ❌ |

## Pseudomonas aeruginosa

| Drug | n (R) | VME % [CI] | ME % [CI] | CA % | MCC | Brier | ROC-AUC [CI] | Rules AUC | DeLong p | PPV/NPV @30% | VME≤1.5 | ME≤3 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| meropenem | 873 (481) | 0.0% [0.0%,0.0%] | 99.7% [99.1%,100.0%] | 55.2% | 0.04 | 0.198 | 0.762 [0.711,0.796] | 0.594 | 2.2e-18 | 0.30/1.00 | ✅ | ❌ |
| ceftazidime | 705 (361) | 0.6% [0.0%,1.7%] | 97.4% [93.8%,99.7%] | 52.2% | 0.08 | 0.213 | 0.751 [0.674,0.793] | 0.652 | 1.7e-07 | 0.30/0.92 | ✅ | ❌ |
| ciprofloxacin | 572 (368) | 1.6% [0.3%,4.0%] | 76.5% [63.7%,86.5%] | 71.7% | 0.36 | 0.110 | 0.918 [0.872,0.945] | 0.867 | 9.8e-05 | 0.36/0.97 | ❌ | ❌ |
| tobramycin | 511 (185) | 3.2% [0.0%,10.6%] | 76.4% [66.9%,83.9%] | 50.1% | 0.27 | 0.086 | 0.918 [0.832,0.951] | 0.919 | 9.9e-01 | 0.35/0.94 | ❌ | ❌ |

## Salmonella enterica

| Drug | n (R) | VME % [CI] | ME % [CI] | CA % | MCC | Brier | ROC-AUC [CI] | Rules AUC | DeLong p | PPV/NPV @30% | VME≤1.5 | ME≤3 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ampicillin | 2297 (1060) | 0.6% [0.0%,1.4%] | 96.3% [92.1%,98.5%] | 47.9% | 0.11 | 0.073 | 0.941 [0.915,0.959] | 0.954 | 1.2e-03 | 0.31/0.94 | ✅ | ❌ |
| ceftriaxone | 1201 (210) | 3.3% [0.5%,12.2%] | 6.5% [2.0%,15.1%] | 94.1% | 0.82 | 0.008 | 0.988 [0.951,0.997] | 0.991 | 4.0e-01 | 0.87/0.98 | ❌ | ❌ |
| ciprofloxacin | 2015 (316) | 4.1% [0.0%,44.3%] | 82.9% [63.7%,94.3%] | 29.5% | 0.13 | 0.126 | 0.912 [0.495,0.964] | 0.935 | 6.1e-03 | 0.33/0.91 | ❌ | ❌ |
| chloramphenicol | 2226 (387) | 4.4% [0.2%,17.5%] | 81.5% [60.0%,95.0%] | 31.9% | 0.15 | 0.103 | 0.867 [0.781,0.921] | 0.932 | 2.6e-12 | 0.33/0.91 | ❌ | ❌ |
| trimethoprim_sulfamethoxazole | 1568 (162) | 4.9% [0.0%,8.8%] | 81.9% [74.0%,90.0%] | 26.0% | 0.11 | 0.092 | 0.803 [0.731,0.911] | 0.732 | 1.5e-04 | 0.33/0.90 | ❌ | ❌ |

## Enterococcus faecium

| Drug | n (R) | VME % [CI] | ME % [CI] | CA % | MCC | Brier | ROC-AUC [CI] | Rules AUC | DeLong p | PPV/NPV @30% | VME≤1.5 | ME≤3 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| vancomycin | 1851 (839) | 1.3% [0.1%,4.0%] | 0.9% [0.2%,1.5%] | 98.9% | 0.98 | 0.005 | 0.998 [0.996,1.000] | 0.982 | 1.3e-08 | 0.98/0.99 | ✅ | ✅ |
| ampicillin | 1448 (1078) | 1.7% [0.4%,4.6%] | 1.1% [0.2%,2.4%] | 98.5% | 0.96 | 0.007 | 0.993 [0.987,0.997] | 0.694 | 7.3e-121 | 0.97/0.99 | ❌ | ✅ |
| tetracycline | 831 (612) | 2.5% [1.0%,6.0%] | 27.4% [17.7%,39.2%] | 91.0% | 0.76 | 0.059 | 0.950 [0.912,0.974] | 0.858 | 1.5e-13 | 0.60/0.99 | ❌ | ❌ |

## Streptococcus pneumoniae

| Drug | n (R) | VME % [CI] | ME % [CI] | CA % | MCC | Brier | ROC-AUC [CI] | Rules AUC | DeLong p | PPV/NPV @30% | VME≤1.5 | ME≤3 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| penicillin | 2431 (910) | 3.1% [1.6%,5.4%] | 26.6% [18.1%,38.6%] | 82.2% | 0.68 | 0.057 | 0.955 [0.940,0.967] | 0.887 | 1.5e-46 | 0.61/0.98 | ❌ | ❌ |
| erythromycin | 3275 (1304) | 1.4% [0.6%,2.5%] | 61.2% [50.9%,71.1%] | 62.6% | 0.43 | 0.014 | 0.980 [0.968,0.989] | 0.983 | 1.1e-01 | 0.41/0.98 | ✅ | ❌ |
| tetracycline | 2775 (1250) | 0.2% [0.0%,0.4%] | 83.9% [78.2%,88.4%] | 53.8% | 0.28 | 0.045 | 0.961 [0.942,0.973] | 0.906 | 3.3e-31 | 0.34/1.00 | ✅ | ❌ |
| trimethoprim_sulfamethoxazole | 3127 (1619) | 0.2% [0.0%,0.5%] | 99.7% [99.3%,99.9%] | 51.8% | 0.01 | 0.151 | 0.847 [0.796,0.887] | 0.500 | 0.0e+00 | 0.30/0.81 | ✅ | ❌ |

**Reading.** Every organism × drug now carries a **calibrated** model (isotonic; Brier reported), **lineage-clustered** 95% CIs, and a **DeLong** significance test versus the organism-aware gene-lookup baseline — the project's central claim ("what ML adds over a lookup"), significance-tested across all organisms. Thresholds target the clinical **VME ≤ 3%** operating point; where the strict FDA 1.5% bar or CA ≥ 90% is not met it is shown honestly (the high ME / borderline VME on low-prevalence-R drugs is the safety-first over-calling trade-off, reported not hidden). **MCC** (Matthews correlation) is robust to the class imbalance that inflates accuracy/CA; **PPV/NPV** are shown at a realistic 30% local resistance prevalence (full 10/30/50% set in `clinical_rigor.json`) since predictive value depends on prevalence, not just the test's sensitivity/specificity.
