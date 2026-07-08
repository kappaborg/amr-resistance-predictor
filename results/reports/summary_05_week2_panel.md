# Summary #5 — Week 2: Full-Panel Per-Drug Models

One model per antibiotic, reusing the already-annotated 1472-genome set (no new download/annotation). Each drug: lineage-aware split (StratifiedGroupKFold on MLST), expert rules baseline vs logistic regression vs calibrated XGBoost. **VME/ME lead.**

VME = resistant called susceptible (dangerous miss); ME = susceptible called resistant.

---
## Interpretation (the honest story)

**The headline result — where ML earns its keep: cefoxitin.** The known-gene rules baseline is
*useless* here (ROC 0.53, VME **0.924** — it misses 92% of resistant strains). Reason: cefoxitin
resistance in K. pneumoniae is driven largely by **porin loss (ompK35/36)**, not just plasmid AmpC
genes — a single-gene lookup can't see it. Logistic regression captures it (ROC **0.957**, VME
0.114). This is the project's thesis proven on one drug: *what the model adds over the lookup.*

**ML clearly beats the baseline on 3/5 drugs** (cefoxitin, TMP-SMX, ciprofloxacin) — higher ROC and
large cuts in major errors. **On meropenem and gentamicin the rules baseline is already strong**
(carbapenemase / aminoglycoside-modifying-enzyme presence is a direct signal), and ML roughly matches
it on discrimination — an honest "the model doesn't dominate everywhere" finding, not spin.

**⚠️ A real problem to fix (Phase 8): the calibrated XGBoost has high VME at the default 0.5
threshold** (meropenem 0.233, gentamicin 0.265, cefoxitin 0.304) — it under-calls resistance, the
clinically dangerous error. This is a *threshold* artifact on imbalanced data, not a discrimination
failure (its ROC/PR are fine). **Logistic regression (balanced class weight) is currently the
better-behaved model** on VME. The fix is VME-bounded threshold selection (choose the operating point
that caps VME at a clinical tolerance), which is the immediate next step.

**Best model per drug is not always XGBoost** — for these sparse binary determinant features,
logistic regression is competitive or better. We report all three transparently.

---


## meropenem  (test 338: 73R/265S, 99 unseen lineages)

| Model | ROC-AUC | PR-AUC | VME | ME | Brier |
|---|---|---|---|---|---|
| Rules baseline | 0.949 | 0.872 | 0.082 | 0.019 | — |
| Logistic regression | 0.947 | 0.922 | 0.082 | 0.045 | — |
| XGBoost (calibrated) | 0.951 | 0.879 | 0.233 | 0.015 | 0.050 |

## gentamicin  (test 351: 102R/249S, 98 unseen lineages)

| Model | ROC-AUC | PR-AUC | VME | ME | Brier |
|---|---|---|---|---|---|
| Rules baseline | 0.929 | 0.861 | 0.118 | 0.024 | — |
| Logistic regression | 0.962 | 0.910 | 0.196 | 0.032 | — |
| XGBoost (calibrated) | 0.942 | 0.825 | 0.265 | 0.080 | 0.097 |

## ciprofloxacin  (test 368: 184R/184S, 106 unseen lineages)

| Model | ROC-AUC | PR-AUC | VME | ME | Brier |
|---|---|---|---|---|---|
| Rules baseline | 0.897 | 0.834 | 0.022 | 0.185 | — |
| Logistic regression | 0.976 | 0.976 | 0.049 | 0.043 | — |
| XGBoost (calibrated) | 0.975 | 0.967 | 0.038 | 0.120 | 0.059 |

## trimethoprim_sulfamethoxazole  (test 320: 179R/141S, 88 unseen lineages)

| Model | ROC-AUC | PR-AUC | VME | ME | Brier |
|---|---|---|---|---|---|
| Rules baseline | 0.895 | 0.865 | 0.039 | 0.170 | — |
| Logistic regression | 0.969 | 0.963 | 0.045 | 0.092 | — |
| XGBoost (calibrated) | 0.957 | 0.952 | 0.045 | 0.121 | 0.066 |

## cefoxitin  (test 168: 79R/89S, 47 unseen lineages)

| Model | ROC-AUC | PR-AUC | VME | ME | Brier |
|---|---|---|---|---|---|
| Rules baseline | 0.532 | 0.500 | 0.924 | 0.011 | — |
| Logistic regression | 0.957 | 0.961 | 0.114 | 0.079 | — |
| XGBoost (calibrated) | 0.875 | 0.873 | 0.304 | 0.135 | 0.140 |
