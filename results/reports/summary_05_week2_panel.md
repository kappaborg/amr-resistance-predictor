# Summary #5 — Week 2: Full-Panel Per-Drug Models

One model per antibiotic, reusing the already-annotated 1472-genome set (no new download/annotation). Each drug: lineage-aware split (StratifiedGroupKFold on MLST), expert rules baseline vs logistic regression vs calibrated XGBoost. **VME/ME lead.**

VME = resistant called susceptible (dangerous miss); ME = susceptible called resistant.


## How to read this

- **ROC-AUC / PR-AUC** are threshold-independent — they measure discrimination (can the model rank R above S on unseen lineages?).
- **VME / ME** are at a chosen operating point: the threshold is picked on out-of-fold TRAIN predictions to cap VME at the clinical target, then evaluated on the held-out TEST lineages. Test VME may exceed the target when data is thin (few resistant genomes → noisy threshold) — that generalization gap is reported honestly, not hidden.
- Lowering the threshold to cut VME **necessarily raises ME** — this is the clinical trade-off, not a defect.


## meropenem  (test 881: 355R/526S, 132 unseen lineages)

| Model | ROC-AUC | PR-AUC | VME | ME | Brier |
|---|---|---|---|---|---|
| Rules baseline | 0.934 | 0.865 | 0.054 | 0.078 | — |
| Logistic regression | 0.982 | 0.976 | 0.017 | 0.293 | — |
| XGBoost (calibrated) | 0.980 | 0.975 | 0.014 | 0.262 | 0.051 |

## gentamicin  (test 929: 401R/528S, 141 unseen lineages)

| Model | ROC-AUC | PR-AUC | VME | ME | Brier |
|---|---|---|---|---|---|
| Rules baseline | 0.967 | 0.938 | 0.030 | 0.036 | — |
| Logistic regression | 0.968 | 0.957 | 0.032 | 0.129 | — |
| XGBoost (calibrated) | 0.973 | 0.955 | 0.030 | 0.127 | 0.078 |

## ciprofloxacin  (test 793: 568R/225S, 139 unseen lineages)

| Model | ROC-AUC | PR-AUC | VME | ME | Brier |
|---|---|---|---|---|---|
| Rules baseline | 0.904 | 0.930 | 0.005 | 0.187 | — |
| Logistic regression | 0.973 | 0.986 | 0.033 | 0.089 | — |
| XGBoost (calibrated) | 0.971 | 0.984 | 0.033 | 0.098 | 0.051 |

## trimethoprim_sulfamethoxazole  (test 738: 512R/226S, 123 unseen lineages)

| Model | ROC-AUC | PR-AUC | VME | ME | Brier |
|---|---|---|---|---|---|
| Rules baseline | 0.875 | 0.902 | 0.029 | 0.221 | — |
| Logistic regression | 0.970 | 0.981 | 0.039 | 0.128 | — |
| XGBoost (calibrated) | 0.975 | 0.985 | 0.031 | 0.155 | 0.056 |

## cefoxitin  (test 634: 381R/253S, 77 unseen lineages)

| Model | ROC-AUC | PR-AUC | VME | ME | Brier |
|---|---|---|---|---|---|
| Rules baseline | 0.517 | 0.615 | 0.966 | 0.000 | — |
| Logistic regression | 0.906 | 0.934 | 0.018 | 0.850 | — |
| XGBoost (calibrated) | 0.935 | 0.965 | 0.013 | 0.854 | 0.127 |
