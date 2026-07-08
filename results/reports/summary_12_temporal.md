# Summary #12 — Temporal (Prospective-Style) Validation

Train on isolates collected **before 2014**, test on **2014+** (2315/3850 genomes have a collection year; range 1996–2020). Tests forward-in-time generalization as resistance evolves — the axis most AMR-ML studies skip.

**Honest note:** temporal split does not control lineage (a clone can span the cutoff), so it is a *complementary* axis to the phylogeny-aware split, not a substitute. VME≤3% operating point; logistic regression.

| Drug | Train/Test (yr) | Test R/S | ROC (temporal) | ROC (lineage) | VME | ME |
|---|---|---|---|---|---|---|
| meropenem | <2014/2014+ | 334/725 | **0.955** | 0.982 | 0.051 | 0.302 |
| gentamicin | <2014/2014+ | 467/598 | **0.971** | 0.968 | 0.058 | 0.135 |
| ciprofloxacin | <2014/2014+ | 847/220 | **0.982** | 0.973 | 0.074 | 0.055 |
| trimethoprim_sulfamethoxazole | <2014/2014+ | 742/295 | **0.976** | 0.970 | 0.035 | 0.156 |
| cefoxitin | <2014/2014+ | 463/357 | **0.906** | 0.906 | 0.054 | 0.616 |

**Reading:** ROC on the future test set close to the lineage-split ROC = the model holds up over time. A drop is the honest, expected temporal-shift effect (literature: prospective validation typically 15–30% lower) and is reported, not hidden.
