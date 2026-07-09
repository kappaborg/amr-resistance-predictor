# Summary #19 — MIC Regression (predicting the resistance *level*)

Per-drug regression of continuous **log2(MIC, mg/L)** from determinant features, under the phylogeny-aware split. **Essential Agreement (EA)** = predictions within ±1 two-fold dilution (the CLSI/FDA metric for genotype-based MIC prediction); higher is better. RMSE in doubling dilutions; Pearson r.

| Drug | n (MIC) | test | EA (±1 dilution) | RMSE (dilutions) | Pearson r |
|---|---|---|---|---|---|
| meropenem | 2526 | 620 | **59.4%** | 1.38 | 0.811 |
| gentamicin | 2503 | 718 | **81.6%** | 0.98 | 0.753 |
| ciprofloxacin | 2498 | 623 | **88.9%** | 0.63 | 0.847 |
| cefoxitin | 2252 | 562 | **83.6%** | 0.76 | 0.822 |

**Reading:** EA is the clinical gold-standard for MIC prediction — an EA around or above the ~90% CLSI expectation means the predicted MIC lands within one doubling dilution of the lab value, i.e. the model predicts not just *whether* a strain resists but *how strongly*. This is a richer, quantitative output than binary R/S, evaluated by the metric regulators actually use.
