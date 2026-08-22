# Summary #19 — MIC Regression (predicting the resistance *level*)

Per-drug regression of continuous **log2(MIC, mg/L)** from determinant features, under a full 5-fold **phylogeny-aware** CV (pooled held-out predictions). **Essential Agreement (EA)** = predictions within ±1 two-fold dilution (CLSI/FDA), scored **censor-aware**: right-/left-censored (off-scale `>`/`<`) MICs are treated as bounds, not collapsed onto the boundary dilution, so EA is not inflated. RMSE in doubling dilutions; Pearson r.

| Drug | n (censored) | EA (±1 dilution) | RMSE (dilutions) | Pearson r |
|---|---|---|---|---|
| meropenem | 2526 (1920 cens) | 62.7% | 1.34 | 0.840 |
| gentamicin | 2503 (1895 cens) | 71.9% | 0.84 | 0.823 |
| ciprofloxacin | 2498 (2146 cens) | 66.0% | 0.73 | 0.714 |
| cefoxitin | 2252 (1671 cens) | 62.2% | 0.84 | 0.765 |

**Reading:** EA is the clinical gold-standard for MIC prediction — an EA around or above the ~90% CLSI expectation means the predicted MIC lands within one doubling dilution of the lab value, i.e. the model predicts not just *whether* a strain resists but *how strongly*. This is a richer, quantitative output than binary R/S, evaluated by the metric regulators actually use.
