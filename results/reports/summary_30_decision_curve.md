# Summary #30 — Decision-Curve Analysis (clinical net benefit)

**Date:** 2026-07-17 · *Klebsiella pneumoniae* · net benefit vs threshold probability p_t from the pooled lineage-held-out **calibrated** predictions. A low p_t encodes "missing resistance is far costlier than over-treating" — the project's VME-first stance. The model has clinical utility over the p_t range where its curve sits highest.

| Drug | model > rules over p_t | model > treat-all over p_t | NB(model) @ p_t=0.1 |
|---|---|---|---|
| meropenem | 0.01–0.50 | 0.03–0.50 | 0.373 |
| gentamicin | 0.01–0.50 | 0.01–0.50 | 0.411 |
| ciprofloxacin | 0.33–0.50 | 0.03–0.50 | 0.706 |
| trimethoprim_sulfamethoxazole | 0.01–0.50 | 0.02–0.50 | 0.678 |
| cefoxitin | 0.01–0.50 | 0.18–0.50 | 0.555 |

**Reading.** Where the **model** curve is above both **treat-all** and **rules**, using the model to guide therapy yields more net benefit (more true resistant caught per false-positive over-treatment) than treating everyone or trusting a gene lookup. The advantage concentrates in the **low-p_t** region — exactly where a clinician who fears missing resistance operates — which is the clinical case for the model beyond raw AUC.
