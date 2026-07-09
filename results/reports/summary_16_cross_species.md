# Summary #16 — Cross-Species Transfer (zero-shot across organisms)

Train a per-drug model on one organism, test on the other, over the **159 shared determinants**. Zero target-organism training. ROC-AUC.

| Drug | K.pneu→E.coli | E.coli→K.pneu | (within K.pneu) | (within E.coli) |
|---|---|---|---|---|
| ciprofloxacin | **0.979** | **0.974** | 0.989 | 0.995 |
| gentamicin | **0.954** | **0.962** | 0.975 | 0.970 |
| trimethoprim_sulfamethoxazole | **0.947** | **0.968** | 0.970 | 0.967 |

**Reading:** a cross-species ROC-AUC well above 0.5 means the determinant→phenotype mapping the model learned on one organism **transfers to a different species without retraining** — strong evidence it captures real mechanism, not organism-specific lineage structure. Gaps vs the within-organism reference quantify the species shift.
