# Summary #25 — Population-Structure Leakage, Quantified

**Date:** 2026-07-10 · *Klebsiella pneumoniae* · same calibrated XGBoost, three evaluation protocols · 5-fold where applicable · temporal cutoff 2014.

CLAUDE.md non-negotiable #1 is *split by lineage, never randomly*. Here is what that discipline is worth: the **same model** scored three ways. A random split lets a clone appear in train and test (leakage); the lineage-grouped split forces generalization to unseen lineages; the temporal split is a prospective-style check.

| Drug | n (lineages) | Random split AUC | **Lineage-held-out AUC** | Temporal AUC | ΔAUC inflation (random − lineage) |
|---|---|---|---|---|---|
| meropenem | 3529 (526) | 0.978 | **0.968** | 0.948 | **+0.010** |
| gentamicin | 3713 (543) | 0.983 | **0.981** | 0.973 | **+0.002** |
| ciprofloxacin | 3176 (508) | 0.989 | **0.983** | 0.982 | **+0.006** |
| trimethoprim_sulfamethoxazole | 2966 (460) | 0.978 | **0.977** | 0.972 | **+0.001** |
| cefoxitin | 2537 (337) | 0.935 | **0.903** | 0.896 | **+0.032** |

**Reading.** Every drug's random-split AUC is **inflated** relative to the honest lineage-held-out AUC (mean inflation **+0.010**). That gap is exactly the population-structure leakage the field warns about (*Biased sampling driven by bacterial population structure confounds ML prediction of AMR*, PLOS Biology 2024): because resistance is lineage-associated and bacteria are clonal, a random split lets the model memorize lineage identity instead of learning mechanism. Our headline numbers everywhere else use the **lineage-held-out** column — the smaller, honest one. The temporal column shows the model also holds up predicting *future* isolates from *past* ones. A project reporting only the random-split number would look better and mean less.
