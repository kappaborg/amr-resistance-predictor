# Summary #25 — Population-Structure Leakage, Quantified

**Date:** 2026-07-10 · *Klebsiella pneumoniae* · same calibrated XGBoost, three evaluation protocols · 5-fold where applicable · temporal cutoff 2014.

CLAUDE.md non-negotiable #1 is *split by lineage, never randomly*. Here is what that discipline is worth: the **same model** scored three ways. A random split lets a clone appear in train and test (leakage); the lineage-grouped split forces generalization to unseen lineages; the temporal split is a prospective-style check.

| Drug | n (lineages) | Random split AUC | **Lineage-held-out AUC** | Temporal AUC | ΔAUC inflation (random − lineage) |
|---|---|---|---|---|---|
| meropenem | 3529 (526) | 0.978 | **0.968** | 0.962 | **+0.010** |
| gentamicin | 3713 (543) | 0.983 | **0.981** | 0.977 | **+0.002** |
| ciprofloxacin | 3176 (508) | 0.989 | **0.983** | 0.986 | **+0.006** |
| trimethoprim_sulfamethoxazole | 2966 (460) | 0.978 | **0.977** | 0.973 | **+0.001** |
| cefoxitin | 2537 (337) | 0.935 | **0.903** | 0.899 | **+0.032** |

**Reading.** Every drug's random-split AUC is **inflated** relative to the honest lineage-held-out AUC (mean inflation **+0.010**, up to **+0.032** for cefoxitin). That gap is the population-structure leakage the field warns about (*Biased sampling driven by bacterial population structure confounds ML prediction of AMR*, PLOS Biology 2024): because resistance is lineage-associated and bacteria are clonal, a random split lets the model memorize lineage identity instead of learning mechanism.

**Why our inflation is small — and why that is the point.** k-mer / whole-genome models in that literature can lose **0.1–0.2 AUC** when moved to a lineage-held-out split, because their features partly encode lineage. Ours barely move because **determinant presence/absence encodes mechanism, not ancestry** — a genome resists meropenem because it carries blaKPC, not because of its clade. So the small ΔAUC is itself evidence the model learned biology rather than population structure: a validation of the determinant-feature choice, not a disappointment. The leakage is still real (nonzero on every drug, largest where the signal is weakest — cefoxitin), so we still correctly report the **lineage-held-out** column as the headline everywhere. The temporal column confirms the model also holds up predicting *future* isolates from *past* ones (cutoff 2014). A project reporting only the random-split number would look better and mean less.
