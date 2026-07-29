# Summary #21 — ESM-2 650M Protein Language Model for MIC Prediction

**Date:** 2026-07-09 · ESM-2 650M (`esm2_t33_650M_UR50D`, 1280-d) on MPS · 5-fold lineage-grouped CV · metric = **Essential Agreement** (±1 doubling dilution, CLSI/FDA).

Each resistance-gene protein is extracted from the genome (AMRFinderPlus coordinates → translate) and embedded by ESM-2, so **different alleles of the same gene get different vectors** — the allele-level signal that presence/absence throws away. Per genome we mean-pool the embeddings (all AMR proteins + the drug-class subset). MIC is regressed from three feature sets on identical folds.

| Drug | n (MIC) | Determinants (baseline) | ESM-2 PLM only | Determinants + PLM | Δ EA |
|---|---|---|---|---|---|
| meropenem | 2526 | 68.6% | 72.5% | **73.5%** | +4.9% |
| cefoxitin | 2252 | 78.8% | 81.6% | **79.3%** | +0.5% |
| gentamicin | 2503 | 84.7% | 78.3% | **80.2%** | -4.5% |
| ciprofloxacin | 2498 | 87.6% | 87.1% | **87.1%** | -0.6% |

**Reading.** Δ EA = (determinants + PLM) − determinants. A positive Δ on the hard drugs (meropenem, cefoxitin) means ESM-2's allele-resolution adds MIC signal beyond presence/absence. If Δ ≈ 0, the curated determinants already carry the phenotype and allele sequence adds little — an honest, useful finding either way. PLM-only shows how much the embeddings carry alone (no gene-identity features).
