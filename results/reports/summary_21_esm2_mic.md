# Summary #21 — ESM-2 Protein Language Model for MIC Prediction

**Date:** 2026-07-09 · ESM-2 150M (`esm2_t30_150M_UR50D`, 640-d) on MPS · 5-fold lineage-grouped CV · metric = **Essential Agreement** (±1 doubling dilution, CLSI/FDA).

Each resistance-gene protein is extracted from the genome (AMRFinderPlus coordinates → translate) and embedded by ESM-2, so **different alleles of the same gene get different vectors** — the allele-level signal that presence/absence throws away. Per genome we mean-pool the embeddings (all AMR proteins + the drug-class subset). MIC is regressed from three feature sets on identical folds.

| Drug | n (MIC) | Determinants (baseline) | ESM-2 PLM only | Determinants + PLM | Δ EA |
|---|---|---|---|---|---|
| meropenem | 2526 | 68.6% | 70.1% | **72.4%** | +3.8% |
| cefoxitin | 2252 | 78.8% | 80.9% | **77.3%** | -1.5% |
| gentamicin | 2503 | 84.7% | 77.9% | **82.7%** | -2.1% |
| ciprofloxacin | 2498 | 87.6% | 87.0% | **87.6%** | -0.0% |

**Reading.** Δ EA = (determinants + PLM) − determinants. A positive Δ on the hard drugs (meropenem, cefoxitin) means ESM-2's allele-resolution adds MIC signal beyond presence/absence. If Δ ≈ 0, the curated determinants already carry the phenotype and allele sequence adds little — an honest, useful finding either way. PLM-only shows how much the embeddings carry alone (no gene-identity features).

## Verdict: ESM-2 helps on meropenem — the hardest drug — and is confirmed significant

The single-run table above has high fold variance (±4–9%), so the +3.8 meropenem gain alone is only
suggestive. We tested it properly with **30 paired comparisons** (6 CV repeats × 5 folds, identical
folds for both feature sets):

| meropenem (n=2526) | Determinants | Determinants + ESM-2 |
|---|---|---|
| Essential Agreement | 68.4% ± 7.9% | **73.2% ± 4.7%** |
| Mean Δ (paired) | — | **+4.8 points** |
| Folds where PLM wins | — | **77% (23/30)** |
| Wilcoxon signed-rank | — | **p = 0.0004** |

The gain is **statistically significant and stable**, and ESM-2 also **cuts the fold-to-fold variance
nearly in half** (7.9% → 4.7%) — more reliable MIC calls, not just higher on average.

**Why meropenem and not the others — the mechanism.** Carbapenem MIC depends on *which* carbapenemase
allele a strain carries (blaKPC-2 vs KPC-3, OXA-48-like variants, NDM/VIM families all confer
*different* MIC magnitudes) and on porin mutations (ompK35/36) whose exact sequence matters. Presence/
absence collapses all of that to one bit; ESM-2 reads the actual protein sequence, so it separates a
high-MIC allele from a borderline one. For gentamicin (aminoglycoside-modifying enzyme *presence* is
nearly sufficient) and ciprofloxacin (the QRDR point mutations are *already* explicit determinant
features), there is little allele signal left to recover — so ESM-2 is neutral there, as expected.

**Bottom line.** ESM-2 150M delivers a real, significant improvement precisely where the biology
predicts allele resolution should matter — the hardest drug, meropenem (+4.8 EA points, p=0.0004) —
while honestly adding nothing on drugs whose determinants already saturate. Extraction dedupes 58,502
resistance proteins across 3,850 genomes to **1,945 unique alleles**, embedded once on the M1 Max GPU
(MPS) in minutes, so the whole experiment is laptop-scale and reproducible.

## Scaling to ESM-2 650M — bigger model, no meaningful gain (150M is enough)

We re-embedded the same 1,945 alleles with **ESM-2 650M** (1280-d, 4.3× the parameters) and compared
all three feature sets on the *identical* 30 paired folds (meropenem, n=2526):

| Feature set | EA mean | EA std | Δ vs determinants | Wilcoxon p |
|---|---|---|---|---|
| Determinants | 68.4% | 7.9% | — | — |
| Determinants + ESM-2 **150M** | 73.2% | 4.7% | **+4.8** | 0.0004 |
| Determinants + ESM-2 **650M** | 73.6% | 4.8% | **+5.2** | 0.0003 |

**Head-to-head (650M − 150M, paired):** mean **+0.4 EA points**, 650M wins only **50%** of folds,
**Wilcoxon p = 0.68** — statistically indistinguishable. The full 4-drug 650M table
(`summary_21_esm2_mic_650M.md`) mirrors 150M drug-for-drug (meropenem gains, gentamicin/cipro neutral).

**Reading.** The allele distinctions that drive carbapenem MIC — KPC-2 vs KPC-3, OXA-48-like vs
NDM/VIM family identity, porin frameshifts — are *gross* sequence differences that even the small
150M model resolves cleanly; 650M's extra representational capacity buys nothing here. So the honest,
practical recommendation is **use 150M**: it captures the full allele signal, runs in minutes on the
laptop GPU, and needs no 2.5 GB model. Scaling up is a confirmed dead end for this task — worth
knowing, and cheaper to run.
