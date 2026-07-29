# Summary #21 — ESM-2 150M for MIC Prediction · Acinetobacter baumannii

**Date:** 2026-07-10 · ESM-2 150M (`esm2_t30_150M_UR50D`, 640-d) on MPS · *Acinetobacter baumannii* (taxon 470) · 5-fold lineage-grouped CV · metric = **Essential Agreement** (±1 doubling dilution, CLSI/FDA).

Each resistance-gene protein is extracted from the genome (AMRFinderPlus coordinates → translate) and embedded by ESM-2, so **different alleles of the same gene get different vectors** — the allele-level signal that presence/absence throws away. Per genome we mean-pool the embeddings (all AMR proteins + the drug-class subset). MIC is regressed from three feature sets on identical folds.

| Drug | n (MIC) | Determinants (baseline) | ESM-2 PLM only | Determinants + PLM | Δ EA |
|---|---|---|---|---|---|
| meropenem | 449 | 52.8% | 56.6% | **60.0%** | +7.2% |
| imipenem | 959 | 53.5% | 57.8% | **59.8%** | +6.3% |

**Reading.** Δ EA = (determinants + PLM) − determinants. The single-run table above has very high fold
variance (±11–21%) because A. baumannii carbapenem MIC is genuinely hard (baseline EA ~53%, far below
K. pneumoniae's 68%) and the labelled sets are smaller. So these deltas must be confirmed with paired CV.

## Verdict (30 paired folds, 6 repeats × 5) — significant on imipenem, borderline on meropenem

| Drug | n | lineages | Determinants | + ESM-2 150M | Δ EA (median) | +PLM wins | Wilcoxon p |
|---|---|---|---|---|---|---|---|
| **imipenem** | 959 | 97 | 53.8% ± 15.8% | **60.5% ± 15.6%** | **+6.7% (+5.6)** | 80% | **0.0001** |
| meropenem | 449 | 67 | 44.8% ± 18.5% | 48.4% ± 16.2% | +3.6% (+4.9) | 67% | 0.051 |

- **Imipenem — a clear, significant win (+6.7 EA points, p=0.0001, wins 80% of folds).** This is the
  strongest carbapenem result across both organisms, and it lands exactly where the biology predicts:
  *A. baumannii* carbapenem resistance is driven by which OXA-carbapenemase variant a strain carries
  (OXA-23 / OXA-24/40 / OXA-58 vs the intrinsic OXA-51-like, each with different hydrolytic strength).
  Presence/absence sees "an OXA gene"; ESM-2 reads the actual variant and separates high- from
  low-MIC alleles.
- **Meropenem — suggestive but not conclusive (p=0.051).** The gain is the same direction and size,
  but with only 449 genomes across 67 lineages the fold variance is too high to cross significance.
  Reported honestly as borderline, not claimed. More labelled meropenem genomes would settle it.

**Cross-organism pattern.** ESM-2 allele resolution helps carbapenem MIC in **both** species —
decisively for *K. pneumoniae* meropenem (+4.8, p=0.0004) and *A. baumannii* imipenem (+6.7, p=0.0001),
suggestively for *A. baumannii* meropenem (+3.6, p=0.051). The effect is specific to carbapenems, whose
MIC magnitude depends on carbapenemase *variant identity* — precisely the signal presence/absence
discards. The pipeline is organism-agnostic: A. baumannii was added as one `ORGANISMS` registry entry,
reusing the shared allele-embedding cache (only 586 of its 648 alleles were new).
