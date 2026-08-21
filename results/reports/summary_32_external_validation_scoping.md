# Summary #32 — External validation, Step 1: scoping the cohort (metadata only)

**Date:** 2026-08-22 · **Status:** scoping complete → **GO**, pending approval for Step 2
**Downloaded:** 179 MB of metadata. **No genomes downloaded. No models run. No numbers produced yet.**

---

## Why this exists

Every genome and phenotype in this project came from a single aggregator, **BV-BRC**. The
lineage-held-out split makes the evaluation honest *within* that corpus, but it cannot show that
performance holds on isolates curated by someone else. That gap is the project's principal
limitation (`REPORT.md` §7, and the N/A evaluation track in Appendix A).

## The trap we had to avoid first

BV-BRC's own documentation states that it collects AMR phenotype data **"from NCBI BioSample and
Antibiogram records as part of the genome ingestion process."** NCBI Pathogen Detection is therefore
**upstream of** our training data. Validating on it naively would mean **testing on isolates we
trained on** — a false external validation, and worse than reporting none at all.

The fix is an accession-level anti-join. This summary is the measurement of whether enough
genuinely-unseen data survives it to be worth downloading genomes for.

## What we pulled

| Item | Source | Size |
|---|---|---|
| Klebsiella AMR metadata (AST phenotypes + accessions) | `ftp.ncbi.nlm.nih.gov/pathogen/Results/Klebsiella/PDG000000012.2502/AMR/PDG000000012.2502.amr.metadata.tsv` | **179 MB** (187,325,514 bytes, verified) |
| BioSample/assembly accessions for our 3,850 training genomes | BV-BRC `genome` API (metadata only) | < 1 MB |

## Result

```
Klebsiella isolates in NCBI PD release PDG000000012.2502 : 172,568
  species_taxid == 573 (K. pneumoniae)                   : 161,779
  ...carrying any AST phenotype                          :   2,529
  ...with >= 1 of our five panel drugs                   :   1,630
```

**AST phenotypes are rare in NCBI PD** — 2,529 of 161,779 K. pneumoniae isolates (1.6%). The
bottleneck is phenotype availability, not genomes.

### Overlap with our training set

| | count | % of 1,630 |
|---|---|---|
| Matched by BioSample accession | 262 | 16.1% |
| Matched by assembly accession | 44 | (subset of the above) |
| **Matched by either — excluded** | **262** | **16.1%** |
| **Genuinely unseen — retained** | **1,368** | **83.9%** |

The 16.1% overlap confirms the ingestion concern was real. It also confirms it is **not fatal**:
five out of six of these isolates never reached our training set.

### The practical cohort

Of the 1,368 unseen isolates, **1,038** carry a downloadable assembly accession (the remainder are
SRA-only and would require assembling from reads — out of scope).

| Drug | R | S | R+S usable | balance |
|---|---|---|---|---|
| meropenem | 314 | 462 | 776 | 0.40 |
| gentamicin | 365 | 359 | 724 | 0.50 |
| ciprofloxacin | 399 | 381 | 780 | 0.49 |
| trimethoprim-sulfamethoxazole | 401 | 346 | 747 | 0.46 |
| cefoxitin | 162 | 334 | 496 | 0.33 |

400 isolates carry R/S calls on **all five** drugs; 760 on at least three. Class balance
(0.33–0.50) is **better than our training data** on several drugs.

Collection years span 2013–2023 with a mode in 2022–2023 — a **later and differently-distributed
sampling frame** than our BV-BRC pull, which strengthens the case that this is a genuine
distribution shift and not a re-slice of the same population.

Cohort manifest written to `data/raw/external_validation/external_cohort_kpneu.csv`.

## Honest caveats — stated now, before any result exists

1. **129 of our 3,850 training genomes (3.4%) carry neither a BioSample nor an assembly accession in
   BV-BRC**, so they could not participate in the anti-join. Worst case — assuming every one of them
   also sits in the NCBI AST set, which is highly unlikely — undetected overlap is bounded above by
   **129 isolates, ≤12.4% of the external cohort**. This bound will be reported alongside the result.
2. **This is "independent-source", not "independent-population" validation.** NCBI Pathogen Detection
   and BV-BRC draw on an overlapping submitter ecosystem. The isolates are new to the model and
   reached us by a different curation path; they are not from a sealed, separately-recruited cohort.
   The write-up must say exactly this and no more.
3. **NCBI does not vet submitted AST methods.** Per NCBI's own AST Browser documentation, "NCBI staff
   do not vet the methods used or values supplied for AST data." Label noise is therefore plausibly
   *higher* here than in our curated training labels — which biases the external result **downward**.
   That is the correct direction for an honesty check, but it must be disclosed.
4. `Intermediate` calls are excluded, consistent with the training convention (still ⚕ open).

## Pre-registered analysis plan (Step 2)

Committed **before** any genome is downloaded, so the result cannot be cherry-picked:

- **Models are frozen.** The existing `results/models/kpneu/*.joblib` are applied as-is. No
  retraining, no re-tuning, no threshold re-selection — the VME≤3% operating points chosen on
  training out-of-fold predictions are carried over unchanged.
- **All five drugs will be reported**, whatever the numbers, with the same metric set used
  internally: ROC-AUC, PR-AUC, VME, ME, and lineage-clustered bootstrap 95% CIs.
- The **rules baseline is scored on the same cohort**, so the "what ML adds" comparison is external too.
- The internal (lineage-held-out) and external numbers will be shown **side by side**, and any
  degradation reported as the headline of that section rather than a footnote.
- If the external result is materially worse, that is the finding and it will be written as such.

## Cost of Step 2 (for approval)

| | estimate |
|---|---|
| Assembly downloads | **~5.1 GB** (1,038 genomes × ~5 MB) |
| AMRFinderPlus annotation | **~13.6 h** at the measured 47 s/genome — resumable, overnight |
| MLST typing (for lineage-clustered CIs) | ~1 h |
| New model training | **none** |

## Recommendation

**GO.** 1,038 unseen isolates with 496–780 usable R/S calls per drug is enough for tight confidence
intervals on all five drugs, and it closes the single largest gap in the submission. The 179 MB spent
here bought the answer to "is there enough genuinely independent data?" without touching a genome.
