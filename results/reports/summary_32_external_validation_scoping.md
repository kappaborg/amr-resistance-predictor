# Summary #32 — External validation, Step 1: scoping the cohort (metadata only)

**Date:** 2026-08-22 · **Status:** scoping complete → **GO**, pending approval for Step 2
**Downloaded:** ~196 MB of metadata. **No genomes downloaded. No models run. No numbers produced yet.**

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

The fix is a two-layer exclusion: a **provenance filter** where the source exposes one, plus an
**accession-level anti-join** in every case.

## Sources scoped

| # | Source | Downloaded | Size |
|---|---|---|---|
| 1 | **EMBL-EBI AMR Portal / CABBAGE**, release `2026-07` | `phenotype.csv.gz` (md5-verified) | **16.6 MB** |
| 2 | **NCBI Pathogen Detection**, Klebsiella release `PDG000000012.2502` | `*.amr.metadata.tsv` | **179 MB** |
| — | BV-BRC accessions for our 3,850 training genomes | `genome` API, metadata only | < 1 MB |

Source 1 (CABBAGE; Dickens et al., bioRxiv 2025, doi:10.1101/2025.11.12.688105, **CC BY 4.0**, no
login) is the better of the two because it carries an explicit **`database` provenance column**. For
*K. pneumoniae* its 133,136 phenotype rows break down as PATRIC 55,263, PATRIC;CABBAGE_PubMed 40,626,
NCBI_antibiogram 13,590, CABBAGE_PubMed 7,004, and so on — so **every row PATRIC/BV-BRC touched can be
dropped by construction**, not merely by hoping the accessions don't collide.

## Result

Both sources were filtered (Source 1: drop any row whose provenance contains `PATRIC`) and then
anti-joined against our training accessions.

| | isolates |
|---|---|
| Source 1 — AMR Portal, PATRIC-filtered + anti-joined, with assembly | **800** |
| Source 2 — NCBI PD, anti-joined, with assembly | **1,036** |
| intersection | 524 |
| Portal-only | 276 |
| NCBI-only | 512 |
| **UNION** | **1,312** |
| ...with ≥1 usable R/S call on our panel | **1,144** |

For Source 2, of 1,630 NCBI PD *K. pneumoniae* isolates carrying ≥1 panel drug, **262 (16.1%)
overlapped our training set and were excluded**. AST is the bottleneck, not genomes: only 2,529 of
161,779 NCBI PD *K. pneumoniae* (1.6%) carry any AST at all.

### Merged cohort, per drug

| Drug | R | S | R+S usable | balance |
|---|---|---|---|---|
| ciprofloxacin | 451 | 493 | 944 | 0.48 |
| gentamicin | 412 | 478 | 890 | 0.46 |
| meropenem | 339 | 514 | 853 | 0.40 |
| trimethoprim-sulfamethoxazole | 423 | 415 | 838 | 0.50 |
| cefoxitin | 167 | 339 | 506 | 0.33 |

406 isolates carry R/S on **all five** drugs. Balance (0.33–0.50) is **better than our training data**
on several drugs. Collection years span 2013–2023 with a mode in 2022–23 — a later and
differently-distributed sampling frame than our July 2026 BV-BRC pull.

### A free consistency check

The two sources overlap on 524 isolates, curated by different paths. On our five drugs they produced
**zero conflicting R/S calls**. That is real evidence the phenotype labels are stable across curation
pipelines, and it is worth reporting.

## ⚠️ The methodological gap this scoping exposed

**Accession-level de-duplication is not sufficient.** Removing isolates that share a BioSample with
our training set does *not* remove **near-clonal siblings** — an external ST258 isolate is a different
BioSample but effectively the same lineage the model already learned. This is precisely the failure
mode of non-negotiable #1, reappearing at the external-validation boundary, and neither source carries
an ST column.

**Step 2 must therefore run `mlst` on the external assemblies and report two numbers:**

1. performance on **all** external isolates, and
2. performance restricted to **sequence types absent from training**.

The gap between those two is the honest headline. Reporting only (1) would repeat, on external data,
the exact leakage this project exists to avoid.

## Honest caveats — stated now, before any result exists

1. **129 of our 3,850 training genomes (3.4%) carry neither a BioSample nor an assembly accession** in
   BV-BRC, so they could not participate in the anti-join. Worst case — assuming every one also sits
   in the external set, which is highly unlikely — undetected overlap is bounded above by **129
   isolates, ≤11.3% of the 1,144-isolate cohort**. This bound will be reported with the result.
2. **This is "independent-source", not "independent-population" validation.** The Portal's provenance
   filter makes Source 1 genuinely independent of PATRIC curation, but the submitter ecosystems still
   overlap. The write-up must say exactly this and no more.
3. **NCBI does not vet submitted AST methods** ("NCBI staff do not vet the methods used or values
   supplied for AST data"). Label noise is plausibly higher than in our curated training labels, which
   biases the external result **downward** — the right direction for an honesty check, but disclosable.
4. **Features must be recomputed, not reused.** Both sources ship precomputed AMRFinderPlus calls, and
   skipping the download is tempting — but NCBI rows carry refgene DB `2026-01-21.1` while our models
   were trained on **`2026-05-15.1`**, and the Portal does not disclose its DB version at all. Only 266
   of our 688 K. pneumoniae feature columns appear in the Portal's vocabulary. A feature-extraction
   mismatch would depress external performance **for reasons unrelated to generalization, and we could
   not tell the two apart.** Step 2 therefore re-runs AMRFinderPlus 4.2.7 / DB 2026-05-15.1 ourselves.
5. `Intermediate` calls are excluded, consistent with the training convention (still ⚕ open).

## Sources that would be FALSE external validation — checked and rejected

Useful defense material. Each is already inside BV-BRC with laboratory-evidence phenotypes:

| Candidate | Why rejected |
|---|---|
| David 2019 **EuSCAPE** *K. pneumoniae* | 1,717/1,717 in BV-BRC |
| **Nguyen 2018/2019** (*K. pneumoniae* / *Salmonella*) | in BV-BRC; the papers state PATRIC genome IDs are published |
| **Moradigaravand 2018** *E. coli* (BSAC) | 1,906/1,936 with lab phenotypes in BV-BRC |
| Gordon 2014 / Bradley 2015 *S. aureus* | 885 + 100, all with lab phenotypes |
| **Pathogenwatch** | AMR is genotypically *predicted*, not lab-measured — would compare our model to another model |
| hzi-bifo `AMR_benchmarking` | PATRIC-derived and PATRIC-keyed; use for method comparison only |
| AllTheBacteria / 661k | no phenotypes at all |
| **CRyPTIC** | genuinely independent and gold-standard, but *M. tuberculosis* — not one of our eight organisms. Cite as the standard we emulate; do not use |
| CDC & FDA AR Isolate Bank | genuinely independent phenotypes (CDC reference broth microdilution) but no bulk download and low yield — much of its Klebsiella already sits in Source 2 |

## Pre-registered analysis plan (Step 2)

Committed **before** any genome is downloaded, so the result cannot be cherry-picked:

- **Models are frozen.** `results/models/kpneu/*.joblib` applied as-is. No retraining, no re-tuning,
  no threshold re-selection — the VME≤3% operating points chosen on training out-of-fold predictions
  carry over unchanged.
- **Features recomputed** with AMRFinderPlus 4.2.7 / DB 2026-05-15.1, matching training exactly.
- **All five drugs reported**, whatever the numbers, with ROC-AUC, PR-AUC, VME, ME and
  lineage-clustered bootstrap 95% CIs.
- **Two cohorts reported side by side:** all external isolates, and external isolates whose **ST is
  absent from training**.
- The **rules baseline is scored on the same cohort**, so "what ML adds" is measured externally too.
- Internal (lineage-held-out) vs external shown side by side; any degradation is the headline of that
  section, not a footnote.
- If the external result is materially worse, that is the finding and it will be written as such.

## Cost of Step 2 (for approval)

| | estimate |
|---|---|
| Assembly downloads | **~1.8 GB** (1,144 × **1.61 MB**, measured on a real NCBI datasets fetch) |
| AMRFinderPlus annotation | **~15 h** at the measured 47 s/genome — resumable, overnight |
| `mlst` typing (required, see above) | ~1 h |
| New model training | **none** |

## Recommendation

**GO.** 1,144 unseen isolates with 506–944 usable R/S calls per drug supports tight confidence
intervals on all five drugs, and closes the single largest gap in the submission. ~196 MB of metadata
bought the answer to "is there enough genuinely independent data, and which source is actually
independent?" without touching a genome.
