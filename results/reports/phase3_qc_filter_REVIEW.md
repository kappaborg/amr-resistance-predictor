# Genome Quality-Control Filter — For Microbiology Review & Sign-off

**Project:** Reading Resistance (K. pneumoniae AMR predictor) · **Phase:** 3 (QC & label harmonization)
**Prepared for:** co-builder (M.Sc. Microbiology) · **Date:** 2026-07-07
**Decision requested:** approve the QC thresholds below (or adjust) so we can filter the genomes and
proceed to feature extraction. **Nothing is deleted** — filtered genomes are excluded from modeling
but kept on disk.

---

## 1. Why we filter at all
We predict resistance from the **genes and mutations present in each genome**. If an assembly is
**contaminated** (contains DNA from a second organism) or **incomplete/fragmented**, the gene list is
wrong — a contaminant can add resistance genes the strain doesn't really have, and a broken assembly
can miss genes that are there. Either way the model learns from corrupted features. QC removes the
genomes we can't trust *before* they reach the model. This is standard practice and protects every
downstream number.

## 2. What we measured (per genome, from BV-BRC / CheckM)
| Metric | Plain meaning | What bad looks like |
|---|---|---|
| **Completeness %** | How much of the expected single-copy gene set is present | Low → missing genes → false "gene absent" |
| **Contamination %** | Fraction of duplicated/foreign marker genes → mixed or contaminated DNA | High → foreign resistance genes → false "gene present" |
| **Contigs** | Number of assembly fragments | Very high → fragmented → genes split/missed |
| **Genome length** | Total assembled bp | Far from the species norm (~5.3 Mbp) → contamination or breakage |

## 3. What our data looks like (thin slice: 1,500 K. pneumoniae genomes, ciprofloxacin)
The dataset is **high quality overall** — BV-BRC labels **1,485 / 1,500 (99%) as "Good."**
- Completeness: median **100%**, 5th percentile **99.7%** (only 3 genomes below 90%).
- Contamination: median **0%**, 95th percentile **1.0%** (14 genomes above 5%, one as high as 51%).
- Contigs: median **92**, 95th percentile **309**, worst case **3,848**.
- Genome length: median **5.50 Mbp**; a handful inflated to 8–11 Mbp (a contamination signature).

## 4. Proposed thresholds (the "moderate" tier)
| Criterion | Threshold | Why this value |
|---|---|---|
| Completeness | **≥ 95%** | Our genomes sit at ~100%; 95% still admits all real assemblies, cuts only broken ones |
| Contamination | **≤ 5%** | Standard "acceptable" ceiling; our clean genomes are ~0%, so this only catches true mixtures |
| Contigs | **≤ 500** | Removes severely fragmented assemblies (median is 92) |
| Genome length | **4.8 – 6.5 Mbp** | Species sanity window around K. pneumoniae ~5.3 Mbp; catches contamination/breakage |

## 5. Result of applying the filter
- **Kept: 1,469 genomes · Dropped: 31 (2%).**
- **The drop is class-balanced — 16 Resistant / 15 Susceptible** — so filtering does **not** bias the
  resistant/susceptible ratio. This matters: a filter that removed mostly one class would distort results.
- Reasons genomes were dropped (some hit more than one): contamination > 5% → 14; contigs > 500 → 23;
  length outside window → 11; completeness < 95% → 5.

## 6. The genomes being dropped — please spot-check
The worst offenders are clearly **contaminated/mixed assemblies** (high contamination *and* inflated
length — a classic signature of two genomes merged). All happen to be labeled Resistant, which is
exactly why they're dangerous to keep: a contaminant's resistance genes would teach the model the
wrong lesson.

| genome_id | completeness % | contamination % | contigs | length Mbp | label |
|---|---|---|---|---|---|
| 573.31080 | 99.4 | **51.3** | 1715 | **11.15** | Resistant |
| 573.12987 | 100 | **48.7** | 3848 | 8.35 | Resistant |
| 573.13822 | 100 | 40.7 | 1747 | 8.61 | Resistant |
| 573.13345 | 100 | 40.3 | 1209 | 8.86 | Resistant |
| 573.15492 | 100 | 35.9 | 276 | 9.47 | Resistant |
| 573.14450 | 100 | 27.1 | 1367 | 7.94 | Resistant |

(Look these up in BV-BRC by genome_id if you'd like to confirm they're genuinely contaminated.)

## 7. Why we are NOT filtering more strictly
A stricter filter (e.g. completeness ≥ 98%, length 5.0–6.2) would drop ~98 genomes. We deliberately
avoid that: **every genome removed is a lineage we can no longer test on**, and the project's central
claim is generalizing to *unseen bacterial lineages* (the phylogeny-aware split). Removing clean
genomes just to look tidy would shrink lineage diversity and make the honest evaluation weaker, not
stronger. We remove genuine junk, not borderline-good genomes.

## 8. Decisions we need from you (microbiology sign-off)
1. **Contamination ceiling — 5%?** (Standard; our data is ~0% so this only removes true mixtures.)
2. **Length window — 4.8–6.5 Mbp** for K. pneumoniae? (Adjust if you'd prefer a tighter/wider range.)
3. **Completeness ≥ 95% and contigs ≤ 500** — acceptable?
4. **"Intermediate" phenotype handling** — currently the labels use only Resistant/Susceptible and
   exclude Intermediate. Confirm this is the right call for these drugs, or tell us how to treat it.

## 9. Sign-off
- [ ] Thresholds approved as-is
- [ ] Thresholds approved with changes: ______________________________________________
- [ ] "Intermediate" handling confirmed: exclude / treat as ____________
- Reviewer: ____________________  Date: __________

Once approved, we apply the filter, finalize the R/S labels, and move to feature extraction
(AMRFinderPlus) — the next step in the pipeline.
