# Data Manifest

Every data source used in this project: source, version/date, exact query, and provenance.
`data/raw/` is git-ignored and regenerates from this manifest + code. **No entry, no data.**

## Sources
| Role | Source | Version / date | Access | Status |
|---|---|---|---|---|
| Genomes + phenotypes | BV-BRC (formerly PATRIC) | pulled 2026-07-07 … 2026-07-09 | REST Data API (`requests`) | ✅ pulled |
| Reference determinants | CARD (via AMRFinderPlus catalog) | catalog 2026-05-15.1 | bundled in AMRFinderPlus DB | ✅ |
| Annotation DB | AMRFinderPlus Reference Gene Catalog | **2026-05-15.1** | `amrfinder -u` | ✅ installed |
| Lineage typing | `mlst` schemes (PubMLST) | bundled w/ mlst 2.33.1 | bundled w/ tool | ✅ |
| Protein language model | ESM-2 (`esm2_t30_150M_UR50D`, `esm2_t33_650M_UR50D`) | fair-esm 2.0.0 weights | auto-download (fair-esm) | ✅ cached |

## Acquisition log — genomes + phenotypes (one row per organism)
All via BV-BRC. Phenotype filter is **laboratory-confirmed only** (`evidence="Laboratory Method"`);
genomes are assembled contigs (`genome_sequence`, `limit(25000)`, `application/dna+fasta`). QC per
`docs/decisions.md` #19 (completeness ≥90%, contamination ≤5%, contigs ≤500, organism-specific length).

| Organism | Taxon | AMRFinderPlus `--organism` | MLST scheme | Genomes (post-QC) | Determinant features | Panel drugs |
|---|---|---|---|---|---|---|
| *Klebsiella pneumoniae* | 573 | Klebsiella_pneumoniae | klebsiella | 3,850 | 688 | meropenem, gentamicin, ciprofloxacin, TMP-SMX, cefoxitin |
| *Escherichia coli* | 562 | **Escherichia** | ecoli | 3,035 | 595 | ciprofloxacin, gentamicin, TMP-SMX, ampicillin, ceftriaxone |
| *Acinetobacter baumannii* | 470 | Acinetobacter_baumannii | abaumannii | 1,395 | 267 | meropenem, imipenem, ciprofloxacin, gentamicin, amikacin |
| *Staphylococcus aureus* | 1280 | Staphylococcus_aureus | saureus | 3,532 | 149 | oxacillin, cefoxitin, ciprofloxacin, erythromycin, clindamycin |
| *Pseudomonas aeruginosa* | 287 | Pseudomonas_aeruginosa | paeruginosa | 1,057 | 798 | meropenem, ceftazidime, ciprofloxacin, tobramycin |
| *Salmonella enterica* | 28901 | Salmonella | salmonella | 2,321 | 159 | ampicillin, ceftriaxone, ciprofloxacin, chloramphenicol, TMP-SMX |
| *Enterococcus faecium* | 1352 | Enterococcus_faecium | efaecium | 1,890 | 99 | vancomycin, ampicillin, tetracycline |
| *Streptococcus pneumoniae* | 1313 | Streptococcus_pneumoniae | spneumoniae | 3,857 | 20 | penicillin, erythromycin, tetracycline, TMP-SMX |

**Evaluated and excluded (data-limited):** *Enterobacter cloacae* (taxon 550/complex 354276 — ≤109
genomes/drug) and *Campylobacter jejuni* (taxon 197 — ~489 genomes, only tetracycline near-viable).
Documented in `docs/decisions.md` #46/#49 — honest data-driven organism selection.

**Phenotype query (per drug):** `genome_amr` · `and(eq(taxon_id,T),eq(evidence,"Laboratory Method"),`
`eq(antibiotic,DRUG))` · `select(genome_id,resistant_phenotype)` · `limit(25000)`; consistent-label
genomes only (a genome with conflicting R/S records is dropped).
**MIC query (regression):** same filter with `select(genome_id,measurement,measurement_unit)`; keep
`measurement_unit="mg/L"`, strip censor signs (`<`,`>`,`=`), take median log2(MIC) per genome.
**Files:** `data/raw/genomes/<genome_id>.fna.gz` · features `data/processed/<tag>_features.csv` ·
lineages `data/processed/<tag>_lineages.csv` · K. pneumoniae labels `data/processed/panel_labels.csv`
(other organisms fetched live). Tags: `thin_slice_cipro` (K.pneu), `ecoli`, `abaumannii`, `saureus`.

## Derived data — ESM-2 allele embeddings (`data/interim/esm/`)
Resistance-gene proteins extracted from each genome via AMRFinderPlus coordinates → BioPython
translate (`src/models/esm2_mic.py`), deduplicated by sequence, embedded once with ESM-2 (Apple MPS).
| Artifact | Contents |
|---|---|
| `proteins.pkl` | K. pneumoniae: 3,850 genomes → 58,502 AMR proteins → **1,945 unique alleles** |
| `proteins_abaum.pkl` | A. baumannii: 1,395 genomes → 17,052 proteins → 648 unique (586 new) |
| `embeddings_150M.pkl` | 2,531 unique-allele vectors (640-d), ESM-2 150M |
| `embeddings_650M.pkl` | 1,945 unique-allele vectors (1280-d), ESM-2 650M |

## Tool / library versions (env `amr-resistance-predictor`, osx-arm64 / Apple M1 Max, MPS backend)
| Tool / lib | Version | | Tool / lib | Version |
|---|---|---|---|---|
| python | 3.11 | | ncbi-amrfinderplus | 4.2.7 (DB 2026-05-15.1) |
| numpy / pandas | 2.x | | mlst | 2.33.1 |
| scikit-learn / xgboost | 1.x / 2.x | | torch | 2.13.0 (MPS) |
| scipy / shap | 1.17.1 / 0.5x | | fair-esm | 2.0.0 |
| matplotlib / requests | 3.x / 2.34.2 | | biopython | 1.87 |

## Notes
- Phenotypes are **laboratory-confirmed only** — computational predictions in BV-BRC are excluded.
- Splits are **lineage-grouped by MLST sequence type** (zero train/test lineage overlap; `tests/test_split.py`).
- AMRFinderPlus reference DB pulled via `amrfinder -u` (the CARD/NCBI determinant catalog, 2026-05-15.1).
- ESM-2 weights auto-download on first use and cache under the fair-esm torch hub directory.
- No genome, label, or MIC value is synthesized; missing data is dropped and reported (never imputed).
