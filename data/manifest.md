# Data Manifest

Every data source used in this project is recorded here: source, version/date, exact query, and
checksum. `data/raw/` is git-ignored and regenerates from this manifest. **No entry, no data.**

## Sources
| Role | Source | Version / date | Access | Status |
|---|---|---|---|---|
| Genomes + phenotypes | BV-BRC (formerly PATRIC) | TBD | CLI/API | not yet pulled |
| Reference determinants | CARD | TBD | download | not yet pulled |
| Annotation DB | AMRFinderPlus Reference Gene Catalog | **2026-05-15.1** | `amrfinder -u` | ✅ installed |
| Lineage typing | `mlst` schemes (PubMLST) | bundled w/ mlst 2.33.1 | bundled w/ tool | ✅ available |

## Acquisition log
_(one row per pull)_

| Date | Organism | Scope | #Genomes | Size (gz) | Selection |
|---|---|---|---|---|---|
| 2026-07-07 | K. pneumoniae (573) | thin-slice ciprofloxacin | 1499 (750 R / 749 S) | 2.49 GB | balanced, seed 42; 1 genome (573.46086) failed download |

**Query (phenotypes):** `genome_amr` · `taxon_id=573` · `evidence="Laboratory Method"` ·
`antibiotic=ciprofloxacin` → 3728 consistent-label genomes (2548 R / 1180 S) → balanced sample.
**Query (genomes):** `genome_sequence` · `eq(genome_id,X)&limit(25000)` · `application/dna+fasta`.
**Files:** `data/raw/genomes/<id>.fna.gz` · labels `data/raw/thin_slice_cipro_labels.csv` ·
checksums `data/raw/thin_slice_cipro_checksums.csv`. **Integrity:** verified PASS
(mean 5.64 Mbp, 0 checksum/gzip/FASTA errors). Two genomes flagged for Phase-3 QC
(one 11.5 Mbp = possible contamination; one 3848 contigs = fragmented).

## Tool versions (env `amr-resistance-predictor`, built 2026-07-07, osx-arm64, conda 24.1.2)
| Tool / lib | Version |
|---|---|
| python | 3.11.15 |
| numpy / pandas / scikit-learn | 2.4.6 / 3.0.3 / 1.9.0 |
| xgboost / lightgbm / shap | 3.2.0 / 4.6.0 / 0.51.0 |
| matplotlib / requests / pytest | 3.11.0 / 2.34.2 / 9.1.1 |
| ncbi-amrfinderplus | 4.2.7 |
| mlst | 2.33.1 |

## Notes
- Phase 1 records **counts only** (metadata queries) — no bulk download until the data plan is approved.
- Flag any download > `data.max_download_gb` (config) before pulling.
- ⚠️ **AMRFinderPlus reference database not yet downloaded.** Run `amrfinder -u` in Phase 4 before
  annotation (network download — flag/record size when pulled). CARD data also pulled in Phase 4.
- `mlst` bundles its PubMLST schemes with the package build.
