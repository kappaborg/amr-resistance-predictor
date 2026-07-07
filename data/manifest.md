# Data Manifest

Every data source used in this project is recorded here: source, version/date, exact query, and
checksum. `data/raw/` is git-ignored and regenerates from this manifest. **No entry, no data.**

## Sources
| Role | Source | Version / date | Access | Status |
|---|---|---|---|---|
| Genomes + phenotypes | BV-BRC (formerly PATRIC) | TBD | CLI/API | not yet pulled |
| Reference determinants | CARD | TBD | download | not yet pulled |
| Annotation DB | AMRFinderPlus Reference Gene Catalog | TBD | `amrfinder -u` | not yet pulled |
| Lineage typing | `mlst` schemes (PubMLST) | TBD | bundled w/ tool | not yet pulled |

## Acquisition log
_(one row per pull — filled in Phase 2)_

| Date | Organism | Query | #Genomes | Size | SHA256 (index) |
|---|---|---|---|---|---|
| — | — | — | — | — | — |

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
