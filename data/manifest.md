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

## Notes
- Phase 1 records **counts only** (metadata queries) — no bulk download until the data plan is approved.
- Flag any download > `data.max_download_gb` (config) before pulling.
- Pin exact tool + database versions here once the environment is built (Phase 0/2).
