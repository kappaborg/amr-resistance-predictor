# Phase 2 — Data Acquisition Plan (Klebsiella pneumoniae)

**Status:** DRAFT — **download size EXCEEDS the 5 GB ceiling; approval required before any pull.**
Organism/panel locked (decisions #11–13). This phase turns the lab labels into genomes on disk.

## 1. What we acquire
For *K. pneumoniae* (taxon 573), lab phenotypes only (`evidence="Laboratory Method"`) for the panel:
meropenem · gentamicin · ciprofloxacin · trimethoprim/sulfamethoxazole · cefoxitin.
Two artifacts:
1. **Phenotype table** — (genome_id, antibiotic, R/S/I, method) — tiny, pull in full always.
2. **Assembled genome FASTA** — one nucleotide contigs file per genome — this is the heavy part.

## 2. Real footprint (measured, metadata-only)
| Scope | Distinct genomes | ~gzipped | ~uncompressed |
|---|---|---|---|
| Union of all 5 drugs | **6,597** | **~12.2 GB** | ~36.6 GB |
| ciprofloxacin only (thin slice) | 4,809 | ~8.6 GB | ~26 GB |
| meropenem only | 6,074 | ~10.9 GB | ~33 GB |

Per-drug distinct genomes: meropenem 6074 · gentamicin 4890 · ciprofloxacin 4809 · TMP-SMX 3746 ·
cefoxitin 3199. Mean assembly 5.54 Mbp (min 4.96, max 8.35); 96% WGS / 4% Complete (300-genome sample).

**→ Full pull (~12 GB gzipped) is >2× the `data.max_download_gb: 5` ceiling. Decision needed (§5).**

## 3. Phased download strategy (recommended)
Aligns with the sprint's "thin slice first, then broaden":
- **Phase 2a — thin slice (now):** a **balanced ~1,500-genome ciprofloxacin subsample**
  (~750 R / ~750 S), **~2.7 GB gzipped** — under the ceiling. Enough to run the whole Week-1
  pipeline (features → lineage split → model → metrics) and clear the go/no-go.
- **Phase 2b — full set (after the gate):** the remaining genomes for all 5 drugs (~12 GB gzipped),
  pulled only if Week-1 succeeds and you approve the size.
- Subsampling caveat: pre-MLST we can't stratify by lineage, so the thin-slice subsample is **random
  within class**. Final models (Phase 2b/7) use the full set so no lineage is silently dropped.

## 4. Mechanism & reproducibility
- **Phenotypes:** BV-BRC `genome_amr` API → `data/raw/kpneumoniae_amr.csv` (deduped to one R/S per
  genome×drug in Phase 3).
- **Genomes:** BV-BRC genome-sequence API, one gzipped contigs FASTA per genome_id →
  `data/raw/genomes/<genome_id>.fna.gz`. Resumable (skip already-present ids), rate-limited,
  SHA256 per file.
- **Manifest:** every pull logged in `data/manifest.md` — date, query, #genomes, total bytes,
  index checksum. `data/raw/` stays git-ignored.
- **No fabrication / no silent capping:** any subsample size is logged; skipped/failed downloads reported.

## 5. Decision required before pulling
Pick one (see the follow-up question):
- **A. Phased (recommended):** thin-slice cipro ~1.5k now (~2.7 GB) → full after gate.
- **B. Capped full:** balanced cap at ~3,000 genomes across the panel (~5.5 GB) — one pull, near ceiling.
- **C. Full now:** all 6,597 (~12 GB gzipped) — raise the ceiling explicitly.

## 6. Alternative worth knowing (not recommended as primary)
BV-BRC ships **precomputed** AMR gene annotations (`sp_gene` core) — using them would avoid the genome
download entirely, but the proposal calls for running **AMRFinderPlus ourselves** (reproducibility +
control over the determinant vocabulary). Recommendation: run our own annotation; optionally use
`sp_gene` later only as a cross-check. Flagged so the trade-off is on the record.

## 7. Compute note (for Phase 4, not now)
AMRFinderPlus ≈ 15–40 s/genome. Thin slice (~1.5k) ≈ 1–3 h multicore; full (~6.6k) ≈ 6–12 h multicore.
Laptop-feasible; the genome download and annotation are the real time costs, not modeling.
