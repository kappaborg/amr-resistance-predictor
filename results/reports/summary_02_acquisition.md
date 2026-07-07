# Summary #2 — Phase 2: Data Acquisition (thin slice)

**Date:** 2026-07-07 · **Status:** ✅ Complete · Integrity **PASS**.

## What was done
- Locked target (decisions #11–13): **K. pneumoniae**, 5-drug panel, thin-slice drug **ciprofloxacin**.
- Pulled cipro lab phenotypes: 3,728 genomes with consistent R/S (2,548 R / 1,180 S, 0 conflicts).
- Downloaded a **balanced 750 R / 750 S** genome subsample (seed 42) as gzipped contigs FASTA.
- Installed the **AMRFinderPlus reference DB 2026-05-15.1** (for Phase 4).
- Built a reusable integrity checker (`src/data/verify_download.py`).

## Result
- **1,499 genomes, 2.49 GB gzipped** (1 genome `573.46086` failed download after retries → 750 R / 749 S).
- Integrity **PASS**: mean assembly **5.64 Mbp**, 0 checksum mismatches, 0 gzip failures, 0 non-FASTA,
  0 undersized, every genome labeled + checksummed.

## Two bugs caught by verification (both fixed) — the value of checking
1. **Silent API truncation.** `genome_sequence` defaults to `limit(25)` → assemblies cut to ~2 Mbp
   (25 contigs). Would have dropped genes on missing contigs → false-absent features → corrupted
   labels. **Fix:** `limit(25000)` + a <2.5 MB reject guard. All genomes re-downloaded.
2. **No timeout resilience.** A single network `ReadTimeoutError` aborted the run at genome 306.
   **Fix:** 4-attempt retry with backoff; run is resumable (skips existing files).

## Flagged for Phase 3 QC (not errors — expected)
- One genome 11.5 Mbp (K. pneumoniae is ~5.3 Mbp → possible contamination/mixed assembly).
- One genome 3,848 contigs (highly fragmented). The completeness/contamination filter handles both.

## Artifacts
- `data/raw/genomes/*.fna.gz` (1,499) · `data/raw/thin_slice_cipro_labels.csv` ·
  `data/raw/thin_slice_cipro_checksums.csv` · manifest updated.

## Next → Phase 3 (QC & label harmonization) then Phase 4 (feature extraction)
- QC: apply completeness/contamination + contig thresholds; drop the flagged outliers.
- Then run AMRFinderPlus over the QC'd set → determinant feature matrix (Phase 4).
- Note: the full 5-drug panel (~6.6k genomes, ~12 GB) is deferred to Phase 2b **after the Week-1 gate**.
