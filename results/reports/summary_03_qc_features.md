# Summary #3 — Phase 3 (QC filter) complete · Phase 4 (features) started

**Date:** 2026-07-07 · **Status:** Phase 3 ✅ · Phase 4 ⏸ (compute decision — see below).

## Phase 3 — QC & label harmonization ✅
- Applied the ⚕-approved, literature-aligned QC filter (`src/data/qc_filter.py`):
  completeness ≥90%, contamination ≤5% (MIMAG), contigs ≤500, length 4.5–7.5 Mbp (K. pneumoniae).
- **Kept 1,472 genomes for modeling**, dropped 27 (failed QC) + 1 (download failure).
- **Class balance 735 R / 737 S** (minority fraction 0.50 — essentially perfect).
- "Intermediate" excluded (approved). Output: `data/processed/thin_slice_cipro_labels_qc.csv`.

## Phase 4 — feature extraction (in progress)
- Built `src/features/build.py`: AMRFinderPlus per genome with `--organism Klebsiella_pneumoniae`
  (enables gyrA/parC point mutations — the key cipro determinants), resumable TSV cache, parallel,
  then a binary genome × determinant matrix (Type==AMR: acquired genes + point mutations).
- **Verified on 8 genomes:** 76 determinants, correct signal — `parC_S80I` in 62%, oqxA/oqxB efflux,
  blaCTX-M-15/SHV, sul1, aadA, mph(A). AMRFinderPlus DB 2026-05-15.1.

## ⚠️ Compute flag (why Phase 4 is paused)
- Measured throughput: **~47 s/genome** at 8 workers on 10 cores.
- **Full 1,472 genomes ≈ 15–19 h of CPU-pegged runtime** on the laptop.
- Options: (A) full run overnight (resumable background); (B) balanced ~600-genome subset first
  (~7 h) to get Week-1 numbers fast, annotate the rest later; (C) offload to Colab/other machine.
- Job is **resumable** (skips cached TSVs), so an interrupted run continues without rework.

## Next
Await compute choice → run annotation → build matrix → Phase 5 (lineage split) → Phase 6 (thin-slice
model + first per-drug numbers, the Week-1 go/no-go).
