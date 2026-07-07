# Summary #0 — Phase 0: Project Setup & Scaffold

**Date:** 2026-07-07 · **Status:** ✅ Complete · **Gate:** none (internal prep)

## What was done
- Reorganized loose files per CLAUDE.md §5: brief → `CLAUDE.md` at root; proposal docx + 3 diagrams
  → `proposal/` (original brief preserved as `proposal/CLAUDE_AMR_original.md`).
- Built the full directory tree: `src/{data,features,split,models,evaluation,interpret,app}`,
  `config/`, `data/{raw,interim,processed}`, `notebooks/`, `tests/`,
  `results/{figures,metrics,models,reports}`, `docs/`. Python package markers added.
- Authored config & build files: `environment.yml`, `config/config.yaml`, `Makefile`,
  `README.md`, `.gitignore`, `data/manifest.md`, `docs/decisions.md`.
- Wrote the **leakage test** (`tests/test_split.py`) — the project's core guardrail — running now
  on synthetic fixtures so the no-shared-lineage invariant is locked before real data exists.
- Initialized git; committed 27 tracked files (caches/raw data correctly ignored).

## Verification
- `pytest`: **2 passed, 1 skipped** (skipped = real-split test, enabled in Phase 5).
- Git clean: no `__pycache__`, `.pytest_cache`, or `.pyc` tracked.

## Decisions logged (docs/decisions.md)
1. Scaffold per CLAUDE.md §5. · 2. Global seed = 42. · 3. AMRFinderPlus + mlst in env; RGI as commented alternative.

## Key config defaults set (all overridable, several ⚕ pending)
- `min_genomes_per_class: 100`, `max_download_gb: 5` (hard stop-and-flag ceiling)
- QC: completeness ≥ 0.90, contamination ≤ 0.05
- split: `mlst`, 25% of *lineages* held out · model: LogReg baseline → XGBoost → isotonic calibration

## Environment — BUILT & VERIFIED ✅
- Created env `amr-resistance-predictor` (conda 24.1.2, `/opt/homebrew/anaconda3`), **native osx-arm64**
  (no Rosetta). Python 3.11.15.
- Verified imports: numpy 2.4.6 · pandas 3.0.3 · scikit-learn 1.9.0 · xgboost 3.2.0 · lightgbm 4.6.0 ·
  shap 0.51.0 · matplotlib 3.11.0 · requests 2.34.2 · pytest 9.1.1.
- Bio tools run: **AMRFinderPlus 4.2.7**, **mlst 2.33.1**. Versions pinned in `data/manifest.md`.
- `pytest` inside the env: **2 passed, 1 skipped**.
- Config change: dropped non-existent `bvbrc-cli` pip dep → use BV-BRC REST API via `requests`.

## Blocked / caveats
- ⚠️ **AMRFinderPlus reference DB not downloaded yet** — run `amrfinder -u` in Phase 4 (network
  download; record size). Not needed until feature extraction.

## Next → Phase 1 (GATE, ⚕ microbiology sign-off)
BV-BRC organism/drug **selection**: write the query plan, run metadata-only counts, present a
class-balanced shortlist. **No downloads until the data plan is approved.**
