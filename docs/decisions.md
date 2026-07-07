# Decision Log (Lab Notebook)

Every non-trivial scientific or engineering choice, one line each, with a rationale. This is the
basis for the write-up and defense. ⚕ = microbiology sign-off required.

| # | Date | Phase | Decision | Rationale | Owner |
|---|---|---|---|---|---|
| 1 | 2026-07-07 | 0 | Scaffold repo per CLAUDE.md §5; brief → `CLAUDE.md`, docx+diagrams → `proposal/`. | Reproducible skeleton before any science. | SWE |
| 2 | 2026-07-07 | 0 | Fixed global seed = 42 in `config/config.yaml`. | Reproducibility mandate. | SWE |
| 3 | 2026-07-07 | 0 | Both AMRFinderPlus + mlst in env; RGI commented as alternative. | Keep env ready; annotator chosen in Phase 4. | SWE |

## Open decisions (pending)
- ⚕ **Phase 1:** organism + drug shortlist (data-driven).
- ⚕ **Phase 3:** MIC/SIR → binary R/S mapping; handling of "intermediate" (config default: drop).
- ⚕ **Phase 5:** lineage method — MLST vs. phylogenetic cluster — and granularity.
- **Phase 4:** AMRFinderPlus vs. RGI as the annotation engine.
- **Phase 3:** genome QC thresholds (config defaults: completeness ≥ 0.90, contamination ≤ 0.05).
- **Phase 1/3:** minimum genomes per class to keep a drug (config default: 100).
