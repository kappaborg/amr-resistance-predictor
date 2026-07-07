# Decision Log (Lab Notebook)

Every non-trivial scientific or engineering choice, one line each, with a rationale. This is the
basis for the write-up and defense. ⚕ = microbiology sign-off required.

| # | Date | Phase | Decision | Rationale | Owner |
|---|---|---|---|---|---|
| 1 | 2026-07-07 | 0 | Scaffold repo per CLAUDE.md §5; brief → `CLAUDE.md`, docx+diagrams → `proposal/`. | Reproducible skeleton before any science. | SWE |
| 2 | 2026-07-07 | 0 | Fixed global seed = 42 in `config/config.yaml`. | Reproducibility mandate. | SWE |
| 3 | 2026-07-07 | 0 | Both AMRFinderPlus + mlst in env; RGI commented as alternative. | Keep env ready; annotator chosen in Phase 4. | SWE |
| 4 | 2026-07-07 | 0 | Dropped `bvbrc-cli` pip dep; use BV-BRC REST Data API via `requests`. | No maintained pip CLI package; API is sufficient and lighter. | SWE |
| 5 | 2026-07-07 | 0 | Env built natively on osx-arm64 (no Rosetta); bio tools have arm64 builds. | Confirmed amrfinderplus 4.2.7 + mlst 2.33.1 available for osx-arm64. | SWE |
| 6 | 2026-07-07 | 1 | Survey via BV-BRC `genome_amr`; counts by antibiotic × phenotype using facet + Content-Range. | Exact counts with no bulk download; metadata only. | SWE |
| 7 | 2026-07-07 | 1 | ~~min(R,S)≥100, balance≥0.20~~ → **superseded by #9**. | Initial guess before data. | SWE |
| 8 | 2026-07-07 | 1 | Count **lab phenotypes only** (`evidence="Laboratory Method"`); exclude computational predictions. | Using BV-BRC computational AMR calls as labels would be circular. E.coli: 6.3M rows → 243K lab. | SWE |
| 9 | 2026-07-07 | 1 | **Tuned thresholds from live data: min(R,S)≥150, balance≥0.30, total(R+S)≥300.** | Data abundant (median min(R,S)=191); balance is the binding quality lever (median 0.15). Yields 51 balanced drugs / 7 orgs. | SWE (⚕ review) |
| 10 | 2026-07-07 | 1 | genome_amr has no `taxon_lineage_ids`; filter by species `taxon_id`. | Verified against live API schema. | SWE |
| 11 | 2026-07-07 | 1 | **LOCKED: organism = _Klebsiella pneumoniae_ (taxon 573).** | Best fit to thesis: carbapenem resistance is where single-gene rules fail (porin+β-lactamase combos → catchable VMEs); strong MLST clonality makes the phylogeny-aware split bite. | SWE + ⚕ |
| 12 | 2026-07-07 | 1 | **LOCKED: 5-drug panel** = meropenem, gentamicin, ciprofloxacin, TMP-SMX, cefoxitin. | 5 distinct determinant families → convincing per-drug SHAP; all clear thresholds. | SWE + ⚕ |
| 13 | 2026-07-07 | 1 | **Thin-slice drug = ciprofloxacin** (meropenem = Week-2 flagship). | Chromosomal gyrA/parC point mutations: clean direct signal, lineage-spread, exercises point-mutation feature path → de-risks Week-1 go/no-go. | SWE + ⚕ |

## Open decisions (pending)
- ⚕ **Phase 1:** organism + drug shortlist (data-driven).
- ⚕ **Phase 3:** MIC/SIR → binary R/S mapping; handling of "intermediate" (config default: drop).
- ⚕ **Phase 5:** lineage method — MLST vs. phylogenetic cluster — and granularity.
- **Phase 4:** AMRFinderPlus vs. RGI as the annotation engine.
- **Phase 3:** genome QC thresholds (config defaults: completeness ≥ 0.90, contamination ≤ 0.05).
- **Phase 1/3:** minimum genomes per class to keep a drug (config default: 100).
