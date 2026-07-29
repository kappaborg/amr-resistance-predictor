# Round-4 Expansion — Adding 6 Organisms (scaffolded 2026-07-10)

Goal: extend the pipeline from 4 to **10 WHO-priority pathogens**. Config is fully scaffolded; data
acquisition is **pending BV-BRC availability** (503 outage at scaffold time) and runs **one organism
at a time**, ~3,000+ genomes each, starting with *P. aeruginosa*. Each drug panel below is
**PROVISIONAL** — validate R/S balance against live BV-BRC counts (and ⚕ microbiology sign-off) before
training; drop any drug with too few/too-imbalanced labels (same rule as Rounds 1–3).

## Candidates (tooling verified locally)

| Organism | Taxon | AMRFinderPlus `--organism` | MLST scheme | Genome len (Mbp) | Provisional panel |
|---|---|---|---|---|---|
| *Pseudomonas aeruginosa* | 287 | Pseudomonas_aeruginosa | paeruginosa | 5.5–7.5 | meropenem, ceftazidime, ciprofloxacin, tobramycin, amikacin |
| *Salmonella enterica* | 28901 | Salmonella | salmonella | 4.4–5.4 | ampicillin, ceftriaxone, ciprofloxacin, chloramphenicol, TMP-SMX |
| *Enterococcus faecium* | 1352 | Enterococcus_faecium | efaecium | 2.4–3.4 | vancomycin, ampicillin, ciprofloxacin, erythromycin, tetracycline |
| *Enterobacter cloacae* | 550 | Enterobacter_cloacae | ecloacae | 4.6–6.0 | meropenem, ceftazidime, ciprofloxacin, gentamicin, TMP-SMX |
| *Streptococcus pneumoniae* | 1313 | Streptococcus_pneumoniae | spneumoniae | 1.9–2.4 | penicillin, erythromycin, tetracycline, TMP-SMX, ciprofloxacin |
| *Campylobacter jejuni* | 197 | Campylobacter | campylobacter | 1.5–2.0 | ciprofloxacin, tetracycline, erythromycin, gentamicin |

## Per-organism runbook (when BV-BRC is up)
For `<org>` in the order above:
```bash
# 1) acquire ~3,000+ balanced genomes + labels (FLAG the download size first)
python -m src.data.acquire_organism --organism <org>
# 2) annotate (AMRFinderPlus) + MLST + feature matrix + per-drug models
python -m src.organism_pipeline --organism <org>
# 3) validate: check per-drug R/S counts; drop imbalanced drugs; ⚕ review determinants
# 4) train deployable models + add to the demo registry + species markers
python -m src.models.save_models --organism <org>
```
Then add the organism to `src/app/registry.py` (display, schemes, intrinsic markers mined from its
feature matrix) so it appears in the web demo and species check. **Gate: approve each organism's
results before starting the next.**

## Honest flags / risks (per drug/organism)
- **S. pneumoniae penicillin** is driven by *pbp1a/2b/2x* mosaic genes — AMRFinderPlus may call these
  weakly; expect a modest model and report it honestly (candidate to drop).
- **P. aeruginosa ceftazidime / carbapenems** involve efflux (mexAB) + AmpC (blaPDC) derepression +
  porin (oprD) loss — partly regulatory, so determinant features will under-capture them (like the
  A. baumannii carbapenem case). ML-vs-rules gap may be the interesting story.
- **Enterococcus / Streptococcus** are Gram-positive — different determinant vocab; markers &
  MLST already scheme-supported.
- Panels use laboratory phenotypes only; any drug with <~150 R or S after acquisition is dropped.

## Status — COMPLETE (2026-07-15)
Final outcome: the demo went from **4 → 8 organisms**. Added, validated, and in the demo:
| Organism | Genomes (QC) | Drugs | Note |
|---|---|---|---|
| **P. aeruginosa** | 1,057 | 4 | ceftazidime ML-adds (regulatory); fixed intrinsic-gene rule artifact |
| **Salmonella enterica** | 2,321 | 5 | ceftriaxone 0.996; chloramphenicol honest weak spot |
| **Enterococcus faecium** | 1,890 | 3 | vancomycin 0.998 (VRE); ampicillin pbp5 (org-specific rule) |
| **Streptococcus pneumoniae** | 3,857 | 4 | penicillin pbp-mosaic; TMP-SMX co-selection caveat |

**Excluded (data-limited, honest):**
- **E. cloacae** — ≤109 genomes/drug across species/complex/genus taxa (#46).
- **C. jejuni** — ~489 genomes, only tetracycline barely qualifies (#49).

**Engineering improvements made during the expansion:** organism-specific rule baselines
(`ORG_RULES`, #45), per-organism download truncation floor (#47), markerless-organism guard fix (#43).
All 8 organisms in the demo + species check; 13 tests pass.
