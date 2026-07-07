# Summary #1 — Phase 1: Organism & Drug Selection (BV-BRC survey)

**Date:** 2026-07-07 · **Status:** ⏸ Awaiting ⚕ microbiology sign-off (GATE) · Counts complete.

## What was done
- Verified the BV-BRC Data API contract live, then surveyed `genome_amr` for 9 candidate organisms.
- **Counted lab phenotypes only** (`evidence="Laboratory Method"`) — computational AMR predictions
  excluded, because using them as labels would be circular. (E. coli: 6.3M rows → 243K lab.)
- Collected **221 (organism, drug) count rows** → `results/metrics/phase1_counts.csv`.
- Tuned thresholds from the actual distribution and generated `results/reports/phase1_shortlist.md`.
- **No genome downloads** — metadata counts only.

## Threshold tuning (data-driven)
Distribution across 221 drug-counts: min(R,S) median 191 / p75 501 / max 6060; balance median 0.15.
→ Data is abundant, so **balance is the binding quality lever, not raw count.**
**Chosen: min(R,S) ≥ 150 · balance ≥ 0.30 · total(R+S) ≥ 300** → 51 balanced drugs across 7 organisms.

## Shortlist — qualifying drugs per organism
| Organism | # qualifying | Notable balanced drugs |
|---|---|---|
| Escherichia coli | 11 | ampicillin, amox/clav, TMP-SMX, ceftriaxone |
| Klebsiella pneumoniae | 11 | meropenem, imipenem, gentamicin, ciprofloxacin, cefoxitin, TMP-SMX |
| Mycobacterium tuberculosis | 6 | isoniazid, streptomycin, moxifloxacin, ethionamide, rifabutin |
| Staphylococcus aureus | 6 | oxacillin, cefoxitin, methicillin, erythromycin, ciprofloxacin |
| Acinetobacter baumannii | 5 | imipenem, meropenem, amikacin, tobramycin |
| Pseudomonas aeruginosa | 5 | meropenem, ceftazidime, ciprofloxacin, levofloxacin |
| Salmonella enterica | 3 | ampicillin, tetracycline, sulfisoxazole |
| Neisseria gonorrhoeae | 2 | ciprofloxacin, penicillin |
| Streptococcus pneumoniae | 2 | TMP-SMX, tetracycline |

## Recommendation (⚕ = co-builder's call)
**Primary organism: _Klebsiella pneumoniae_.** Clinically urgent carbapenem resistance
(meropenem/imipenem), well-defined MLST lineages (ST258/ST11/etc. → strong phylogeny-aware-split
story), and clear determinants (blaKPC/NDM/OXA, aminoglycoside + fluoroquinolone mechanisms) for a
rich SHAP/interpretability narrative. E. coli is an equally strong runner-up (most data).
- **Proposed 5-drug panel** (diverse mechanisms): meropenem · gentamicin · ciprofloxacin ·
  trimethoprim-sulfamethoxazole · cefoxitin.
- **Week-1 thin-slice drug:** ciprofloxacin (exercises point-mutation features gyrA/parC) *or*
  meropenem (clinically weightiest; VME matters most).

## Blocked / needs decision (GATE)
1. ⚕ **Organism** — confirm K. pneumoniae (or pick another from the shortlist).
2. ⚕ **Drug panel** (3–6) and the single **thin-slice drug**.
3. ⚕ **Intermediate** handling — currently excluded from counts; final rule set in Phase 3.
Once decided → write `organism`/`drugs` into config, log it, then Phase 2 (download, size flagged first).

## Caveats
- Counts are lab-phenotype **rows**; a genome may be tested more than once → dedup in Phase 2/3.
- Balance/abundance are shortlist proxies; true **lineage diversity** (can we split without leakage?)
  is only confirmed after MLST in Phase 5.
