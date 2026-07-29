# Summary #27 — Salmonella enterica added (organism #6)

**Date:** 2026-07-13 · Second Round-4 organism. Same pipeline; registry entry only.

## Data
- **3,109 genomes** (4.3 GB, 0 failures; download paused/resumed cleanly through a BV-BRC slow spell) →
  **2,321 QC-passed** (length 4.4–5.4 Mbp). AMRFinderPlus `--organism Salmonella` + MLST `salmonella`.
- **Panel: 5 drugs** — ampicillin, ceftriaxone, ciprofloxacin, chloramphenicol, TMP-SMX. All cleared
  ≥150 R and ≥150 S, though several are **low-R-prevalence** (cipro 9%, TMP-SMX 11%, ceftriaxone 17%).

## Results (unseen-lineage test, VME ≤ 3% operating point)
| Drug | ROC-AUC | PR-AUC | VME | ME | Rules ROC | Reading |
|---|---|---|---|---|---|---|
| **ceftriaxone** | **0.996** | 0.953 | 0.000 | 0.044 | 0.983 | excellent; ESBL/AmpC directly detectable |
| ciprofloxacin | 0.978 | 0.963 | 0.030 | 0.482 | 0.975 | strong; ML ≈ rules (gyrA/parC + qnr) |
| ampicillin | 0.947 | 0.931 | 0.025 | 0.832 | 0.952 | strong; ML ≈ rules (blaTEM etc.) |
| **TMP-SMX** | 0.832 | 0.441 | 0.037 | 0.722 | **0.565** | **ML adds real value** (+0.27 over a weak sul/dfr baseline) |
| chloramphenicol | 0.880 | 0.632 | **0.109** | 0.537 | 0.918 | **honest weak spot** — ML < rules, VME misses target |

**Honesty checks (done, not skipped):**
- **No intrinsic-gene artifact.** Unlike P. aeruginosa (blaPDC), every Salmonella rule gene is at
  modest prevalence (blaTEM 29%, sul2 33%, floR 13%) — so the rules baselines are genuine. The low
  TMP-SMX rules AUC (0.565) is real (sul/dfr only weakly predictive here), so ML's 0.832 is a genuine
  gain; the high chloramphenicol rules AUC (0.918) is real too, so ML underperforming it is an honest
  weakness, not an artifact.
- **chloramphenicol flagged**: ML (0.880) is *worse* than the gene lookup (0.918) **and** VME 10.9%
  exceeds the 3% target — for this drug a determinant lookup is the better tool; reported, not hidden.
- **No reliable Salmonella species marker** (Enterobacterales, like E. coli) — species relies on its
  distinct MLST `salmonella` scheme; no collision with other organisms' markers.

## Demo-behaviour fix (markerless organisms + imbalanced drugs)
Adding Salmonella surfaced two demo issues, now fixed:
1. **Zero-determinant guard was too aggressive.** It rejected a *susceptible* Salmonella genome (0
   acquired determinants) as "not Salmonella". Fixed: the reject-on-zero rule now applies only to
   organisms with **near-universal intrinsic markers** (kpneu/paeruginosa/abaumannii/saureus). For
   markerless organisms (E. coli, Salmonella) 0 determinants = a genuinely possible pan-susceptible
   isolate, so it predicts instead of rejecting.
2. **VME≤3% over-calls hard on low-prevalence drugs.** Salmonella thresholds fall to ~0.02–0.05, so
   the safety-first operating point labels almost everything RESISTANT — even a determinant-free
   genome (P≈0.05) reads RESISTANT. The calibrated **probability** stays honest (low = susceptible);
   the demo note now tells the reader to read P(resistant), not just the label. *Open design question:
   whether the demo should show a balanced threshold instead of the clinical VME≤3% one.*

## Status
- ✅ In the demo (6 organisms: + **Salmonella enterica**) + species check; 5 models saved; 13 tests pass.
- Next queued: E. faecium.
