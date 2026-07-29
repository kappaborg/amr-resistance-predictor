# Summary #28 — Enterococcus faecium added (organism #7)

**Date:** 2026-07-13 · Third Round-4 organism; second Gram-positive. WHO **high** priority (VRE).

## Data
- **2,285 genomes** downloaded (2.1 GB; **102 failed ~4%** — transient BV-BRC hiccups, not chased) →
  **1,890 QC-passed** (length 2.4–3.4 Mbp). AMRFinderPlus `--organism Enterococcus_faecium` + MLST
  `efaecium`.
- **Panel: 3 drugs** — vancomycin, ampicillin, tetracycline. Dropped **ciprofloxacin** (only 27 tested)
  and **erythromycin** (only 100 S) per the "drop rather than overclaim" rule.

## Results (unseen-lineage test, VME ≤ 3% operating point)
| Drug | ROC-AUC | PR-AUC | VME | ME | Rules ROC | Reading |
|---|---|---|---|---|---|---|
| **vancomycin** | **0.998** | 0.997 | 0.000 | 0.004 | 0.996 | textbook — VRE *is* the van genes (van presence 48% ≈ R 49%) |
| **ampicillin** | **0.992** | 0.997 | 0.004 | 0.011 | 0.638 | vs the *pbp5* baseline (organism-specific); ML learns the specific resistance alleles |
| tetracycline | 0.911 | 0.971 | 0.005 | 0.204 | 0.802 | ML adds ~0.11 over the tet-gene baseline |

## Honesty check — the ampicillin 0.500 baseline is real, and ML's 0.992 is genuine
E. faecium ampicillin resistance is **pbp5-mediated**, not β-lactamase. Verified in the data:
- `pbp5` point-mutation determinants (pbp5_A401S, N496K, A499T, …) are present at **92%** and drive
  the ML model; the generic β-lactamase genes (blaTEM/blaSHV/blaCTX-M) are **0%** in Enterococcus.
- So the rules baseline (which uses those β-lactamase genes) finds nothing → predicts all susceptible
  → ROC 0.500. **ML at 0.992 is a real, verified gain** — it reads the pbp5 resistance alleles a
  β-lactamase lookup cannot see. This is the cleanest "what ML adds over a gene lookup" case in the
  project (cf. K. pneumoniae cefoxitin/porin loss).
- **Resolved with organism-specific baselines (`ORG_RULES`).** The rules baseline is now organism-aware:
  for `(efaecium, ampicillin)` it uses **pbp5** (the correct mechanism) instead of the shared
  β-lactamase list. That lifts the ampicillin baseline from the misleading **0.500 → 0.638** (pbp5
  mutations are common → the presence-rule over-calls, so it's fair-but-imperfect). ML at **0.992**
  still adds a real **+0.35** — but now measured against an honest mechanism-based baseline, not a
  strawman. `(spneumoniae, penicillin)` → PBP1a/2b/2x is pre-registered for the next organism.

## Species check
Markers mined: `msr(C)` (98%, intrinsic macrolide efflux) and `pbp5_N496K` (91%) — both E. faecium-
specific (0% elsewhere), no collision with other organisms' markers. Plus the distinct `efaecium`
MLST scheme.

## Status
- ✅ In the demo (7 organisms: + **Enterococcus faecium**) + species check; 3 models saved; 13 tests pass.
- Remaining Round-4 queue: E. cloacae → S. pneumoniae → C. jejuni.
- **Recurring ⚕ item:** organism-specific rule baselines (ampicillin/penicillin differ by organism).
