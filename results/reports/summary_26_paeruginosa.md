# Summary #26 — Pseudomonas aeruginosa added (organism #5)

**Date:** 2026-07-11 · First of the Round-4 expansion. Same organism-general pipeline; only the
registry entry changed. WHO **critical** priority pathogen, complex/regulatory resistance.

## Data
- Acquired **1,160 genomes** (2.27 GB, 0 download failures) from BV-BRC; **1,057 QC-passed**
  (length 5.5–7.5 Mbp, contamination/contigs filters). AMRFinderPlus `--organism Pseudomonas_aeruginosa`
  + MLST `paeruginosa` scheme. Lab-confirmed phenotypes only.
- **Panel: 4 drugs** — meropenem, ceftazidime, ciprofloxacin, tobramycin. **Amikacin dropped**
  (126 R / 21% — too imbalanced), per the "drop rather than overclaim" rule.

## Results (unseen-lineage test, VME ≤ 3% operating point)
| Drug | ROC-AUC | PR-AUC | VME | ME | Rules ROC | What ML adds |
|---|---|---|---|---|---|---|
| meropenem | 0.841 | 0.878 | 0.025 | 0.90 | 0.740 | +0.10 |
| ceftazidime | 0.873 | 0.844 | 0.000 | 0.83 | 0.702 | **+0.17** |
| ciprofloxacin | 0.925 | 0.954 | 0.011 | 0.85 | 0.872 | +0.05 |
| tobramycin | 0.982 | 0.960 | 0.000 | 0.40 | 0.971 | +0.01 |

**Honest reading.** Strong discrimination (ROC 0.84–0.98 on unseen lineages). ML adds the most on
**ceftazidime** (+0.17) and **meropenem** (+0.10) — the drugs where P. aeruginosa resistance is
partly *regulatory* (AmpC/blaPDC derepression, MexAB efflux, OprD porin loss), which a determinant-
presence lookup cannot fully see. On **ciprofloxacin** (gyrA/parC) and **tobramycin** (acquired
aminoglycoside enzymes) the rules baseline is already strong and ML only matches it — reported
honestly, same pattern as gentamicin in K. pneumoniae. VME met (≤3%) on all four; high ME reflects the
safety-first over-calling on resistant-heavy drugs (same trade-off as elsewhere).

## Honesty checkpoint — corrected a baseline artifact (⚕ follow-up)
The first run showed rules ROC = **0.500** for ceftazidime and tobramycin. This was **an artifact of
provisional rule-gene lists**, not real: `blaPDC` (intrinsic P. aeruginosa AmpC, ~100% of genomes) and
`aph(3')-IIb` (intrinsic aminoglycoside phosphotransferase, ~99%, doesn't confer tobramycin R) were in
the baselines, so the rule called *everything* resistant → no discrimination. Removing these intrinsic/
universal genes gave the honest baselines above (0.702 and 0.971). **The corrected numbers are what is
reported; the inflated 0.98-vs-0.50 "gap" was not.** P. aeruginosa rule genes flagged for microbiology
sign-off.

## Species check
Intrinsic markers mined from the data: `aph(3')-IIb`, `catB7`, `nalC_G71E` (~95–99% in P. aeruginosa,
0% in the other panel organisms). **Fixed a `fosA` collision** — `fosA` is ~99% in *both* K. pneumoniae
and P. aeruginosa, so it was dropped from the K. pneumoniae marker set (oqxA/oqxB retained). Verified: a
P. aeruginosa genome selected as K. pneumoniae is now correctly flagged as a mismatch (via MLST and via
markers on protein input).

## Status
- ✅ In the demo (5 organisms now: K. pneumoniae, E. coli, A. baumannii, S. aureus, **P. aeruginosa**)
  + species check; 4 deployable models saved; 13 tests pass.
- Next in the one-at-a-time queue: **Salmonella enterica**.
