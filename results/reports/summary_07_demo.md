# Summary #7 — Week 3: Demo Interface (Phase 11)  ✅

**Date:** 2026-07-08 · **Status:** working demo — genome in → per-drug R/S + confidence + determinants.

## What it does
`src/app/predict.py` takes a K. pneumoniae genome and returns, for each of the 5 drugs:
- **Resistant / susceptible call** at the clinical VME≤3% operating threshold,
- **Calibrated P(resistant)** (isotonic-calibrated), with an `[uncertain]` flag when 0.2 < P < 0.8,
- **The determinants that drove the call**, SHAP-ranked (the "why this strain resists drug X").

Two input modes:
- `--genome path.fna[.gz]` — runs AMRFinderPlus live (~1–3 min).
- `--genome-id 573.xxxx` — reuses a cached annotation (instant) for genomes we already have.

Deployment models (`results/models/<drug>.joblib`): XGBoost fit on all labeled data + isotonic
calibrator + VME-tuned threshold + feature vocabulary (via `src/models/save_models.py`).

## Verified behaviour
- **MDR genome 573.12772** → all 5 drugs RESISTANT, drivers cite the correct mechanisms:
  meropenem→blaKPC-3, gentamicin→aac(3)-IVa, ciprofloxacin→parC_S80I/gyrA_S83I, TMP-SMX→dfrA12/sul1.
- **Pan-susceptible genome 573.12859** → 4/5 susceptible with low P(resistant) (0.00–0.19).
- Live-annotation path reproduces the cached-path result exactly.

## Honest limitation surfaced by the demo
The susceptible genome gets a **cefoxitin "RESISTANT [uncertain]" at P=0.21** — a direct consequence
of cefoxitin's data-limited, low VME-threshold (0.08) from Week 2. The demo faithfully shows this
rather than hiding it; the Phase-2b top-up (downloading now) should stabilise the cefoxitin operating
point.

## Ethics guard (shown in output)
"Research/decision-support only — not a substitute for phenotypic susceptibility testing."
