# Summary #14 — Interactive Web Demo (Streamlit)

**Date:** 2026-07-08 · **Status:** ✅ built & boots (verified HTTP 200 headless).

## What it is
`src/app/streamlit_app.py` — a one-screen web app for the poster/defense. Pick a cached genome
(instant) or upload a FASTA (live AMRFinderPlus, ~1–3 min), and see per drug:
- resistant/susceptible **call** (colour-coded) at the clinical VME≤3% threshold,
- **calibrated P(resistant)** with a progress bar and an `uncertain` flag,
- the **determinants** that drove the call (SHAP-ranked chips),
- the **lab result** side-by-side when known (shows model-vs-truth agreement ✓/✗),
- a one-click **Claude clinical narrative** (explanation layer; falls back to a templated report
  offline).

```bash
streamlit run src/app/streamlit_app.py
# optional: export ANTHROPIC_API_KEY=sk-ant-...  for the AI narrative
```

## Why it matters (competition)
Communication is an explicitly judged criterion. A clickable app that shows a genome → calibrated,
uncertainty-aware, interpretable calls (with the lab result next to the prediction) lands far better
in a poster session than a CLI, and ties together every feature — calibration, the conformal
`uncertain` flag, SHAP interpretability, and the AI narrative — in one place. Reuses the exact
prediction machinery (`src/app/predict.py`, `src/app/report.py`); no duplicate logic. Streamlit
pinned in `environment.yml`.
