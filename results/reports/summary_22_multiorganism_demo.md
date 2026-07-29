# Summary #22 — Multi-Organism Web Demo (all four pathogens live)

**Date:** 2026-07-10 · The demo now serves **all four organisms** we modelled, not just K. pneumoniae.

## What changed
Previously the website (`streamlit_app.py`) and CLI predictor were hardcoded to *K. pneumoniae*
(5 drugs). The other three validated organisms existed only as evaluation numbers. Now:

- **Shared registry** (`src/app/registry.py`): one source of truth mapping each organism → display
  name, NCBI taxon, AMRFinderPlus `--organism` string (E. coli = "Escherichia", the earlier gotcha),
  feature/lineage files, label source, and per-drug BV-BRC encodings. `save_models`, `predict`, and
  the Streamlit app all import it, so they cannot drift.
- **20 deployable models** trained + saved under `results/models/<organism>/<drug>.joblib` — each a
  bundle of {calibrated XGBoost, isotonic calibrator, VME-capped threshold, feature vocabulary,
  AMRFinderPlus organism string, known-gene rules}.
- **`predict.py` / `report.py` / `streamlit_app.py`** all take `--organism` (app: a dropdown). The
  app boots clean (Streamlit 1.59.1, headless verified).

## Coverage now live

| Organism | Drugs served | Notable interpretable calls (spot-checked) |
|---|---|---|
| *Klebsiella pneumoniae* | meropenem, gentamicin, ciprofloxacin, TMP-SMX, cefoxitin | blaKPC-2 → meropenem R; dfrA/sul → TMP-SMX R |
| *Escherichia coli* | ciprofloxacin, gentamicin, TMP-SMX, ampicillin, ceftriaxone | folate-pathway genes → TMP-SMX R |
| *Acinetobacter baumannii* | meropenem, imipenem, ciprofloxacin, gentamicin, amikacin | blaOXA carbapenemases drive carbapenem calls |
| *Staphylococcus aureus* | oxacillin, cefoxitin, ciprofloxacin, erythromycin, clindamycin | **mecA → oxacillin R**; gyrA_S84L/parC_S80F → cipro R |

Every call still comes with: calibrated P(resistant), the VME-capped R/S decision, the SHAP-ranked
determinants that drove it, an uncertainty flag, and an optional Claude narrative (explanation-only).
Lab-truth overlay is shown for K. pneumoniae (saved panel); the others predict live.

## How to run
```
python -m src.models.save_models              # (re)train all 20 models
python -m src.app.predict --organism saureus --genome-id 1280.10000
python -m streamlit run src/app/streamlit_app.py
```

## Why this matters for the competition
A judge can now pick any of four WHO-priority pathogens across the Gram divide and get an
interpretable, calibrated, determinant-backed prediction from a genome — a four-pathogen live demo,
not a single-organism toy. The generalization claim (Summary #18) is now something you can *show*,
not just report.
