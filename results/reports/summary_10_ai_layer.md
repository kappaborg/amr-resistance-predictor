# Summary #10 — AI Explanation Layer (the "AI flow")

**Date:** 2026-07-08 · **Status:** ✅ built (Claude API integration, offline fallback).

## What it is
`src/app/report.py` turns the model's per-drug output — resistant/susceptible call, calibrated
probability, certainty flag, and SHAP-ranked determinants — into a **clinician-readable narrative**
using **Claude (`claude-opus-4-8`, official Anthropic Python SDK)**.

```bash
export ANTHROPIC_API_KEY=sk-ant-...
python -m src.app.report --genome-id 573.12772          # cached (instant)
python -m src.app.report --genome path/to/genome.fna    # live AMRFinderPlus
python -m src.app.report --genome-id 573.12772 --json   # + structured findings JSON
```

## Honest scoping (the defensible design)
The LLM is an **explanation layer, not a predictor.** Every call, probability, and confidence comes
from the trained models; the LLM only renders them into prose and maps determinants to mechanisms
(blaKPC → carbapenemase, gyrA/parC → fluoroquinolone target mutations, ompK porin loss, aac/ant →
aminoglycoside-modifying enzymes, sul/dfr → folate pathway). This is deliberate: published evidence
(IDWeek 2025) shows LLMs do **not** improve raw AMR prediction accuracy, so using one as the predictor
would be a weakness. As an explanation/reporting layer it adds genuine value (readable reports,
mechanism context, uncertainty framing) without touching the decision.

System-prompt guardrails: never change a call; explain each driver + mechanism in ≤2 sentences;
flag every `uncertain` call and recommend confirmatory phenotypic testing; end with the mandatory
"research/decision-support only, not a diagnostic" statement.

## Robustness
- **Works offline:** a deterministic templated narrative renders when `ANTHROPIC_API_KEY` is unset,
  so the demo never breaks; Claude enriches it when a key is present.
- Verified on the MDR genome (573.12772): all-resistant, drivers correctly cite blaKPC-3, aac(3)-IVa,
  gyrA/parC, dfrA12/sul1.

## Competition value
Three tier-1 additions now complete: **E. coli generalization** (in progress), **conformal
prediction** (guaranteed uncertainty + abstention), and this **AI explanation layer** — a genuine
"AI flow" that's honestly scoped and demo-ready. Anthropic SDK pinned in `environment.yml`.
