"""Phase 3 feature — AI explanation layer: turn the model's per-drug predictions + SHAP
determinants into a clinician-readable narrative report, using Claude (Anthropic API).

IMPORTANT — scope: the LLM is an EXPLANATION layer, not a predictor. All resistance calls,
probabilities, and confidence come from the trained models (`src/app/predict.py`); the LLM only
renders them into prose and situates the determinants in known biology. (Published evidence shows
LLMs don't improve raw AMR prediction accuracy, so we deliberately keep them out of the decision.)

Falls back to a deterministic templated narrative when ANTHROPIC_API_KEY is not set, so the demo
works offline; the LLM enriches it when a key is present.

Usage:
    export ANTHROPIC_API_KEY=sk-ant-...
    python -m src.app.report --genome-id 573.12772
    python -m src.app.report --genome path/to/genome.fna     # live AMRFinderPlus
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import joblib
import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from src.app.predict import DRUG_ORDER, MODELS, explain, load_determinants  # noqa: E402

MODEL = "claude-opus-4-8"  # Anthropic default; explanation quality matters here


def build_findings(name: str, dets: set[str]) -> dict:
    bundles = {d: joblib.load(MODELS / f"{d}.joblib") for d in DRUG_ORDER
               if (MODELS / f"{d}.joblib").exists()}
    cols = next(iter(bundles.values()))["feature_cols"]
    x = np.array([1 if c in dets else 0 for c in cols], dtype=int)
    findings = []
    for drug in DRUG_ORDER:
        if drug not in bundles:
            continue
        b = bundles[drug]
        raw = b["model"].predict_proba(x.reshape(1, -1))[0, 1]
        prob = float(b["calibrator"].transform([raw])[0])
        call = "Resistant" if prob >= b["threshold"] else "Susceptible"
        certainty = "high" if (prob >= 0.8 or prob <= 0.2) else "uncertain"
        drivers = [{"determinant": d, "direction": arrow}
                   for d, arrow, _ in explain(b, x, dets)]
        findings.append({"drug": drug.replace("_", "/"), "call": call,
                         "p_resistant": round(prob, 2), "certainty": certainty,
                         "drivers": drivers})
    return {"genome": name, "organism": "Klebsiella pneumoniae",
            "n_determinants": len(dets), "findings": findings}


SYSTEM = (
    "You are a clinical-microbiology reporting assistant. You are given the OUTPUT of a validated, "
    "interpretable genomic antibiotic-resistance model for a Klebsiella pneumoniae isolate: for each "
    "antibiotic, a resistant/susceptible call, a calibrated probability, a certainty flag, and the "
    "SHAP-ranked resistance determinants that drove the call. Write a concise clinician-facing report.\n"
    "STRICT RULES:\n"
    "- Do NOT change, second-guess, or recompute any call or probability — report them exactly as given.\n"
    "- For each drug, explain in one or two sentences which determinant(s) drive the call and the "
    "mechanism (e.g. blaKPC = carbapenemase; gyrA/parC = fluoroquinolone target mutations; ompK "
    "porin loss; aac/ant = aminoglycoside-modifying enzymes; sul/dfr = folate pathway).\n"
    "- Explicitly flag any call marked 'uncertain' and recommend confirmatory phenotypic testing for it.\n"
    "- End with the mandatory statement that this is research/decision-support only, not a diagnostic, "
    "and does not replace phenotypic susceptibility testing.\n"
    "- Be factual and terse. No preamble, no hedging beyond the uncertainty flags."
)


def llm_narrative(findings: dict) -> str | None:
    if not os.environ.get("ANTHROPIC_API_KEY"):
        return None
    import anthropic
    client = anthropic.Anthropic()
    msg = client.messages.create(
        model=MODEL, max_tokens=1500,
        system=SYSTEM,
        messages=[{"role": "user", "content":
                   "Model output (JSON):\n" + json.dumps(findings, indent=2) +
                   "\n\nWrite the clinician report."}],
    )
    return "".join(b.text for b in msg.content if b.type == "text")


def templated_narrative(findings: dict) -> str:
    lines = [f"AMR PREDICTION REPORT — {findings['genome']} ({findings['organism']})",
             f"{findings['n_determinants']} resistance determinants detected.\n"]
    for f in findings["findings"]:
        tag = "  [UNCERTAIN — recommend phenotypic testing]" if f["certainty"] == "uncertain" else ""
        drv = ", ".join(f"{d['determinant']} {d['direction']}" for d in f["drivers"]) or "no model-relevant determinants"
        lines.append(f"- {f['drug']}: {f['call'].upper()} (P={f['p_resistant']}){tag}\n    drivers: {drv}")
    lines.append("\nResearch/decision-support only — not a diagnostic; does not replace phenotypic "
                 "susceptibility testing.")
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser(description="AI clinical narrative for a genome's AMR prediction.")
    ap.add_argument("--genome")
    ap.add_argument("--genome-id")
    ap.add_argument("--json", action="store_true", help="also print the structured findings JSON")
    args = ap.parse_args()
    if not (args.genome or args.genome_id):
        sys.exit("provide --genome or --genome-id")

    name, dets = load_determinants(args.genome_id, args.genome)
    findings = build_findings(name, dets)
    if args.json:
        print(json.dumps(findings, indent=2) + "\n")

    narrative = llm_narrative(findings)
    source = "Claude (claude-opus-4-8)" if narrative else "templated fallback (set ANTHROPIC_API_KEY for AI narrative)"
    print("=" * 70)
    print(f"  AI EXPLANATION LAYER — source: {source}")
    print("=" * 70 + "\n")
    print(narrative or templated_narrative(findings))
    return 0


if __name__ == "__main__":
    sys.exit(main())
