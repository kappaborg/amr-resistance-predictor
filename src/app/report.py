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
import time
from pathlib import Path

import joblib
import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from src.app.predict import MODELS, explain, load_determinants  # noqa: E402
from src.app.registry import ORGANISMS  # noqa: E402

MODEL = "claude-opus-4-8"  # Anthropic default; explanation quality matters here


def build_findings(name: str, dets: set[str], org_key: str = "kpneu") -> dict:
    org = ORGANISMS[org_key]
    mdir = MODELS / org_key
    drug_order = list(org["drugs"])
    bundles = {d: joblib.load(mdir / f"{d}.joblib") for d in drug_order
               if (mdir / f"{d}.joblib").exists()}
    cols = next(iter(bundles.values()))["feature_cols"]
    n_in_vocab = len(set(dets) & set(cols))   # detected determinants the models actually use
    x = np.array([1 if c in dets else 0 for c in cols], dtype=int)
    findings = []
    for drug in drug_order:
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
    return {"genome": name, "organism": org["display"], "organism_key": org_key,
            "n_determinants": len(dets), "n_in_vocab": n_in_vocab, "findings": findings}


def system_prompt(organism: str) -> str:
    return (
        "You are a clinical-microbiology reporting assistant. You are given the OUTPUT of a validated, "
        f"interpretable genomic antibiotic-resistance model for a {organism} isolate: for each "
        "antibiotic, a resistant/susceptible call, a calibrated probability, a certainty flag, and the "
        "SHAP-ranked resistance determinants that drove the call. Write a concise clinician-facing report.\n"
        "STRICT RULES:\n"
        "- Do NOT change, second-guess, or recompute any call or probability — report them exactly as given.\n"
        "- For each drug, explain in one or two sentences which determinant(s) drive the call and the "
        "mechanism (e.g. blaKPC/blaNDM/blaOXA-23/blaOXA-48 = carbapenemases; gyrA/parC/grlA = "
        "fluoroquinolone target mutations; ompK/porin loss; mecA/mecC = PBP2a methicillin resistance; "
        "erm/msr/mph = macrolide-lincosamide resistance; aac/ant/aph/armA/rmt = aminoglycoside enzymes/"
        "16S methyltransferases; sul/dfr = folate pathway).\n"
        "- Explicitly flag any call marked 'uncertain' and recommend confirmatory phenotypic testing for it.\n"
        "- End with the mandatory statement that this is research/decision-support only, not a diagnostic, "
        "and does not replace phenotypic susceptibility testing.\n"
        "- Be factual and terse. No preamble, no hedging beyond the uncertainty flags."
    )


# --- Provider selection: Claude (preferred), DeepSeek (OpenAI-compatible), or Google Gemini (native) ---
# The narrative is an EXPLANATION layer only, so the provider is scientifically interchangeable. Pick
# by whichever API key is present; force one with AMR_LLM_PROVIDER=anthropic|deepseek|gemini.
# provider id -> (label, kind, default model, env keys). Order = preference (Claude > DeepSeek > Gemini).
_PROVIDERS = {
    "anthropic": ("Claude", "anthropic", MODEL, ("ANTHROPIC_API_KEY",)),
    "deepseek": ("DeepSeek", "openai", "deepseek-chat", ("DEEPSEEK_API_KEY",)),
    # gemini-flash-latest is a moving alias -> always a current model (avoids "no longer available" 404s)
    "gemini": ("Google Gemini", "gemini", "gemini-flash-latest", ("GEMINI_API_KEY", "GOOGLE_API_KEY")),
}
_MODEL_ENV = {"deepseek": "DEEPSEEK_MODEL", "gemini": "GEMINI_MODEL"}
_DEEPSEEK_URL = "https://api.deepseek.com/v1"


def _key_for(keys):
    return next((os.environ[k] for k in keys if os.environ.get(k)), None)


def _provider():
    """First configured provider in preference order (Claude > DeepSeek > Gemini), or None.
    Returns (id, label, model)."""
    forced = os.environ.get("AMR_LLM_PROVIDER", "").strip().lower()
    opts = []
    for pid, (label, _kind, default_model, keys) in _PROVIDERS.items():
        if any(os.environ.get(k) for k in keys):
            model = os.environ.get(_MODEL_ENV.get(pid, ""), default_model)
            opts.append((pid, label, model))
    if forced:
        opts = [o for o in opts if o[0] == forced]
    return opts[0] if opts else None


def _call_anthropic(key, system, user, model):
    import anthropic
    msg = anthropic.Anthropic().messages.create(
        model=model, max_tokens=1500, system=system,
        messages=[{"role": "user", "content": user}])
    return "".join(b.text for b in msg.content if b.type == "text")


def _call_openai_compat(base_url, key, system, user, model):
    import requests
    r = requests.post(f"{base_url}/chat/completions",
                      headers={"Authorization": f"Bearer {key}"},
                      json={"model": model, "max_tokens": 1500,
                            "messages": [{"role": "system", "content": system},
                                         {"role": "user", "content": user}]},
                      timeout=120)
    r.raise_for_status()
    return r.json()["choices"][0]["message"]["content"]


def _call_gemini(key, system, user, model):
    """Gemini NATIVE generateContent endpoint (works with all AI-Studio key types, unlike the
    OpenAI-compat layer). maxOutputTokens is generous because flash models spend some on reasoning."""
    import requests
    r = requests.post(
        f"https://generativelanguage.googleapis.com/v1beta/models/{model}:generateContent",
        headers={"x-goog-api-key": key, "Content-Type": "application/json"},
        json={"system_instruction": {"parts": [{"text": system}]},
              "contents": [{"role": "user", "parts": [{"text": user}]}],
              "generationConfig": {"maxOutputTokens": 4096}}, timeout=120)
    r.raise_for_status()
    cand = (r.json().get("candidates") or [{}])[0]
    text = "".join(p.get("text", "") for p in cand.get("content", {}).get("parts", []))
    if not text:
        raise RuntimeError(f"empty response (finishReason={cand.get('finishReason')})")
    return text


_CALLERS = {"anthropic": _call_anthropic, "openai": _call_openai_compat, "gemini": _call_gemini}


def llm_narrative(findings: dict) -> tuple[str | None, str]:
    """Return (narrative_text_or_None, source_label). None text -> caller shows templated fallback."""
    prov = _provider()
    if not prov:
        return None, ("no LLM API key set — using the deterministic templated report "
                      "(set ANTHROPIC_API_KEY, DEEPSEEK_API_KEY, or GEMINI_API_KEY to enable an AI narrative)")
    pid, label, model = prov
    _lbl, kind, _dm, keys = _PROVIDERS[pid]
    key = _key_for(keys)
    system = system_prompt(findings.get("organism", "bacterial"))
    user = ("Model output (JSON):\n" + json.dumps(findings, indent=2) +
            "\n\nWrite the clinician report.")
    # Retry TRANSIENT server errors (503/502/500/429/timeout) with backoff — these are common on the
    # free Gemini tier under load — before giving up to the deterministic templated report.
    transient = ("503", "502", "500", "429", "overloaded", "unavailable", "timeout", "timed out")
    last = None
    for attempt in range(3):
        try:
            if kind == "openai":
                text = _call_openai_compat(_DEEPSEEK_URL, key, system, user, model)
            else:
                text = _CALLERS[kind](key, system, user, model)
            return text, f"{label} ({model})"
        except Exception as e:  # noqa: BLE001
            last = e
            if attempt < 2 and any(t in str(e).lower() for t in transient):
                time.sleep(2 * (attempt + 1))   # 2s, 4s
                continue
            break
    return None, f"{label} call failed: {str(last)[:140]} — using templated fallback"


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
    ap.add_argument("--organism", choices=list(ORGANISMS), default="kpneu")
    ap.add_argument("--genome")
    ap.add_argument("--genome-id")
    ap.add_argument("--json", action="store_true", help="also print the structured findings JSON")
    args = ap.parse_args()
    if not (args.genome or args.genome_id):
        sys.exit("provide --genome or --genome-id")

    name, dets = load_determinants(args.genome_id, args.genome, ORGANISMS[args.organism]["amrfinder"])
    findings = build_findings(name, dets, args.organism)
    if args.json:
        print(json.dumps(findings, indent=2) + "\n")

    narrative, source = llm_narrative(findings)
    print("=" * 70)
    print(f"  AI EXPLANATION LAYER — source: {source}")
    print("=" * 70 + "\n")
    print(narrative or templated_narrative(findings))
    return 0


if __name__ == "__main__":
    sys.exit(main())
