"""Phase 11 — interactive web demo (Streamlit).

Genome in -> per-drug resistant/susceptible + calibrated confidence + the determinants behind each
call + (optional) a Claude-generated clinical narrative. Pick a cached genome for an instant demo,
or upload a FASTA to run AMRFinderPlus live.

Run:
    streamlit run src/app/streamlit_app.py
    # optional AI narrative:  export ANTHROPIC_API_KEY=sk-ant-...  first
"""
from __future__ import annotations

import csv
import sys
import tempfile
from pathlib import Path

import joblib
import numpy as np
import streamlit as st

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from src.app.predict import DRUG_ORDER, MODELS, explain, load_determinants  # noqa: E402
from src.app import report as report_mod  # noqa: E402

AMR_CACHE = REPO / "data/interim/amrfinder"
PANEL = REPO / "data/processed/panel_labels.csv"

st.set_page_config(page_title="Reading Resistance", page_icon="🧬", layout="wide")


@st.cache_resource
def load_bundles():
    return {d: joblib.load(MODELS / f"{d}.joblib") for d in DRUG_ORDER
            if (MODELS / f"{d}.joblib").exists()}


@st.cache_data
def cached_genome_ids() -> list[str]:
    have = {p.stem for p in AMR_CACHE.glob("*.tsv")}
    if PANEL.exists():
        labelled = [r["genome_id"] for r in csv.DictReader(PANEL.open())]
        return [g for g in labelled if g in have]
    return sorted(have)


@st.cache_data
def lab_labels() -> dict[str, dict]:
    if not PANEL.exists():
        return {}
    return {r["genome_id"]: r for r in csv.DictReader(PANEL.open())}


def compute(name: str, dets: set[str]) -> dict:
    """Structured per-drug findings — same logic as the CLI report."""
    return report_mod.build_findings(name, dets)


# ---------------- UI ----------------
st.title("🧬 Reading Resistance")
st.caption("Interpretable antibiotic-resistance prediction from a *Klebsiella pneumoniae* genome — "
           "per-drug call, calibrated confidence, and the determinants behind every call.")

with st.sidebar:
    st.header("Input")
    mode = st.radio("Genome source", ["Cached genome (instant)", "Upload FASTA (live, ~1–3 min)"])
    genome_id, upload = None, None
    if mode.startswith("Cached"):
        ids = cached_genome_ids()
        genome_id = st.selectbox(f"Genome ({len(ids)} available)", ids)
    else:
        upload = st.file_uploader("Assembled genome FASTA", type=["fna", "fasta", "fa", "gz"])
    run = st.button("Predict", type="primary")
    st.markdown("---")
    st.caption("Research / decision-support only — **not a diagnostic**. Complements, does not "
               "replace, phenotypic susceptibility testing.")

if run:
    try:
        if genome_id:
            name, dets = load_determinants(genome_id, None)
        elif upload:
            suffix = ".fna.gz" if upload.name.endswith(".gz") else ".fna"
            tmp = Path(tempfile.gettempdir()) / f"upload_{upload.name}{'' if upload.name.endswith(('.gz','.fna','.fasta','.fa')) else suffix}"
            tmp.write_bytes(upload.getvalue())
            with st.spinner("Running AMRFinderPlus (~1–3 min)…"):
                name, dets = load_determinants(None, str(tmp))
        else:
            st.warning("Choose a cached genome or upload a FASTA.")
            st.stop()
    except SystemExit as e:
        st.error(str(e))
        st.stop()

    findings = compute(name, dets)
    truth = lab_labels().get(genome_id or "", {})

    st.subheader(f"Isolate {name}")
    c1, c2, c3 = st.columns(3)
    c1.metric("Organism", "K. pneumoniae")
    c2.metric("Resistance determinants", findings["n_determinants"])
    n_res = sum(1 for f in findings["findings"] if f["call"] == "Resistant")
    c3.metric("Predicted resistant", f"{n_res}/{len(findings['findings'])} drugs")

    st.markdown("### Per-drug predictions")
    for f in findings["findings"]:
        drug = f["drug"]
        col = "🔴" if f["call"] == "Resistant" else "🟢"
        flag = "  ⚠️ uncertain — recommend phenotypic testing" if f["certainty"] == "uncertain" else ""
        lab_col = drug.replace("/", "_")
        lab_val = truth.get(lab_col, "")
        lab_str = f"  ·  lab result: **{lab_val}**" if lab_val else ""
        match = ""
        if lab_val:
            match = " ✓" if lab_val == f["call"] else " ✗ (disagrees with lab)"
        with st.container(border=True):
            top = st.columns([3, 2])
            top[0].markdown(f"**{col} {drug} — {f['call'].upper()}**{flag}")
            top[1].markdown(f"P(resistant) = **{f['p_resistant']:.2f}**{lab_str}{match}")
            top[0].progress(min(max(f["p_resistant"], 0.0), 1.0))
            if f["drivers"]:
                chips = "  ".join(f"`{d['determinant']} {d['direction']}`" for d in f["drivers"])
                top[0].markdown("determinants: " + chips)
            else:
                top[0].caption("no model-relevant determinants detected for this drug")

    st.markdown("### AI clinical narrative")
    if st.button("Generate narrative (Claude)"):
        with st.spinner("Asking Claude…"):
            narrative = report_mod.llm_narrative(findings)
        if narrative:
            st.info("Generated by Claude (claude-opus-4-8) — explanation layer only; it does not "
                    "change the model's calls.")
            st.markdown(narrative)
        else:
            st.warning("ANTHROPIC_API_KEY not set — showing the deterministic templated report.")
            st.code(report_mod.templated_narrative(findings))
else:
    st.info("Pick a genome in the sidebar and click **Predict**.")
