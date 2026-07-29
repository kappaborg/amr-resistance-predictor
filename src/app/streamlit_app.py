"""Phase 11 — interactive web demo (Streamlit), all eight organisms.

Genome in -> per-drug resistant/susceptible + calibrated confidence + the determinants behind each
call + (optional) a Claude-generated clinical narrative. Pick the organism, then a cached genome for
an instant demo, or upload a FASTA to run AMRFinderPlus live.

Run:
    python -m streamlit run src/app/streamlit_app.py
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
from src.app.predict import MODELS, genome_scheme, load_determinants, species_verdict  # noqa: E402
from src.app.registry import ORGANISMS, ORG_ORDER  # noqa: E402
from src.app import report as report_mod  # noqa: E402

AMR_CACHE = REPO / "data/interim/amrfinder"
PANEL = REPO / "data/processed/panel_labels.csv"   # lab truth available for K. pneumoniae only

st.set_page_config(page_title="Reading Resistance", page_icon="🧬", layout="wide")


@st.cache_data
def cached_genome_ids(taxon: int) -> list[str]:
    """Cached genome_ids for this organism (TSV files named <taxon>.<id>.tsv)."""
    have = sorted(p.stem for p in AMR_CACHE.glob(f"{taxon}.*.tsv"))
    return have


@st.cache_data
def lab_labels() -> dict[str, dict]:
    if not PANEL.exists():
        return {}
    return {r["genome_id"]: r for r in csv.DictReader(PANEL.open())}


# ---------------- UI ----------------
st.title("🧬 Reading Resistance")
st.caption("Interpretable antibiotic-resistance prediction from a bacterial genome — per-drug call, "
           "calibrated confidence, and the determinants behind every call. Eight WHO-priority pathogens.")

with st.sidebar:
    st.header("Input")
    org_key = st.selectbox("Organism", ORG_ORDER,
                           format_func=lambda k: ORGANISMS[k]["display"])
    org = ORGANISMS[org_key]
    mode = st.radio("Genome source", ["Cached genome (instant)", "Upload FASTA (live, ~1–3 min)"])
    genome_id, upload = None, None
    if mode.startswith("Cached"):
        ids = cached_genome_ids(org["taxon"])
        genome_id = st.selectbox(f"Genome ({len(ids)} available)", ids)
    else:
        upload = st.file_uploader("Genome FASTA", type=["fna", "fasta", "fa", "gz", "faa"],
                                  help="Assembled nucleotide genome (.fna/.fasta/.fa, plain or "
                                       "gzipped .gz) — or a protein FASTA (.faa, predicted proteins). "
                                       "Molecule type is auto-detected. Typical size 2–6 MB.")
        st.caption(f"Will run `AMRFinderPlus --organism {org['amrfinder']}` on your genome "
                   "(~1–3 min). Nothing is uploaded to any external service.")
    run = st.button("Predict", type="primary")
    st.markdown("---")
    st.caption("Research / decision-support only — **not a diagnostic**. Complements, does not "
               "replace, phenotypic susceptibility testing.")

# --- On Predict: compute once and STORE in session_state so results survive later reruns
# (clicking "Generate AI narrative" reruns the whole script; without this the results vanish). ---
if run:
    try:
        if genome_id:
            name, dets = load_determinants(genome_id, None, org["amrfinder"])
            scheme = genome_scheme(genome_id, None)
        elif upload:
            # keep the real extension so gzip/plain is detected; sanitise the rest of the name
            safe = "".join(c for c in upload.name if c.isalnum() or c in "._-") or "genome.fna"
            tmp = Path(tempfile.gettempdir()) / f"upload_{safe}"
            tmp.write_bytes(upload.getvalue())
            with st.spinner(f"Running AMRFinderPlus (--organism {org['amrfinder']}, ~1–3 min)…"):
                name, dets = load_determinants(None, str(tmp), org["amrfinder"])
            with st.spinner("Verifying species (MLST)…"):
                scheme = genome_scheme(None, str(tmp))
            tmp.unlink(missing_ok=True)
        else:
            st.warning("Choose a cached genome or upload a FASTA.")
            st.stop()
    except SystemExit as e:
        st.error(str(e))
        st.stop()
    st.session_state["result"] = {
        "name": name, "org_key": org_key, "genome_id": genome_id, "scheme": scheme,
        "dets": sorted(dets),
        "findings": report_mod.build_findings(name, dets, org_key),
    }
    st.session_state.pop("narrative", None)   # drop any stale narrative from a previous genome

res = st.session_state.get("result")
if not res:
    st.info("Pick an organism and genome in the sidebar, then click **Predict**.")
else:
    findings = res["findings"]
    r_org = ORGANISMS[res["org_key"]]
    truth = lab_labels().get(res["genome_id"] or "", {}) if res["org_key"] == "kpneu" else {}

    status, msg = species_verdict(res["org_key"], res.get("scheme"), set(res.get("dets", [])))

    st.subheader(f"Isolate {res['name']}")
    c1, c2, c3 = st.columns(3)
    c1.metric("Organism", r_org["display"])
    c2.metric("Resistance determinants", findings["n_determinants"])
    n_res = sum(1 for f in findings["findings"] if f["call"] == "Resistant")
    c3.metric("Predicted resistant", f"{n_res}/{len(findings['findings'])} drugs")

    # Guard 1 — species mismatch: the genome types (MLST) to a DIFFERENT panel organism.
    if status == "mismatch":
        st.error(f"**Wrong species — {msg}.**\n\n"
                 f"The *{r_org['display']}* models cannot make a valid prediction for a different "
                 "species. Predictions are withheld. Select the matching organism in the sidebar, or "
                 f"upload a *{r_org['display']}* genome.")
        st.stop()

    # Guard 2 — no in-vocabulary determinants. Only a wrong-species signal for organisms that carry
    # near-universal INTRINSIC determinants (kpneu fosA/oqxA, P. aeruginosa blaPDC, etc.). Organisms
    # without such a marker (E. coli, Salmonella) can be genuinely pan-susceptible with 0 determinants,
    # so there we predict "susceptible to all" instead of rejecting.
    if findings.get("n_in_vocab", 0) == 0 and not r_org.get("markers"):
        st.info(
            "**No known resistance determinants detected** — the calibrated P(resistant) is low for "
            f"every drug, so this is most likely a pan-susceptible *{r_org['display']}* isolate. "
            "Read the probabilities below: the R/S labels use a safety-first VME≤3% threshold that "
            "deliberately **over-calls** low-prevalence drugs (so a drug may read RESISTANT at a low "
            "probability)."
            + ("  Species could not be MLST-confirmed for this input." if status == "unverified" else ""))
    elif findings.get("n_in_vocab", 0) == 0:
        st.error(
            f"**No known resistance determinants were found — this does not look like a "
            f"*{r_org['display']}* genome.**\n\n"
            "A genuine strain of this species carries intrinsic determinants, so finding none usually "
            "means the upload is a **different species or organism** (or not an assembled genome). "
            "The models assume the input *is* the selected organism and only see catalogued resistance "
            "genes, so **no meaningful prediction can be made here** — the per-drug calls are withheld. "
            "Upload an assembled genome of the selected species, or switch the organism in the sidebar."
        )
        st.stop()

    if status == "ok":
        st.success(f"✓ Species check: {msg}.")
    else:  # unverified (e.g. protein input, or an ambiguous/unmapped scheme)
        st.caption(f"ℹ️ {msg}. Predictions assume the input is *{r_org['display']}*.")

    st.markdown("### Per-drug predictions")
    for f in findings["findings"]:
        drug = f["drug"]
        col = "🔴" if f["call"] == "Resistant" else "🟢"
        flag = "  ⚠️ uncertain — recommend phenotypic testing" if f["certainty"] == "uncertain" else ""
        lab_val = truth.get(drug.replace("/", "_"), "")
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
    st.caption("Optional. Uses Claude, DeepSeek, or Google Gemini — whichever API key is set "
               "(`ANTHROPIC_API_KEY` / `DEEPSEEK_API_KEY` / `GEMINI_API_KEY`). Explanation layer "
               "only: it never changes the model's calls. With no key, a deterministic report is shown.")
    if st.button("Generate AI narrative"):
        with st.spinner("Asking the language model…"):
            st.session_state["narrative"] = report_mod.llm_narrative(findings)

    if "narrative" in st.session_state:
        narrative, source = st.session_state["narrative"]
        if narrative:
            st.info(f"Generated by **{source}** — explanation layer only; it does not change the "
                    "model's calls.")
            st.markdown(narrative)
        else:
            st.warning(source)
            st.code(report_mod.templated_narrative(findings))
