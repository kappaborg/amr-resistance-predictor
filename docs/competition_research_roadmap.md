# Competition Boost — Research Roadmap (what to add & why)

Deep-search synthesis of what differentiates a winning genome-based AMR-prediction project (2024–2026),
mapped to our pipeline. Ranked by winning-impact ÷ effort. Our current base: 1 organism (K. pneumoniae),
5 drugs, determinant features, phylogeny-aware split, SHAP, calibrated VME≤3% operating point, demo.

Note: genome-AMR prediction is an active *competition* space — e.g. CAMDA 2025 ran an AMR phenotype
prediction challenge — so judges reward generalization, rigor, and novelty, not just accuracy.

---

## TIER 1 — highest impact, feasible in the sprint (recommended)

### 1. Prove generalization on a **second organism (E. coli)**
**Why it wins:** the #1 thing that turns "a model" into "a system" for judges is showing it transfers.
Our pipeline is already organism-agnostic — only `taxon_id`, `--organism`, and the MLST scheme change.
We surveyed E. coli in Phase 1: 11 balanced drugs, the most data of any candidate. Re-running the exact
pipeline on E. coli demonstrates the *method* generalizes without re-engineering.
**Effort:** ~1 day (data pull + annotate + the same scripts). **Impact:** very high.
**Refs:** cross-species AMR is an out-of-distribution problem; species-independent pipelines that hold
up across taxa are the credible frontier (Frontiers cellular-infect-micro 2024; ESKAPE ML 2025).

### 2. **Conformal prediction + selective abstention** (rigorous uncertainty + "defer to lab")
**Why it wins:** extends our existing clinical-safety differentiator (VME focus) with *statistically
guaranteed* uncertainty and an explicit **reject option** — the model abstains ("insufficient
confidence → send for phenotypic testing") instead of guessing. Directly addresses "a wrong AMR call
is life-threatening." Feasible on our existing models (inductive/Mondrian conformal, e.g. MAPIE).
**Effort:** ~1–2 days. **Impact:** high (rigor + a clean clinical-safety story judges love).
**Refs:** conformal prediction gives distribution-free guarantees (Angelopoulos 2021; Frontiers
Bioinformatics 2025); patient-level conformal AMR controls FP/FN per antibiotic (bioRxiv 2023);
selective prediction with cost-aware deferral under distribution shift (Sci Reports 2026).

### 3. **AI/LLM explanation layer** (the "AI flow") — narrative + interactive Q&A
**Why it wins:** turns raw predictions + SHAP determinants into a **clinician-readable report**
("Predicted RESISTANT to meropenem, driven by the KPC-3 carbapenemase and OmpK36 porin loss;
confidence high; recommend confirmatory testing"), and an interactive "ask why" over each call.
High demo wow-factor and on-theme (built with Claude).
**Honest scoping (important):** literature shows LLMs do **not** improve raw AMR *prediction* accuracy
(IDWeek 2025: "LLMs do not improve model performance for AMR prediction"). So we use the LLM as an
**explanation/reporting/triage-narrative layer**, NOT as the predictor — defensible and still
differentiating. **Effort:** ~1 day (Anthropic API over our prediction+SHAP JSON).
**Refs:** generative AI for explainable AMR (PMC 2025); interpretable LLMs for MDR prediction
(Health Info Sci 2025); LLM CDSS in infectious disease (systematic review 2025).

---

## TIER 2 — strong, more effort (pick if time allows)

### 4. **Novel-determinant discovery** (cheap, high novelty)
Flag features the model leans on (high SHAP) that are **absent from the known-gene rules** → candidate
novel/uncharacterised resistance markers for the microbiologist to review. Nearly free — we already
have SHAP + the rules token lists. Ties into the pan-genome "discover new resistance genes" trend
(unitig ML, PMC 2024). **Effort:** ~half day. **Impact:** medium-high novelty.

### 5. **MIC regression** (predict the value, not just R/S)
Predict minimum inhibitory concentration (the gold-standard continuous measure) instead of binary R/S;
avoids outdated breakpoints and captures resistance *degree*. BV-BRC carries `measurement` values.
**Effort:** ~2–3 days (reframe as regression, new metrics). **Impact:** high but heavier.
**Refs:** MIC-as-regression for K. pneumoniae (Bath, bioRxiv 2023); pan-genome MIC feature selection
(PMC 2023).

### 6. **Distribution-shift robustness** (temporal / geographic holdout)
Beyond lineage: train on isolates before a cutoff year (or one region) and test on the rest — the
strongest possible "it generalizes" evidence. BV-BRC has collection date/country metadata.
**Effort:** ~1 day. **Impact:** medium-high (deepens the honest-evaluation story).

---

## TIER 3 — out of scope for a laptop sprint (note, don't build)
- Protein/genomic **language-model features** (ESM/ProtBert) — GPU-heavy (Frontiers Micro 2025).
- Full **cross-species foundation model** — large compute (arXiv 2026 genomic FMs).

---

## Recommended package (best winning-% per unit effort)
**Tier-1 all three** — E. coli second organism (generalization) + conformal/selective prediction
(rigorous safety) + LLM explanation layer (AI flow, wow) — with **novel-determinant discovery** as a
near-free Tier-2 bonus. This adds: proof of generalization, statistically-guaranteed clinical safety,
and a genuine AI-assisted interface, without leaving the determinant-feature + honest-evaluation
philosophy that is already the project's strength.

**Status:** conformal ✅, LLM layer ✅, cross-resistance ✅, E. coli 🔄 (downloading).

---

# Round 2 — deeper search (what still moves the needle)

Grounded in judging-criteria and AMR-ML literature (2024–2026). Judges reward, in order:
**novelty · validation rigor on data the method wasn't built for · real-world impact with real
metrics · reproducibility · clear communication · owning your limitations.** Ranked by boost/effort.

## TIER 1 — highest ROI (strongly recommended)

### A. Temporal (prospective-style) validation — *the single biggest rigor differentiator*
**Why:** literature says only ~2/10 AMR-ML studies do prospective validation, and those show 15–30%
lower performance — reviewers explicitly look for "tested on data it wasn't designed for" and
*temporal* generalizability (evolving resistance). We already have the hardest split (lineage); adding
a **train-on-past / test-on-future** split (BV-BRC has `collection_date`) is the definitive
"it generalizes forward in time" evidence, and honestly reporting any drop is a rigor win, not a loss.
**Effort:** ~1 day. **Impact:** very high.

### B. Clinical-impact quantification — *what separates finalist from winner*
**Why:** "impact with real metrics, calculated or cited" is the finalist→winner line. Genomic
prediction returns a call in **minutes from sequence vs ~48–72 h for phenotypic AST** (cite: mNGS
models cut turnaround ~70 h). Translate our VME/ME into **stewardship terms** (a very-major error =
a missed resistant infection → wrong empiric therapy), and frame the panel against **WHO priority
pathogens** (K. pneumoniae carbapenem-R is #1 critical). **Effort:** ~half day (framing + a small
turnaround calc). **Impact:** very high, near-free.

### C. Interactive web demo (Streamlit) — *communication is a judged criterion*
**Why:** a poster/defense lands far better with a clickable app than a CLI: upload/select a genome →
per-drug R/S, calibrated confidence, **conformal "defer to lab" flag**, SHAP drivers, and the **Claude
narrative** — all on one screen. **Effort:** ~1 day. **Impact:** high (demo wow + clarity).

## TIER 2 — strong if time allows
- **Geographic cross-site validation** — country/region holdout (BV-BRC `geographic_location`);
  a second real distribution-shift test beyond lineage. ~1 day.
- **Benchmark vs an established tool (ResFinder)** — the standard external comparison in the field;
  strengthens the honest benchmark beyond our hand-rolled rules baseline. ~1–2 days.
- **Model card + dataset datasheet** — cheap, directly hits the reproducibility criterion. ~half day.
- **Active-learning / surveillance prioritization** — use conformal abstention to rank *which*
  uncertain isolates to send for phenotypic AST (a real lab-workflow use case). ~1 day.
- **MIC regression** — predict MIC values, not just R/S (gold-standard, richer). ~2–3 days.

## Recommended next package
**A (temporal validation) + B (clinical-impact quantification) + C (Streamlit demo)** — rigor,
impact, and communication, the exact three axes judges score. A and B are near-free given our
infrastructure; C is the presentation multiplier.

**Status:** A ✅, B ✅, C ✅, E. coli generalization ✅.

---

# Round 3 — deeper search (organisms + method upgrades to shine)

Grounded in the **WHO Bacterial Priority Pathogens List 2024** (24 pathogens, 15 families) and the
2025–2026 AMR-ML frontier (multi-species transfer; graph neural nets, e.g. AMR-GNN, *Nat. Commun.*
2026).

## Add organisms (ranked by "shine")
Our two organisms (K. pneumoniae, E. coli) are both **WHO critical-tier Enterobacterales**. Breadth
across families/tiers/Gram-type is the strongest generalization statement.

1. **Staphylococcus aureus — Gram-POSITIVE.** The boldest jump: Gram-negative → Gram-positive is a
   *categorical* change (cell wall, resistance biology). MRSA (mecA → oxacillin/cefoxitin) is iconic
   and WHO high-priority. Proves the method isn't Gram-negative-specific. **Highest shine.**
2. **Acinetobacter baumannii — WHO CRITICAL, non-Enterobacterales Gram-neg.** Completes the
   critical-carbapenem-resistant trio with K. pneumoniae; MDR nosocomial. High clinical relevance.
3. **Pseudomonas aeruginosa — complex intrinsic+acquired resistance;** the exact organism the SOTA
   AMR-GNN (2026) used → a direct positioning point.
4. **M. tuberculosis — mutation-only, no HGT;** a strong *contrast* (different phylogenetic structure)
   but a weaker "ML adds over rules" story (WHO mutation catalogue is strong).

**Recommendation:** add **S. aureus** (Gram-positive jump) and **A. baumannii** (critical trio) →
a 4-organism panel spanning Gram-neg Enterobacterales + Gram-neg non-Enterobacterales + Gram-positive,
across WHO critical & high tiers. Each is ~a few h download (throttled) + ~a few h annotation; the
organism-general runner already supports them (add one `ORGANISMS` entry + the right AMRFinderPlus
organism name, e.g. `Staphylococcus_aureus`, `Acinetobacter_baumannii`).

## Method upgrades (mostly free with existing data)
1. **Cross-species zero-shot transfer** (`src/evaluation/cross_species_transfer.py`, coded): train on
   one organism, test on another over shared determinants. A striking multi-species result — the exact
   frontier the literature flags — computed from data we already have. (Pending BV-BRC uptime.)
2. **Position vs SOTA (AMR-GNN / foundation models):** cite as related work; our edge is
   interpretability + honest phylogeny-aware + temporal evaluation + calibration + conformal
   uncertainty — clinical trustworthiness over raw leaderboard AUC.
3. **Acknowledge XGBoost-SHAP interpretation limits** (*J. Infection* 2025): we independently found
   the mean-|SHAP| co-selection confound and use per-instance SHAP — turn a known field pitfall into a
   demonstrated strength.

**Next package:** S. aureus + A. baumannii (breadth) + cross-species transfer (novel result) +
SOTA positioning (framing). Then poster.
