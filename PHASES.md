# Reading Resistance — Phase Plan

> Working protocol: **every phase ends with a written summary** (what was done, key numbers,
> decisions logged, what's blocked, what's next) before the next phase begins. Phases marked
> **GATE** require explicit sign-off; ⚕ marks microbiology decisions owned by the co-builder.
> Summaries land in `results/reports/`; decisions in `docs/decisions.md`.

Maps onto the 4-week sprint: Phases 0–6 = Week 1 · 7–8 = Week 2 · 9–11 = Week 3 · 12 = Week 4.

---

## Phase 0 — Project setup & scaffold  (Week 1)
**Goal:** a reproducible skeleton before any science.
- Rename brief to `CLAUDE.md` at repo root; move docx + 3 diagrams into `proposal/`.
- `git init`; create tree: `src/{data,features,split,models,evaluation,interpret,app}`,
  `config/`, `data/{raw,interim,processed}`, `notebooks/`, `tests/`, `results/{figures,metrics,models,reports}`, `docs/`.
- Author `environment.yml` (Python 3.11+, pandas, numpy, scikit-learn, xgboost/lightgbm, shap,
  matplotlib; bioconda: ncbi-amrfinderplus or rgi, mlst), `config/config.yaml`, `Makefile`,
  `README.md`, `data/manifest.md`, `docs/decisions.md`.
- **Exit:** environment builds clean; `make` targets stubbed; repo committed.
- **Summary #0** → scaffold state, env versions, any tooling issues.

## Phase 1 — Data survey & organism/drug selection  (Week 1)  ⚕ GATE
**Goal:** pick the target from the data, not blind.
- Write BV-BRC query plan: how to count genomes with R/S phenotypes per organism × antibiotic,
  filtering for **abundant + class-balanced** labels.
- Run counts (metadata only, no bulk download). Produce a shortlist table (organism, drug,
  #R, #S, balance, est. download size).
- **⚕ Decision (co-builder):** choose one organism + a handful of drugs.
- **Exit / GATE:** organism + drug list approved; download size flagged and approved.
- **Summary #1** → shortlist table, the chosen target + rationale, logged decision.

## Phase 2 — Data acquisition  (Week 1)
**Goal:** get exactly the chosen data, versioned.
- Pull assembled genomes + phenotypes for the chosen organism/drugs via BV-BRC CLI/API.
- Record every source, query, version/date, checksum in `data/manifest.md`. `data/raw/` git-ignored.
- **Exit:** raw genomes + phenotype table on disk, manifest complete.
- **Summary #2** → counts pulled, sizes, manifest entries, any access issues.

## Phase 3 — QC & label harmonization  (Week 1)  ⚕
**Goal:** clean genomes, trustworthy binary labels.
- Genome quality filter (completeness / contamination thresholds).
- **⚕ MIC/SIR → binary R/S** mapping; decide how "intermediate" is handled; co-builder audits labels.
- **Exit:** QC'd genome set + harmonized per-drug R/S label table; drop list documented.
- **Summary #3** → genomes kept/dropped, label mapping rules, per-drug R/S counts post-QC.

## Phase 4 — Feature extraction  (Week 1)
**Goal:** the interpretable determinant matrix.
- Run AMRFinderPlus / RGI per genome → gene + point-mutation presence/absence.
- Assemble genome × determinant binary matrix; join to labels; write a builder test.
- **Exit:** feature matrix built; feature-builder test passes.
- **Summary #4** → matrix shape, #determinants, sparsity, join integrity.

## Phase 5 — Lineage assignment & phylogeny-aware split  (Week 1)  ⚕
**Goal:** the project's decisive rigor step.
- **⚕ Choose lineage method** (MLST vs. phylo cluster) + granularity; assign each genome a lineage.
- Build the phylogeny-aware train/test split; **write the leakage test asserting zero shared lineage**.
- **Exit:** split produced; leakage test passes.
- **Summary #5** → lineage method, #lineages, train/test sizes, leakage-test result.

## Phase 6 — Thin end-to-end slice, ONE drug  (Week 1)  GATE — GO/NO-GO
**Goal:** ugly-but-complete pipeline proving the whole chain works.
- Single drug: features → lineage split → logistic-regression baseline → metrics (ROC-AUC,
  PR-AUC, VME, ME) → rules-baseline comparison.
- **Exit / GATE:** working baseline for one drug exists. *If data/feature pipeline is badly
  stuck → invoke BRCA fallback.*
- **Summary #6** → first per-drug numbers vs. rules baseline; explicit GO/NO-GO recommendation.

## Phase 7 — Full multi-drug modeling  (Week 2)
**Goal:** real models across all chosen drugs.
- Feature matrix across all drugs; per-drug phylogeny-aware splits.
- Per drug: baseline → gradient boosting (XGBoost/LightGBM); hyperparameter tuning.
- **Exit:** trained per-drug models for every kept drug.
- **Summary #7** → per-drug model set, tuning ranges, raw discrimination numbers.

## Phase 8 — Calibration & clinical metrics  (Week 2)  GATE
**Goal:** honest, calibrated, clinically-framed results.
- Probability calibration (isotonic/Platt) per drug; reliability curves.
- Full metric table per drug: ROC-AUC, PR-AUC, **VME, ME** first.
- **Exit / GATE:** validated multi-drug models with calibration + clinical errors.
- **Summary #8** → per-drug metric table (VME/ME lead), calibration quality.

## Phase 9 — Honest benchmark  (Week 3)
**Goal:** measured added value, not a headline.
- Compare each model vs. known-gene rules baseline **and** published numbers on the same organism.
- Emphasize what the model adds — especially reduced very-major errors.
- **Exit:** honest benchmark table complete.
- **Summary #9** → benchmark table + one-line "what the model adds" per drug.

## Phase 10 — Interpretability & biological validation  (Week 3)  ⚕
**Goal:** every call explained and biologically checked.
- Global SHAP per drug → top determinants.
- **⚕ Co-builder validates** top features against known mechanisms; flag artifacts vs. leads.
- **Exit:** SHAP rankings + biological validation notes per drug.
- **Summary #10** → top determinants per drug, agreement with known biology, surprises.

## Phase 11 — Demo application  (Week 3)  GATE
**Goal:** genome in → explained prediction out.
- Minimal interface: input genome → per-drug R/S + calibrated confidence + determinants
  ("why this strain resists drug X").
- **Exit / GATE:** results + working demo exist.
- **Summary #11** → demo capabilities, example run, limitations.

## Phase 12 — Deliver  (Week 4)  GATE
**Goal:** submission-ready package.
- Report/paper write-up (methods+results appended incrementally throughout), poster, slides.
- Limitations + ethics; reproducibility pack (pinned env, manifest, one-command build); defense Q&A prep.
- **Exit / GATE:** submission-ready.
- **Summary #12** → final deliverable checklist, reproducibility verification, open limitations.

---

### Cross-cutting (all phases)
- Draft the write-up **incrementally** — never save it for Week 4.
- Log every non-trivial scientific choice in `docs/decisions.md`.
- Surface bad news early (leaky split, unobtainable data, too-few-labels) — immediately, not papered over.
- Never fabricate; stop before any large download/login/heavy compute and flag size first.
