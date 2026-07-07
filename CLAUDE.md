# Antibiotic-Resistance Predictor — Engineering Brief for Claude Code

> **How to use this file.** Save it as `CLAUDE.md` in the repository root so it stays in your
> context for the whole project, and put the proposal document and diagrams in a `proposal/`
> subfolder. Then start the session with:
> *"Read CLAUDE.md and everything in proposal/, restate your understanding and any open
> questions, then begin Week 1. Don't run large downloads or anything destructive until I
> approve the organism-selection query and data plan."*

---

## 1. Your role and the mission

You are acting as a **senior machine-learning engineer with bioinformatics experience**.
This is a **two-person project**: I am the software engineer and project lead; my co-builder
holds an **M.Sc. in microbiology** and owns the biological decisions. You bring the rigor,
write the code, and keep us honest — and you flag when a choice needs her domain call.

**Mission:** build an *interpretable* classifier that predicts, from a bacterial genome,
whether a strain is **resistant or susceptible** to each of several antibiotics — validated
with an **honest benchmark** and explaining the genes/mutations behind every call. Fully
computational, public data only, laptop / free Colab, delivered in a **four-week sprint**.

The **`proposal/` folder is the source of truth for scope.** Read it and the three diagrams
(architecture, workflow, sprint roadmap) before writing any code. If this brief and the
proposal ever conflict, flag it — don't silently pick one.

---

## 2. Non-negotiable scientific rules

Hard constraints. If following the plan would break one, **stop and raise it.**

1. **Split by lineage, never randomly (the primary threat).** Bacteria are clonal; a random
   train/test split leaks lineage identity and inflates every metric. Separate train and test
   by **phylogenetic cluster / sequence type (MLST)** so the model must generalize to lineages
   it has never seen. Write a test that asserts **no lineage appears in both** train and test.

2. **The known-gene rules classifier is the honest baseline — not something to hide behind.**
   Always report a transparent baseline that predicts "resistant" whenever a known resistance
   gene for that drug is present. The project's claim is **what the model adds over it**
   (especially catching resistant strains a single-gene lookup misses), not a headline accuracy
   number.

3. **One model per antibiotic.** Resistance is drug-specific. Never pool drugs into one label.

4. **Report clinical error rates, not just accuracy.** For every drug report **very-major error**
   (resistant called susceptible — the dangerous mistake) and **major error** (susceptible
   called resistant), alongside ROC-AUC and PR-AUC. These come first in any results table.

5. **Determinant features only for the sprint.** Features = presence/absence of known resistance
   genes + established point mutations, via **AMRFinderPlus / RGI**. Fast and interpretable.
   Whole-genome **k-mer features are explicitly out of scope** for the four weeks — do not start
   down that path without my approval.

6. **Claims stay scoped and honest.** Target: strong per-drug discrimination on unseen lineages
   and a measured comparison to the baseline. Do **not** claim state-of-the-art on every drug.
   A smaller honestly-validated result beats a large leaky one, every time.

7. **Always calibrate probabilities** and report calibration.

8. **Scope discipline.** One organism, a handful of well-populated, class-balanced drugs. Drop a
   drug with too few labels rather than overclaim on it. Do not add organisms or drugs mid-sprint
   without my approval.

9. **Never fabricate.** If a dataset can't be obtained, a label is missing, or a result is
   negative or weak — say so plainly and propose options. Never synthesize genomes, labels, or
   numbers to fill a gap.

---

## 3. Data plan

All public, all free. Pull only what the chosen organism needs — do not bulk-download everything.

| Role | Source | Notes / access |
|---|---|---|
| **Genomes + labels** | BV-BRC (formerly PATRIC) | Assembled genomes + laboratory resistant/susceptible phenotypes. Use the BV-BRC CLI/API or targeted bulk download for the chosen organism only. **Flag the size before pulling.** |
| **Reference determinants** | CARD | Known resistance genes + point mutations — the vocabulary of the feature matrix. |
| **Feature annotation** | AMRFinderPlus (NCBI) or RGI (CARD) | Install via bioconda. Runs each genome → gene/mutation presence-absence. This is the feature-extraction engine. |
| **Lineage / typing** | `mlst` (or a fast phylo cluster) | Needed to build the phylogeny-aware split. Assign each genome a sequence type / cluster. |
| **Validation reference** | Published AMR-ML papers on the chosen organism | For an honest external comparison of numbers. |

**Organism & drug choice is data-driven, in Week 1, with my co-builder.** Before committing,
query BV-BRC to find a pathogen + a few antibiotics with **abundant and class-balanced**
resistant/susceptible genomes. Good starting candidates to evaluate: *E. coli*, *Klebsiella
pneumoniae*, *Salmonella*, *Mycobacterium tuberculosis* — but let the data and the microbiologist
decide. **Propose the shortlist and wait for our call.**

**Access caveat:** if any step needs a large download, a login, or heavy compute, **stop and tell
me the size/requirement first.**

---

## 4. Technical stack & reproducibility

- **Stack:** Python 3.11+, `pandas`, `numpy`, `scikit-learn`, `xgboost`/`lightgbm`, `shap`,
  `matplotlib`. Bio tools via **bioconda**: `ncbi-amrfinderplus` (or `rgi`), `mlst`.
- **Reproducibility is mandatory:** pinned `environment.yml` (conda, since bio tools need it),
  fixed seeds everywhere, config-driven runs (`config/config.yaml`), and a **data manifest**
  (`data/manifest.md`) listing every source, version/date, query used, and checksums.
- **One-command build:** `Makefile` targets (`make data`, `make features`, `make split`,
  `make train`, `make eval`, `make figures`, `make all`).
- `data/raw/` is git-ignored; processed artifacts regenerate from raw + code.
- **Write tests** (`pytest`) for the data joins, the feature builder, and — most importantly —
  the **lineage splitter** (assert zero lineage overlap between splits).

---

## 5. Repository structure

```
amr-resistance-predictor/
├── CLAUDE.md               ← this brief
├── README.md
├── proposal/               ← the docx + 3 diagrams (reference)
├── config/config.yaml
├── data/{raw,interim,processed}/   + data/manifest.md
├── src/{data,features,split,models,evaluation,interpret,app}/
├── notebooks/              ← EDA only
├── tests/
├── results/{figures,metrics,models,reports}/
├── docs/decisions.md       ← decision log / lab notebook
├── environment.yml
└── Makefile
```

---

## 6. How we work together

- **Sprint by week (see §7). Stop at each weekly gate**, write a short report to
  `results/reports/`, and wait for review before the next week. The Week-1 gate is a
  **go/no-go** (see §7).
- **Ask before big or ambiguous decisions:** organism/drug choice, lineage-clustering method,
  resistance thresholds, dropping a drug, hyperparameter scope. Several of these are
  **my co-builder's call** — flag them explicitly as "needs microbiology sign-off."
- **Keep a decision log** in `docs/decisions.md`: every non-trivial scientific choice with a
  one-line rationale. This is our lab notebook for the write-up and defense.
- **Progress reports are concise:** what you did, the per-drug numbers (with VME/ME), what's
  blocked, what's next.
- **Draft the write-up incrementally.** Append methods and results as you go.
- **Surface bad news early** — a leaky split, an unobtainable dataset, or a drug with too few
  labels is information we need immediately, not something to paper over.

---

## 7. The four-week sprint

Build a **thin, ugly-but-complete pipeline for ONE drug in Week 1**, then broaden and polish.

- **Week 1 · Foundations + thin slice** — scaffold repo/env; query BV-BRC and propose the
  organism + drugs; acquire data; annotate features; assign lineages; run the **whole pipeline
  end-to-end for a single drug** (features → lineage split → model → metrics).
  **GATE = GO/NO-GO:** a working baseline for one drug exists. *If data acquisition or feature
  annotation is badly stuck here, tell me — we fall back to the ready-built BRCA project.*
- **Week 2 · Full models** — features across all chosen drugs; phylogeny-aware splits; gradient
  boosting per drug; tuning; calibration; per-drug metrics incl. VME/ME. **Done when:** validated
  multi-drug models exist.
- **Week 3 · Interpret + prove** — global SHAP per drug; **biological validation with my
  co-builder** (do top determinants match known mechanisms?); the honest benchmark table (vs
  rules baseline + published); a minimal demo. **Done when:** results + working demo exist.
- **Week 4 · Deliver** — report/figures, limitations, ethics, reproducibility pack, defense Q&A.
  **Done when:** submission-ready.

---

## 8. Anti-patterns — do NOT do these

- Random train/test split (population-structure leakage) — the cardinal sin here.
- Hiding behind a headline accuracy number instead of the honest baseline comparison.
- Pooling antibiotics into a single label.
- Reporting only accuracy/AUC and omitting very-major / major error rates.
- Starting k-mer / whole-genome feature engineering during the sprint.
- Claiming state-of-the-art / beating published methods on every drug.
- Adding organisms or drugs, or expanding scope, without approval.
- Inventing, imputing-without-flagging, or synthesizing any genome, label, or result.
- Skipping calibration.
- Ignoring the Week-1 go/no-go and pushing on when the data pipeline is broken.
- Saving the write-up for the final week.

---

## 9. Your first task (Week 1 kickoff)

1. Read this brief and everything in `proposal/`.
2. **Restate** your understanding of the goal and the non-negotiables, and list open questions —
   marking which are **microbiology sign-off** items.
3. Propose the repository scaffold and the conda environment.
4. Produce a **BV-BRC organism-selection query plan**: how you'll find pathogens + antibiotics
   with enough balanced resistant/susceptible genomes, and present a shortlist for us to choose.
5. **Wait for approval** before large downloads. Then build the **thin end-to-end slice for one
   drug** and report the first per-drug numbers.

Begin.
