# Reading Resistance — Presentation Speaker Guide

**A complete, from-scratch walkthrough for delivering a detailed talk on the project.**

> Companion file: `slides.html` (open in any browser, press **F** for fullscreen, **P** to print to PDF).
> This guide is the *content* — every slide has a **What's on screen** line and **What to say** notes,
> plus the exact numbers, the code behind each step, and the honest positives/negatives.
> Everything here is verified against the code and result files as of the handoff.

**How to use this:** read it once end-to-end to own the story, then present from `slides.html` with this
open beside you. Sections map 1:1 to slide groups. ⚕ marks slides that are *your* home turf (microbiology).

**The one-sentence pitch:** *We read a bacterial genome and predict, per antibiotic, whether the strain
is resistant — honestly validated on lineages the model has never seen, benchmarked against a transparent
gene-lookup, and explaining the genes behind every call — across 8 WHO-priority pathogens.*

---

## PART 0 — Orientation (Slides 1–3)

### Slide 1 — Title
**On screen:** "Reading Resistance — An interpretable, honestly-benchmarked ML classifier for antibiotic
resistance from bacterial genomes." Team of two.

**What to say:** "Antibiotic resistance is projected to become a leading cause of death worldwide.
Deciding whether a strain resists a drug normally takes *days* of lab culturing. But the genome already
carries the answer — so in principle we can read it directly from sequence in minutes. That's what this
project does, and it does it *honestly* — which, as I'll show, is the hard part and the whole point."

### Slide 2 — The team & division of labor
**On screen:** two lanes. SWE (pipeline, ML, features, splitting, calibration, SHAP tooling, demo) ·
Microbiology M.Sc. (organism/drug choice, label QC, biological validation of determinants, clinical
error trade-offs, limitations & defense narrative).

**What to say:** "This is deliberately two-person. The engineering builds the rigor; the microbiology
owns every biological decision and validates that the model learned *real mechanisms*, not statistical
noise. Neither half stands alone — a leaky pipeline with perfect biology is worthless, and perfect code
predicting nonsense genes is worthless."

### Slide 3 — The 60-second story
**On screen:** 4 icons — Genome in → determinant features → per-drug R/S + confidence → the genes behind
the call. Plus the headline stat.

**What to say:** "Here's the whole system in one breath. A genome goes in. We annotate the known
resistance genes and mutations it carries. A per-drug model gives a resistant/susceptible call with a
calibrated probability. And — unlike a black box — it names the genes that drove the call. The headline
result: for **cefoxitin**, a simple gene-lookup misses 96.6% of resistant strains, because resistance
there comes from *losing* a porin, not gaining a gene. Our model lifts discrimination from a coin-flip
(AUC 0.50) to 0.90. That gap *is* the contribution."

---

## PART 1 — The Problem & Why It Matters (Slides 4–7)

### Slide 4 — Why AMR matters
**What to say:** "Bacteria acquire resistance two ways: they carry a *gene* (e.g. an enzyme that destroys
the antibiotic), or they have a *point mutation* in the drug's target. Standard testing grows the strain
against the drug — accurate but slow. Whole-genome sequencing of clinical isolates is now routine, so
huge public collections of genomes paired with lab resistance phenotypes already exist. The raw material
for prediction is sitting there."

### Slide 5 — The two pitfalls (the gap we fill) ⚕
**On screen:** Pitfall 1: clonal population structure → random splits leak. Pitfall 2: "it's just a
database lookup."

**What to say:** "Two things separate a credible project from a naive one.
**First — bacteria are clonal.** They reproduce by copying, so the dataset is full of near-twins from the
same lineage. If you split train/test randomly, the model sees a resistant strain in training and its
near-identical cousin in test, and gets rewarded for *recognizing the clone* — not for understanding
resistance. Every metric inflates. This is the bacterial version of data leakage, and it's the single
biggest way these projects fool themselves.
**Second — if a model just detects a known resistance gene, it's a database query, not science.** The
value is in what the model *adds* over that lookup."

### Slide 6 — Research question & hypothesis
**On screen:** RQ + hypothesis (a–c).

**What to say:** "The question: can an interpretable model predict resistance from genomic features,
*generalize to lineages it has never seen*, and be honestly benchmarked against a gene-lookup and
published work? Hypothesis, three parts: (a) strong per-drug discrimination on unseen lineages;
(b) match or beat the gene-lookup, especially by cutting the *dangerous* errors; (c) its top features
correspond to *established resistance mechanisms* — the biology check."

### Slide 7 — The non-negotiables (our honesty contract)
**On screen:** the 9 rules as a checklist, all ✓.

**What to say:** "We wrote down nine hard rules before touching data, and never broke them. Split by
lineage, never randomly. Always report the gene-lookup baseline. One model per antibiotic. Report the
*clinical* error rates first, not just accuracy. Determinant features only — no black-box k-mers.
Calibrate every probability. Never claim state-of-the-art on everything. Keep scope disciplined. And
never, ever fabricate a genome, a label, or a number. Every one of these is enforced in code or tests —
I'll show you the leakage test that fails the build if a single lineage crosses the split."

---

## PART 2 — Setting Up From Scratch: The Environment (Slides 8–11)

> This part answers "how would someone rebuild this from nothing?" — critical for a reproducibility grade.

### Slide 8 — The stack
**On screen:** Python 3.11 · pandas/numpy/scikit-learn · XGBoost/LightGBM · SHAP · matplotlib/scipy ·
**bioconda: AMRFinderPlus + mlst** · Streamlit · (optional) PyTorch + ESM-2.

**What to say:** "Everything is free and public. Standard Python data-science stack for the ML. Two
*bioinformatics* tools that must come from **bioconda** — that's why we use conda, not plain pip:
**AMRFinderPlus** (NCBI's resistance-gene caller — our feature engine) and **mlst** (lineage typing).
The optional protein-language-model extension adds PyTorch and ESM-2. It all runs on a laptop."

### Slide 9 — Install it in four commands
**On screen (code):**
```bash
# 1. Get the code
git clone https://github.com/kappaborg/amr-resistance-predictor
cd amr-resistance-predictor

# 2. Build the exact pinned environment (conda, because of the bio tools)
conda env create -f environment.yml
conda activate amr-resistance-predictor

# 3. Prove it works — runs the test suite incl. the zero-leakage test
make test          # -> 20 passed

# 4. Try a prediction on a cached genome (instant, no download)
python -m src.app.predict --organism saureus --genome-id 1280.10000
```
**What to say:** "Four commands from a bare machine to a working prediction. `environment.yml` pins
*exact* versions — AMRFinderPlus 4.2.7 with database 2026-05-15, mlst 2.33.1 — so anyone regenerates the
same features. `make test` runs 20 tests; the important one asserts no lineage appears in both train and
test, for all 8 organisms. If that fails, the build fails."

**Note for you:** the bio-tool paths are resolved automatically from the active conda env, so this works
on *her* laptop, Linux, or Colab — not just the original machine. (This was a portability fix made during
the final audit.)

### Slide 10 — Reproducibility scaffolding
**On screen:** fixed seed 42 · `config/config.yaml` (config-driven) · `data/manifest.md` (every source,
version, query, count) · `Makefile` (one-command stages) · `docs/decisions.md` (53 logged decisions) ·
`data/raw/` git-ignored, regenerates from manifest.

**What to say:** "Reproducibility isn't a slide, it's the spine. One fixed seed everywhere. Nothing
scientific is hard-coded — organism, drugs, thresholds, model choice all live in one config file. A data
manifest records every source, the exact BV-BRC query, tool versions, and genome counts. A Makefile
rebuilds any stage with one command. And a decision log — 53 entries — captures *every* non-trivial
choice with its rationale. That log is our lab notebook and our defense script."

### Slide 11 — Repository tour
**On screen:** the tree (`src/{data,features,split,models,evaluation,interpret,app}`, `tests/`,
`results/`, `config/`, `docs/`).

**What to say:** "Clean separation. `data/` acquires and QCs. `features/` annotates. `split/` does the
lineage split. `models/` trains. `evaluation/` measures honestly. `interpret/` explains. `app/` is the
demo. `tests/` guards it. `results/` holds every metric, figure, model, and report. You can point at any
number in the talk and trace it to the file that produced it."

---

## PART 3 — The Data (Slides 12–16)

### Slide 12 — Resources (where everything comes from)
**On screen (table):**

| Role | Source | What it gives |
|---|---|---|
| Genomes + labels | **BV-BRC** (formerly PATRIC) | Assembled genomes + lab R/S phenotypes |
| Reference determinants | **CARD** | Catalog of known resistance genes + mutations |
| Feature annotation | **AMRFinderPlus** (NCBI) | Genome → gene/mutation presence-absence |
| Lineage typing | **mlst** | Each genome → sequence type (lineage) |
| Validation reference | **Published AMR-ML papers** | Honest external comparison |

**What to say:** "Five public resources, no wet lab. BV-BRC gives genomes *paired with laboratory*
resistance phenotypes. CARD and AMRFinderPlus are the vocabulary and the tool that turns a genome into
features. mlst assigns lineages for the honest split. And published papers give us a yardstick."

### Slide 13 — Choosing the organism (data-driven, Week 1) ⚕
**On screen:** shortlist → **Klebsiella pneumoniae (taxon 573)** locked; 5-drug panel.

**What to say (your call to narrate):** "We did *not* pick the organism by gut feeling. In Week 1 we
queried BV-BRC for lab-phenotype counts across candidate pathogens and drugs, and picked for *abundant,
class-balanced* data. **K. pneumoniae** won on the science: its carbapenem resistance is exactly where
single-gene rules fail — resistance comes from combinations of porin loss *and* β-lactamases — so there
are real, catchable dangerous errors to reduce. And it has strong clonal structure (well-defined
sequence types), which makes the lineage-aware split *bite*. The 5-drug panel — meropenem, gentamicin,
ciprofloxacin, TMP-SMX, cefoxitin — spans five *distinct* resistance mechanisms, so the per-drug
explanations are genuinely different from each other."

### Slide 14 — Labels: lab phenotypes only (avoiding circularity)
**On screen:** `evidence = "Laboratory Method"` only; computational AMR calls excluded. MIC/SIR → binary
R/S; "Intermediate" dropped.

**What to say:** "A subtle but critical choice. BV-BRC contains both *laboratory* phenotypes and
*computationally predicted* ones. If we trained on computational labels, we'd be predicting another
model's predictions — circular. So we kept **only lab-measured** phenotypes. We convert the MIC/SIR
readout to binary resistant/susceptible and *drop* the ambiguous 'Intermediate' category rather than
guess. That's a scope decision flagged for microbiology sign-off."

### Slide 15 — Acquisition & the bug we caught
**On screen:** BV-BRC REST API · resumable, checksummed downloader · **caught a truncation bug**.

**What to say:** "We pull via the BV-BRC API — no giant bulk download, size flagged before pulling. Two
things worth telling because they show rigor. **One:** the genome API silently defaulted to returning
only 25 records, truncating assemblies to ~2 Mbp of a ~5.5 Mbp genome. A truncated genome *drops genes*
on the missing contigs → false-absent features → corrupted labels. A file-size integrity check caught
it; we fixed the query and re-downloaded everything. **Two:** the downloader retries with backoff and
resumes, because one network timeout shouldn't kill an overnight run. These are in the decision log as
entries #15 and #16."

### Slide 16 — Quality control (literature-aligned)
**On screen:** MIMAG (Bowers 2017): completeness ≥90%, contamination ≤5%, contigs ≤500, length 4.5–7.5
Mbp. Drops ~1.9%, class-balanced. **Final: 3,850 genomes.**

**What to say (your domain) ⚕:** "QC thresholds aren't guesses — they follow the MIMAG standard and the
Kleborate genome-size range for K. pneumoniae. We widened the length cap to 7.5 Mbp specifically because
MDR plasmid-carrying strains legitimately reach ~6.6 Mbp — a tight cap would have thrown away *real
resistant* strains. QC drops only 1.9%, and crucially the drops are class-balanced, so we don't bias the
labels. We end with **3,850 QC-passed genomes**."

---

## PART 4 — Methods: The Science Core (Slides 17–23)

### Slide 17 — Feature extraction: the determinant matrix
**On screen:** AMRFinderPlus `--organism Klebsiella_pneumoniae` → binary **3,850 × 688** matrix
(genes + point mutations; virulence/stress excluded).

**What to say:** "Every genome runs through AMRFinderPlus, which reports the known resistance genes and
established point mutations it carries. We assemble those into a binary matrix: rows = genomes, columns =
determinants, 1 = present. For K. pneumoniae that's 3,850 × **688** determinants. We keep only true AMR
elements — we exclude virulence and stress genes. We deliberately include *point mutations* like gyrA and
parC, because those are the whole story for fluoroquinolones. This is the 'fast and interpretable' feature
choice from the proposal — **not** whole-genome k-mers."

### Slide 18 — THE crown jewel: the lineage-aware split ⭐
**On screen:** MLST → sequence type per genome → `StratifiedGroupKFold` holds out *whole lineages* →
zero-overlap test. Untypeable genomes get a *unique* synthetic lineage (never pooled).

**What to say:** "This is the most important slide in the talk. We assign each genome its MLST sequence
type — its lineage. Then we split so that **entire lineages** are held out for testing: the model is
tested only on sequence types it has *never seen in training*. We use StratifiedGroupKFold so the
lineages stay disjoint *and* the resistant/susceptible balance is preserved. Genomes that can't be typed
each get their *own* unique synthetic lineage, so they can never secretly pool into one leaky group. And
a **test asserts zero lineage overlap** — if one sequence type ever appears in both train and test, the
build fails. This one design choice is what makes every number afterwards trustworthy."

### Slide 19 — Why it matters: leakage quantified
**On screen:** random vs lineage-held-out AUC. Mean inflation **+0.010**, up to **+0.032** (cefoxitin).

**What to say:** "We *measured* what the discipline buys. Scoring the same model under a naive random
split vs our lineage split, the random split inflates AUC by +0.010 on average, up to +0.032 for
cefoxitin. Now — that's a *small* inflation, and that's actually a beautiful result: it means our features
encode **mechanism, not ancestry**. A genome resists meropenem because it carries blaKPC, not because of
its family tree. K-mer models, which secretly encode lineage, lose 0.1–0.2 AUC on this same test. The
small gap is *evidence we learned biology*. But we still report the honest lineage-held-out number
everywhere — never the inflated one."

### Slide 20 — Models: one per antibiotic, three tiers
**On screen:** (1) known-gene **rules baseline**, (2) logistic regression, (3) **XGBoost** — per drug,
never pooled.

**What to say:** "Resistance is drug-specific, so it's strictly **one model per antibiotic** — we never
pool drugs into a single label. For each drug we build three things: a transparent **rules baseline**
that says 'resistant if any known gene for this drug is present' — that's the honest yardstick; a simple
**logistic regression**; and a gradient-boosted **XGBoost**. Reporting all three keeps us honest about
what the ML actually adds."

### Slide 21 — The rules baseline is organism-aware ⚕
**On screen:** `ORG_RULES` — e.g. ampicillin = β-lactamases in *E. coli* but **pbp5** in *Enterococcus*;
penicillin = **pbp1a/2b/2x** mosaic in *S. pneumoniae*.

**What to say (your domain):** "A fairness point we're proud of. 'Ampicillin resistance' does not mean
the same gene in every bug — it's β-lactamases in E. coli, but altered **pbp5** in Enterococcus. If we'd
used one generic rule everywhere, we'd be comparing our model to a *strawman* baseline in some organisms.
So the rules baseline varies by (organism, drug) to match the *actual* mechanism. That makes the
'what ML adds' comparison genuinely fair — when we beat the baseline, we beat a real one."

### Slide 22 — Calibration
**On screen:** isotonic calibration + **Brier score**; reliability curves.

**What to say:** "A model that says '90% resistant' should be right about 90% of the time. We calibrate
every probability with isotonic regression and report the **Brier score** and reliability curves. A
clinician acting on a probability needs it to *mean* something — an uncalibrated 0.9 that's really a 0.6
is dangerous."

### Slide 23 — The metric that matters: VME / ME ⚕
**On screen:** **VME** = resistant called susceptible (the *dangerous* miss — patient gets a drug that
won't work). **ME** = susceptible called resistant. Operating point tuned to cap **VME ≤ 3%** (CLSI/FDA
style).

**What to say (your domain — this is the clinical heart):** "Accuracy hides the error that kills people.
A **very-major error** is calling a resistant strain susceptible — the patient gets a drug that won't
work. A **major error** is the reverse — an unnecessary alarm. These are not symmetric: a VME can be
fatal, an ME wastes a broader antibiotic. So we report VME and ME *first*, before AUC, and we tune the
decision threshold on out-of-fold *training* data to cap VME at ≤3% — the regulatory-style target — then
measure on the untouched test lineages. The trade-off is explicit and we never hide it: capping VME
raises ME, sharply for hard-to-call drugs. That's a clinical choice, made in the open."

---

## PART 5 — Training & The Pipeline End-to-End (Slides 24–25)

### Slide 24 — The pipeline, one arrow at a time
**On screen (the workflow diagram):** Genome FASTA → QC → AMRFinderPlus → feature matrix → MLST lineage
→ lineage split (zero-leak guard) → per-drug XGBoost → isotonic calibration + VME≤3% threshold →
evaluation (VME/ME/CIs/DeLong) → SHAP → demo app.

**What to say:** "Here's the whole assembly line, and every box is a real file you can run. `make data`,
`make features`, `make split`, `make train`, `make eval`, `make figures` — or `make all` for the lot.
Two decision gates protect the result: 'are there enough resistant *and* susceptible genomes per drug?'
and 'does it generalize to unseen lineages?' If the second fails, we loop back and check the split rather
than publish an inflated number."

### Slide 25 — How a model is actually trained
**On screen:** per-drug loop → fit on train lineages → pick VME≤3% threshold on out-of-fold train preds
→ evaluate on held-out test lineages → save deployable `.joblib` bundle (model + calibrator + threshold
+ feature columns).

**What to say:** "For each drug: fit XGBoost on the training lineages, calibrate it, and choose the
operating threshold using *out-of-fold* predictions on the training data — never peeking at the test set.
Then we evaluate once on the held-out lineages. For deployment, we bundle the model, its calibrator, the
threshold, and the exact feature columns into one file per organism-drug — 36 of them — which is what the
demo loads. 8 organisms × their drug panels = **36 deployable models**."

---

## PART 6 — Results: The Payoff (Slides 26–31)

### Slide 26 — Headline: cefoxitin, where the lookup fails
**On screen:** cefoxitin rules AUC **0.50** (coin flip) → model **0.90** (pooled CV; DeLong p≈0). Driver:
**ompK36 porin loss**.

**What to say:** "This is the money slide. Cefoxitin resistance in K. pneumoniae is often driven by
*losing* a porin — ompK36 — so the drug can't get in. A gene-*lookup* is blind to a *missing* feature,
so it performs at chance: AUC 0.517, and it misses 96.6% of resistant strains (368 of 381). Our model *learns the
absence pattern* and lifts AUC to 0.90, a difference that's statistically overwhelming. This is the
cleanest possible demonstration of what machine learning adds over a database query."

### Slide 27 — K. pneumoniae 5-drug panel (unseen lineages)
**On screen (table, pooled CV / regulator-grade):**

| Drug | AUC (model) | AUC (rules) | VME | ME | DeLong p |
|---|---|---|---|---|---|
| ciprofloxacin | 0.983 | 0.911 | 2.3% | 7.7% | 7e−36 |
| TMP-SMX | 0.977 | 0.885 | 2.3% | 18.3% | 9e−51 |
| meropenem | 0.968 | 0.933 | 1.3% | 55.9% | 5e−24 |
| gentamicin | 0.981 | 0.948 | 1.6% | 32.4% | 1e−30 |
| **cefoxitin** | **0.905** | **0.500** | 0.7% | 94.2% | ~0 |

**What to say:** "The full panel on unseen lineages. Discrimination is strong everywhere — AUC 0.90 to
0.98. VME is at or below the 3% clinical bar on every drug. The DeLong test says the model beats the
gene-lookup *significantly* on all five — these p-values are astronomically small. Notice meropenem and
cefoxitin have high *major* error — that's the VME-first trade-off I mentioned: to almost never miss
resistance, the model over-calls it. We show that with pass/fail flags, we don't bury it."

### Slide 28 — The honest reading (where we DON'T win) ⚕
**On screen:** gentamicin — rules already excellent (AUC 0.948, direct AME signal), ML only ~matches.
Salmonella chloramphenicol — ML *below* rules. S. pneumoniae TMP-SMX — co-selection, not mechanism
(folA/folP not in the catalog).

**What to say (your credibility slide):** "An award isn't won by claiming you win everything — it's won
by knowing exactly where you don't. **Gentamicin:** the resistance is a direct enzyme presence, so the
gene-lookup is already excellent and we only match it — we say so. **Salmonella chloramphenicol:** the
rules baseline actually *beats* our model, reported plainly. **S. pneumoniae TMP-SMX:** the causal folA
and folP mutations aren't even in the AMRFinderPlus catalog, so any signal there is *co-selection*, not
mechanism — we label it as such rather than dress it up as an ML win. This honesty is a feature, not an
apology."

### Slide 29 — Generalization: 8 WHO-priority pathogens
**On screen:** the `generalization_heatmap.png` — 8 organisms × up-to-18 drugs, ROC 0.84–0.998. Gram−:
K. pneumoniae, E. coli, A. baumannii, P. aeruginosa, Salmonella. Gram+: S. aureus, E. faecium,
S. pneumoniae.

**What to say:** "The identical pipeline — only the taxon, the AMRFinderPlus organism string, and the
MLST scheme change — runs on **eight** priority pathogens, spanning the Gram divide and four
phylogenetic classes. Every organism-drug reaches AUC 0.84 to 0.998 on unseen lineages. And the model
earns its keep exactly where the biology is combinatorial: S. aureus cefoxitin, A. baumannii amikacin
(+0.32 over rules), E. faecium ampicillin capturing pbp5 (+0.35), S. pneumoniae penicillin capturing the
PBP mosaic. It matches the lookup on direct single-gene calls — mecA, van, carbapenemases — as it should."

**Note ⚕:** two organisms — *E. cloacae* and *C. jejuni* — were evaluated and **excluded** for
insufficient public lab data. That's honest data-driven selection: we don't build a model on ~100
genomes just to claim a bigger number.

### Slide 30 — Cross-species zero-shot transfer
**On screen:** `cipro_transfer_matrix.png` — train on one species, test on another, no retraining.
Transfers among Gram-negatives (shared gyrA/parC); collapses to Gram-positive S. aureus (uses grlA,
different gyrA numbering).

**What to say:** "A striking test of whether the model learned biology. We trained a ciprofloxacin model
on one species and tested it *zero-shot* on another — no retraining. Among the Gram-negatives it
transfers well, because they share the gyrA/parC target and residue numbering. It *collapses* to
Gram-positive S. aureus — which uses grlA and different numbering. That collapse is exactly what a
microbiologist would predict. A model that had memorized lineage artifacts couldn't produce a
mechanistically-sensible boundary like this."

### Slide 31 — Statistical rigor
**On screen:** DeLong paired test (model vs rules, every drug) · lineage-clustered bootstrap 95% CIs ·
MCC + PPV/NPV at realistic prevalence.

**What to say:** "The central claim — 'the model beats the gene-lookup' — is backed by a **DeLong paired
significance test**, not just a visible gap. Confidence intervals come from a **lineage-clustered**
bootstrap that resamples whole lineages, not genomes, because the resistant class is small and naive CIs
would be too optimistic. And because the resistant class is often the minority, we add MCC and predictive
values at realistic local prevalences — accuracy alone would flatter us."

---

## PART 7 — Advanced & Novel Extensions (Slides 32–36)

### Slide 32 — MIC regression: predicting the *level*
**On screen:** continuous MIC, scored by **Essential Agreement** (within ±1 doubling dilution, the
CLSI/FDA criterion). K. pneumoniae: cipro **92.7%**, gentamicin 85.1%, cefoxitin 79.9%, meropenem 75.4%.

**What to say:** "Binary R/S is coarse. Clinicians think in **MIC** — the minimum concentration that
inhibits growth. So we also predict the continuous MIC and score it the way regulators actually do:
**Essential Agreement**, prediction within one two-fold dilution of the truth. Ciprofloxacin hits 92.7%,
right at the clinical expectation. This is a richer, quantitative output — not just *whether* a strain
resists, but *how strongly*. We also made the scoring **censor-aware**: most real MICs are reported as
'>32' or '≤0.5' bounds, and naively clipping them distorts the metric."

### Slide 33 — ESM-2 protein language model (honest positive)
**On screen:** presence/absence discards *which allele* — but carbapenem MIC depends on it (KPC-2 vs
KPC-3, OXA/NDM/VIM). Embed each resistance protein with ESM-2. **Meropenem EA 68→73% (+4.8, p=0.0004);
A. baumannii imipenem 54→60% (+6.7, p=0.0001).** Neutral elsewhere; 650M gave no extra gain over 150M.

**What to say:** "Here's a genuinely novel extension. Presence/absence throws away *which variant* of a
gene a strain carries — but for carbapenems, the specific allele matters a lot: KPC-2 vs KPC-3, or the
OXA/NDM/VIM families, behave differently. So we extract each resistance protein's actual sequence and
embed it with **ESM-2**, a protein language model, so different alleles get different vectors. On exactly
the hard carbapenem drugs it gives a real, significant gain — meropenem Essential Agreement 68 to 73%,
A. baumannii imipenem 54 to 60% — and it's *neutral* on drugs whose determinants already saturate, which
is precisely what the biology predicts. And scaling to the bigger 650M model gave *no* further gain, so
we honestly use the laptop-scale 150M."

### Slide 34 — GNN: honest negative (we tested SOTA and it lost)
**On screen:** isolate-similarity Graph Neural Network (AMR-GNN style) — **worse than a plain MLP on all
5 drugs.** Reported as a negative, with the reason.

**What to say:** "This is my favorite slide for a defense, because it's a *negative* result. There's a
fashionable graph-neural-network approach where each genome is a node linked by similarity. We
implemented it properly — the actual published architecture — on identical lineage-grouped folds. It was
*worse* than a plain neural net on our features, for every drug. Why? Because our curated determinant
features already capture the mechanistic signal cleanly; a population-structure graph only re-adds the
lineage information we're deliberately trying *not* to lean on, and its smoothing blurs sharp gene→
phenotype boundaries. So we can say: 'we implemented the state-of-the-art graph model, it did not win,
and here is exactly why.' That's a stronger position than assuming it would help."

### Slide 35 — Uncertainty & abstention (conformal prediction)
**On screen:** Mondrian conformal → per-class coverage + an explicit **abstain → defer to phenotypic
testing** option. Confidently classifies ~95% for well-powered drugs; for hard drugs makes confident
calls on ~55–58% with near-zero VME and defers the rest.

**What to say:** "A responsible clinical tool must know when it doesn't know. Conformal prediction gives
each call a validated confidence level and, crucially, an **abstain** option — 'I'm not sure, send this
one to the lab.' For well-powered drugs it confidently handles ~95% of strains; for the hard carbapenems
it confidently calls ~55% with *near-zero* very-major error and defers the rest. And because we test on
*unseen lineages* — breaking the assumption conformal usually relies on — the coverage is *empirically
validated under distribution shift*, which is stricter than the textbook version."

### Slide 36 — Clinical utility: decision curves & risk-coverage
**On screen:** net-benefit vs threshold (model beats rules + treat-all across the VME-averse region);
deferring the least-confident ~30% roughly *halves* VME (meropenem 8.7%→3.5%).

**What to say:** "Two analyses most AMR-ML papers skip. A **decision curve** ties the model to clinical
*net benefit* — and across the low-threshold region where a clinician who fears missing resistance
operates, our model beats both the gene-lookup and 'treat everyone.' And a **risk-coverage** analysis
shows that deferring the least-confident 30% of calls to the lab roughly halves the very-major error
rate. The model knows when it doesn't know, and that translates into measurable clinical value."

---

## PART 8 — Interpretability & Biological Validation (Slides 37–38) ⚕ YOUR SECTION

### Slide 37 — SHAP: the genes behind every call
**On screen:** SHAP beeswarm figures per drug; top determinants.

**What to say (your domain):** "Every prediction is auditable. We use SHAP to rank which determinants
drove each call, globally and per strain. This is where a black box becomes a colleague you can argue
with — if the model leans on a gene that makes no biological sense, we catch it."

### Slide 38 — Do the top features match known biology? (the validation table) ⚕
**On screen (table):**

| Drug | Model's top determinants | Known mechanism | Match? |
|---|---|---|---|
| ciprofloxacin | gyrA, parC (QRDR) + qnr | Target mutation in DNA gyrase/topoisomerase | ✓ |
| gentamicin | aac(3) acetyltransferases | Aminoglycoside-modifying enzyme | ✓ |
| TMP-SMX | sul1, dfrA | Dihydropteroate synthase / DHFR bypass | ✓ |
| meropenem | blaKPC-2, blaKPC-3 | Carbapenemase | ✓ |
| cefoxitin | ompK36 porin loss | Reduced drug entry (loss, not gain) | ✓ |
| S. aureus oxacillin | mecA | Altered PBP2a | ✓ |
| E. faecium ampicillin | pbp5 | Altered penicillin-binding protein | ✓ |
| S. pneumoniae penicillin | pbp1a/2b/2x mosaic | Altered PBPs (mosaic) | ✓ |

**What to say (your headline):** "This is the biological validation — and it's clean. For every drug, the
determinants the model relied on are the *textbook* mechanisms. Fluoroquinolones: gyrA/parC target
mutations. Aminoglycosides: modifying enzymes. Carbapenems: the KPC carbapenemase, and after we enlarged
the dataset, KPC-2 and KPC-3 rose to the *top two* — the added lineage diversity broke an earlier
gene-lineage confound. And cefoxitin correctly keys on porin *loss*. The model learned real microbiology.
**Honest caveat:** co-selection still inflates some co-carried markers — resistance genes travel together
on the same MDR plasmids — so per-strain SHAP, not just the average ranking, is the honest lens. And
because our features are *catalogued* determinants, none of this is novel-gene discovery — that would need
pan-genome features, which are out of scope."

---

## PART 9 — The Demo (Slides 39–40)

### Slide 39 — The interactive app
**On screen:** Streamlit, 8 organisms. Pick organism → cached genome or **upload a FASTA** → per-drug
R/S + calibrated P(resistant) + SHAP drivers. Optional AI clinical narrative.

**What to say:** "It's a working web app, not a mockup. Pick one of the eight organisms, choose a cached
genome or upload your own FASTA, and it runs AMRFinderPlus live and returns per-drug calls with
calibrated probabilities and the genes behind each. There's an optional AI narrative that turns the
output into clinician-readable prose — but it's an *explanation layer only*, it can never change a call.
With no API key it falls back to a deterministic template, so the demo always works."

### Slide 40 — Safety: the species guard
**On screen:** MLST scheme check + intrinsic-marker genes → withholds predictions on wrong-species /
non-bacterial uploads. Every screen: "research/decision-support only — not a diagnostic."

**What to say:** "A safety feature we're proud of. The model is a function of the determinant vector — it
will happily 'predict' on *anything*. Early on, a fungal protein file produced a spurious call. So the app
now verifies the species two independent ways — the MLST scheme *and* intrinsic marker genes — and
**withholds** predictions if the upload doesn't match the selected organism. And every screen repeats the
disclaimer: this is decision-support, it complements phenotypic testing, it does not replace it."

**Live demo script (if you run it):**
```bash
conda activate amr-resistance-predictor
python -m streamlit run src/app/streamlit_app.py
# 1. Organism = S. aureus, genome 1280.10000 -> mecA drives oxacillin RESISTANT (P=0.94)
# 2. Organism = K. pneumoniae, a carbapenem-resistant genome -> blaKPC drives meropenem RESISTANT
# 3. Show the "drivers" line under each call — that's the interpretability
```

---

## PART 10 — Positives & Negatives, Honestly (Slides 41–42)

### Slide 41 — Strengths (the positives)
**On screen:**
- **Honest evaluation by construction** — lineage-held-out split, enforced by a test that fails the build.
- **Beats the gene-lookup where it matters**, proven with DeLong significance (cefoxitin, TMP-SMX,
  pbp-driven β-lactams, amikacin).
- **Generalizes** across 8 pathogens and the Gram divide; transfers zero-shot along mechanistic lines.
- **Interpretable** — every call names its genes, validated against known biology.
- **Clinically framed** — VME-first, calibrated, with abstention and net-benefit analysis.
- **Beyond binary** — MIC regression + ESM-2 allele resolution for the hard carbapenems.
- **Rigorous negatives** — GNN tested and rejected with a mechanistic reason.
- **Reproducible** — pinned env, seeds, manifest, 53-entry decision log, 20 passing tests.

**What to say:** "Pick any three of these to emphasize depending on the audience. For a science panel,
lead with the honest split and the significance testing. For a clinical audience, lead with VME-first and
the abstention option."

### Slide 42 — Limitations (the negatives — say them proudly) ⚕
**On screen:**
- **Scope:** 8 organisms, ~3–5 drugs each; K. pneumoniae is the deep-dive. Not a universal tool.
- **Co-selection confounding** inflates some co-carried determinants (reduced, not eliminated).
- **High major-error cost at VME≤3%** for hard drugs (cefoxitin, meropenem) — a threshold trade-off,
  not a discrimination failure.
- **Reference-database bias** — AMRFinderPlus/CARD/BV-BRC over-represent well-sequenced regions/pathogens.
- **Determinant features only** — no novel-gene discovery (would need pan-genome/k-mer features).
- **A few honestly weak spots** — Salmonella chloramphenicol (rules > ML), A. baumannii carbapenems
  degenerate at strict VME, S. pneumoniae TMP-SMX is co-selection-driven, A. baumannii meropenem MIC gain
  borderline (p=0.051, under-powered).
- **Published comparison is band-level**, not a re-run head-to-head.

**What to say:** "These are on the slide *on purpose*. Every one of them is a question a judge would ask
— so we ask it first, and answer it. A limitation you raise yourself reads as mastery; the same one dragged
out of you reads as a gap. This slide is where the honesty thesis of the whole project pays off."

---

## PART 11 — Ethics & Reproducibility (Slide 43)

### Slide 43 — Ethics & reproducibility
**On screen:** Not a diagnostic device · complements phenotypic testing · public de-identified data only ·
interpretability as an ethical safeguard · full reproducibility pack.

**What to say:** "Ethically: this is a research and decision-support aid, **not a diagnostic device**, and
it does not replace phenotypic susceptibility testing — a missed resistant call has real cost, which is
exactly why we lead with very-major error. Only public, de-identified data. We state the database bias
plainly. And interpretability is itself an ethical safeguard — every call is auditable rather than opaque.
Reproducibility: pinned environment, fixed seed, data manifest, decision log, and a test suite anyone can
run. This is designed to be *rebuilt and checked*, not taken on trust."

---

## PART 12 — Defense Q&A Prep (reference — not slides)

Anticipated questions and crisp answers:

**Q: Isn't this just a database lookup?**
A: For some drugs the lookup is genuinely strong and we say so (gentamicin). But where resistance is
combinatorial or from *loss* of function — cefoxitin porin loss, pbp-driven β-lactams — the lookup is at
chance (AUC 0.50) and the model reaches 0.90, DeLong p≈0. The contribution is the measured *gap*.

**Q: How do you know you're not leaking population structure?**
A: We split by MLST lineage, hold out whole sequence types, and a test fails the build on any overlap. We
also *quantified* the leakage: a random split inflates AUC by only +0.010 on average, because our features
encode mechanism, not ancestry.

**Q: Why only determinant features? Wouldn't k-mers do better?**
A: On a random split they look better — and then lose 0.1–0.2 AUC on a lineage split because they encode
lineage. We chose interpretable, mechanism-based features on purpose. We also tested a graph model that
would exploit structure; it lost.

**Q: Your major-error rates are high on some drugs.**
A: That's a deliberate operating-point choice — we cap the *dangerous* very-major error at ≤3%, which
necessarily raises major errors on hard-to-call drugs. The threshold-free AUC/PR-AUC show the underlying
discrimination is strong, and the conformal/abstention analysis recovers a low-error high-confidence
subset.

**Q: How does this compare to published methods?**
A: Our AUCs (0.84–0.998) sit in or above the published band for these organisms — but measured under a
*stricter* lineage-held-out split, plus interpretability and clinical-error reporting most studies omit.
We don't claim state-of-the-art; we claim honestly-measured.

**Q: What's the single most important thing you did?**
A: The lineage-aware split enforced by a test. Everything downstream is only trustworthy because of it.

**Q (for the microbiologist): Did the model learn real biology?**
A: Yes — the top determinant for every drug is the textbook mechanism (gyrA/parC, aac(3), sul/dfr, KPC,
ompK36, mecA, pbp5, PBP mosaics), and cross-species transfer follows mechanistic lines. Caveat:
co-selection still inflates some co-carried plasmid markers, so we read per-strain SHAP, not just averages.

---

## Appendix — Numbers cheat-sheet (memorize these)

- **3,850** QC-passed K. pneumoniae genomes · **688** determinants · **5** drugs (deep-dive)
- **8** WHO-priority pathogens · **36** deployable models · ROC-AUC **0.84–0.998** on unseen lineages
- Cefoxitin: rules **0.50** → model **0.90** (the headline)
- E. faecium ampicillin: rules 0.69 → model 0.99 (+0.35, captures pbp5)
- Leakage inflation from a random split: mean **+0.010**, max **+0.032**
- MIC Essential Agreement (K. pneumoniae): cipro **92.7%**, gentamicin 85.1%, cefoxitin 79.9%, mero 75.4%
- ESM-2 meropenem EA **68→73%** (p=0.0004); A. baumannii imipenem **54→60%** (p=0.0001)
- GNN: worse than a plain MLP on **all 5** drugs (honest negative)
- **20** passing tests (incl. zero-lineage-leakage, all 8 organisms) · **53** logged decisions
- Excluded for insufficient data: **E. cloacae, C. jejuni** (honest, data-driven)

*Reference files: `results/reports/REPORT.md` (full write-up), `docs/decisions.md` (lab notebook),
`data/manifest.md` (provenance), `results/reports/summary_*.md` (per-phase detail),
`results/poster/poster.html` (A0 poster).*
