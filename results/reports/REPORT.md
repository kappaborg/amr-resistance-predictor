# Reading Resistance — Project Report

*An interpretable, honestly-benchmarked classifier for antibiotic resistance from bacterial genomes —
developed on* Klebsiella pneumoniae *and generalized across eight WHO-priority pathogens.*

> **Status:** complete for the sprint scope; maintained as a living document. Primary results are on
> the full **3,850-genome** *K. pneumoniae* dataset (688 determinants); the identical pipeline runs on
> seven further WHO-priority pathogens (§4.2). Statistical rigor (DeLong, bootstrap CIs, quantified
> leakage) in §4.3; protein-language-model extension in §5.3. Reporting is mapped to the **TRIPOD+AI**
> checklist in Appendix A. Sections marked ⚕ carry claims that are the microbiologist's call and are
> flagged for sign-off. **Known open limitation:** all data derive from a single aggregator (BV-BRC),
> so no *externally-sourced* validation cohort has been run — see §7.

---

## Abstract
We predict, from an assembled *K. pneumoniae* genome, whether a strain is resistant or susceptible
to each of five antibiotics — meropenem, gentamicin, ciprofloxacin, trimethoprim-sulfamethoxazole,
and cefoxitin — using interpretable resistance-determinant features. Every model is validated under a
**phylogeny-aware split** (train and test share no MLST sequence type), reports **clinical error
rates** (very-major / major) alongside ROC-AUC/PR-AUC, and is benchmarked against a transparent
known-gene rules baseline. On unseen lineages the models achieve **ROC-AUC 0.91–0.98**. The headline
result is **cefoxitin**, where the gene-lookup baseline is effectively a coin flip (ROC-AUC 0.517;
it misses **96.6%** of resistant strains — 368 of 381) because resistance is driven by porin loss the
lookup cannot see, while the model reaches **ROC-AUC 0.906** (0.935 for calibrated XGBoost) and cuts
the very-major error rate from 96.6% to 1.8% — a concrete demonstration of what machine learning adds
over a database query.

## 1. Problem & contribution
Genome-based resistance prediction is fast but two pitfalls separate a credible result from a naive
one, and both are central here: (i) bacteria are clonal, so a random train/test split leaks lineage
identity and inflates every metric; (ii) if a model just detects a known resistance gene it is a
database lookup, not a contribution. We address both with a phylogeny-aware evaluation and an explicit
rules baseline, and we report the clinically meaningful errors (missed resistance) that a headline
accuracy hides.

## 2. Data
- **Genomes + phenotypes:** BV-BRC (`genome_amr`, `genome_sequence`), *K. pneumoniae* (taxon 573),
  **laboratory phenotypes only** (`evidence = "Laboratory Method"`) — computational predictions are
  excluded to avoid circular labels.
- **Full set:** **3,850 QC-passed genomes** (initial 1,472 + Phase-2b top-up of ~2,378 additional
  genomes, targeted at resistant strains for meropenem/gentamicin/cefoxitin). Per-drug resistant
  counts rose ~3–5× (e.g. meropenem 294→1,427 R, cefoxitin 319→1,525 R).
- **QC (literature-aligned, MIMAG):** completeness ≥90%, contamination ≤5% (Bowers et al. 2017),
  contigs ≤500, genome length 4.5–7.5 Mbp (Kleborate range). Drops 1.9%, class-balanced.
- Provenance in `data/manifest.md`; every choice in `docs/decisions.md`.

## 3. Methods
- **Features:** AMRFinderPlus (`--organism Klebsiella_pneumoniae`, DB 2026-05-15.1) → binary
  genome × determinant matrix (acquired genes + point mutations; VIRULENCE/STRESS excluded).
  Matrix: **3,850 × 688 determinants**.
- **Lineage-aware split:** MLST (`mlst`, klebsiella scheme) assigns each genome a sequence type
  (400 STs); untypeable genomes get a unique synthetic lineage. **StratifiedGroupKFold** holds out
  whole STs for the test set while balancing R/S. A test asserts zero shared lineage
  (`tests/test_split.py`).
- **Models (one per antibiotic):** transparent known-gene **rules baseline** (drug-specific
  determinants), **logistic regression**, and **XGBoost** with **isotonic calibration**.
- **Hyperparameters — untuned by design.** A single fixed configuration is used for *every* drug and
  *every* organism (`n_estimators=300, max_depth=4, learning_rate=0.1, subsample=0.9,
  colsample_bytree=0.8, scale_pos_weight=<class ratio>, random_state=42`). We deliberately do **not**
  tune per drug: searching hyperparameters inside the same lineage-grouped folds that also select the
  VME-constrained operating point is a realistic route to optimistic bias on the small resistant class.
  The honest description of the method is **"untuned defaults, tuned threshold."**
- **Operating point:** threshold chosen on out-of-fold TRAIN predictions to cap **very-major error at
  ≤3%** (CLSI/FDA-style), evaluated on the held-out TEST lineages.
- **Interpretability:** global SHAP (TreeExplainer) per drug; top determinants checked against known
  mechanisms.

## 4. Results (full 3,850-genome set, VME ≤ 3% operating point)
Per drug, test set = unseen lineages. Logistic regression shown (competitive-or-best on these sparse
features; XGBoost similar — full tables incl. rules + calibrated XGBoost in
`results/reports/summary_05_week2_panel.md`). VME = resistant called susceptible (dangerous);
ME = susceptible called resistant.

| Drug | Test R/S (lineages) | ROC-AUC | PR-AUC | VME | ME | Rules ROC | What ML adds |
|---|---|---|---|---|---|---|---|
| ciprofloxacin | 568/225 (139) | 0.973 | 0.986 | 0.033 | 0.089 | 0.904 | halves major errors |
| TMP-SMX | 512/226 (123) | 0.970 | 0.981 | 0.039 | 0.128 | 0.875 | discrimination + ME |
| meropenem | 355/526 (132) | 0.982 | 0.976 | 0.017 | 0.293 | 0.934 | meets VME target |
| gentamicin | 401/528 (141) | 0.968 | 0.957 | 0.032 | 0.129 | **0.967** | ~matches strong rule |
| **cefoxitin** | 381/253 (77) | **0.906** | 0.934 | 0.018 | 0.850 | **0.517** | **large — porin loss** |

**Honest reading.** Discrimination is strong on unseen lineages (ROC 0.91–0.98). The enlarged data
**stabilised the previously weak operating points** — meropenem and cefoxitin now meet the VME≤3%
target (were data-limited/degenerate before). ML clearly beats the gene-lookup baseline on cefoxitin,
TMP-SMX, ciprofloxacin; on **gentamicin the rules baseline is already excellent** (ROC 0.967 — AME
presence is a direct signal) and ML only matches it — reported honestly.

**The clinical trade-off is explicit:** operating at a strict VME≤3% raises major errors, sharply so
where resistance is hard to call precisely (cefoxitin ME 0.85, meropenem ME 0.29). A screening tool
that refuses to miss resistance necessarily over-calls it; the operating threshold is a clinical
choice (⚕), and ROC/PR-AUC (threshold-free) show the underlying separation is good.

### 4.1 Comparison to published methods
Published genome-based AMR-ML for *K. pneumoniae* reports strong numbers — meropenem ROC-AUC ≈ 0.93,
gentamicin accuracy ≈ 0.91, ciprofloxacin accuracy ≈ 0.96, and large multi-drug studies AUC > 0.9
with ~96% agreement (Nguyen et al. 2018, the *K. pneumoniae* in-silico MIC panel [ref 7]; and
Nguyen-style pipelines more broadly [ref 8]). **But most use random
train/test splits**, which leak population structure and inflate metrics. Our ROC-AUC (0.91–0.98) sits
in or above that published band **while measured under a stricter phylogeny-aware split** — i.e.,
comparable discrimination on a harder, honest evaluation, plus interpretability and clinical-error
reporting that most studies omit. The claim is not "state of the art"; it is "honestly measured, here
is what the model adds over a gene lookup." See `results/figures/benchmark_summary.png` (published band
shown in gold).

### 4.2 Generalization across eight organisms
The identical pipeline runs on **eight WHO-priority pathogens** spanning the Gram divide (organism-agnostic
runner; only taxon / AMRFinderPlus organism / MLST scheme change) — **Enterobacterales:** *K. pneumoniae*,
*E. coli*, *Salmonella enterica*; **Gram-negative non-fermenters (WHO critical):** *A. baumannii*,
*P. aeruginosa*; **Gram-positives:** *S. aureus*, *E. faecium* (VRE), *S. pneumoniae*.
Every organism × drug reaches ROC-AUC **0.84–0.998**
on unseen lineages (`summary_18`, `#26–29`). ML clearly beats the gene-lookup where resistance is
combinatorial or regulatory (S. aureus cefoxitin 0.93 vs 0.80; A. baumannii amikacin 0.90 vs 0.54;
P. aeruginosa ceftazidime 0.87 vs 0.70; E. faecium ampicillin captures *pbp5*; S. pneumoniae penicillin
captures the *pbp1a/2b/2x* mosaic) and matches it on direct single-gene calls (mecA, van, carbapenemases).
The comparison is made fair per organism by **organism-specific rule baselines** (`ORG_RULES`), since
"ampicillin" means β-lactamases in E. coli but *pbp5* in Enterococcus. Honest limits surfaced per
organism: A. baumannii carbapenems degenerate at VME≤3%; Salmonella chloramphenicol (ML < rules);
S. pneumoniae TMP-SMX (causal folA/folP not in the determinant catalog → co-selection-driven).
**E. cloacae and C. jejuni were evaluated and excluded** for insufficient public lab data — honest
data-driven organism selection, not silent omission.

**Cross-species transfer (`cipro_transfer_matrix.png`, Summary #17).** Ciprofloxacin is modelled in
**six** organisms; trained on one and tested **zero-shot** on another (shared determinants only, no
target-organism training), transfer tracks **phylogenetic distance** — which is exactly what a model
that learned mechanism, rather than lineage artifacts, should do:

| Transfer regime | Zero-shot ROC-AUC | Interpretation |
|---|---|---|
| Within **Enterobacterales** (K. pneumoniae ↔ E. coli ↔ Salmonella) | **0.95–0.98** | identical gyrA/parC residue numbering — near-lossless transfer |
| Enterobacterales ↔ **non-fermenters** (A. baumannii, P. aeruginosa) | **0.39–0.89** | homologous targets but divergent sequence context and additional efflux-mediated resistance — degraded and directionally uneven |
| Any ↔ **Gram-positive** *S. aureus* | **0.46–0.69** | uses *grlA* (not *parC*) and different gyrA numbering — collapses, as expected |

The honest reading is that transfer is **excellent within a family, unreliable across families, and
absent across the Gram divide**. The earlier draft quoted a single "0.74–0.98" Gram-negative band; that
was computed when only four organisms were modelled and is superseded by the six-organism matrix above,
which includes genuinely poor pairs (e.g. Salmonella→P. aeruginosa 0.393). Reported as measured.

### 4.3 Statistical rigor — significance, confidence intervals, and quantified leakage
Three additions put the central claim on a defensible statistical footing (`summary_24`, `summary_25`;
`src/evaluation/clinical_rigor.py`, `leakage_delta.py`), *K. pneumoniae*, pooled held-out predictions
from full lineage-grouped CV, in FDA/CLSI vocabulary (bars: **VME ≤ 1.5%** / ≤3% tolerated, **ME ≤ 3%**,
**CA ≥ 90%**, **EA ≥ 90%**):

- **The model significantly beats the rules baseline (DeLong paired test) on every drug.** Pooled
  lineage-held-out ROC-AUC, model vs rules: meropenem 0.968 vs 0.933 (p = 5.3×10⁻²⁴), gentamicin
  0.981 vs 0.948 (p = 1.4×10⁻³⁰), ciprofloxacin 0.983 vs 0.911 (p = 6.9×10⁻³⁶), TMP-SMX 0.977 vs
  0.885 (p = 9.4×10⁻⁵¹), and **cefoxitin 0.905 vs 0.500** (p ≈ 0, below double-precision resolution).
  This is the project's core claim ("what ML adds over a gene lookup") with a significance test behind
  it, not just a gap.
- **Lineage-clustered bootstrap 95% CIs (2,000 resamples)** on VME/ME/CA/AUC — essential because the
  resistant class is the smaller one, so VME point estimates are noisy. **Every drug's VME *point
  estimate* clears the ≤3% bar** (0.66–2.33%), but the honest caveat is the interval: meropenem is
  **1.26% [0.33–3.49%]**, whose upper bound *crosses* 3%, and ciprofloxacin is 2.33% [1.26–5.05%].
  With lineage-clustered resampling the effective sample size is the number of lineages, not genomes,
  so these intervals are wide by construction. **We therefore claim the VME target is met in
  expectation, not guaranteed at the upper confidence bound** — the PASS/FAIL flags in
  `clinical_rigor.json` record both. The high ME on meropenem (55.9%) and cefoxitin (94.2%) is the
  VME-safety trade-off, shown rather than hidden.
- **Population-structure leakage quantified:** the same model scored under random vs lineage-held-out
  vs temporal splits. Random-split AUC is inflated on every drug (mean **+0.010**, up to +0.032 for
  cefoxitin). Notably the inflation is *small* — because determinant features encode **mechanism, not
  ancestry** (a genome resists because it carries blaKPC, not because of its clade), unlike k-mer
  models that lose 0.1–0.2 AUC on a lineage split. The small ΔAUC is itself evidence the model learned
  biology; we still report the honest lineage-held-out column everywhere.

### 4.4 Figure index
All figures live in `results/figures/` and regenerate with `make figures`.

| Figure | Shows | Section |
|---|---|---|
| `generalization_heatmap.png` | unseen-lineage ROC-AUC, **8 organisms × 18 drugs** (36 calibrated models) | §4.2 |
| `cipro_transfer_matrix.png` | 6-organism ciprofloxacin zero-shot cross-species transfer | §4.2 |
| `benchmark_summary.png` | ROC-AUC per drug (ML vs rules vs published band) + VME/ME, *K. pneumoniae* | §4.1 |
| `roc_curves.png` | per-drug ROC on held-out lineages | §4 |
| `calibration.png` | reliability curves after isotonic calibration | §3 |
| `shap_<drug>.png` (×5) | SHAP beeswarm per drug — the determinants behind each call | §5 |
| `conformal.png` | class-conditional conformal coverage / abstention | §5.4 |
| `risk_coverage.png` | risk–coverage curves + AURC for the abstention layer | §5.5 |
| `decision_curve.png` | standardized net benefit vs threshold probability | §5.6 |

## 5. Interpretability & biological validation ⚕

### 5.1 Global SHAP and mechanism recovery
Global SHAP (`results/figures/shap_*.png`) recovers the correct causal mechanisms, and the enlarged
dataset **sharpened the interpretation**: meropenem now ranks **blaKPC-2/blaKPC-3 as the top-2
determinants** (the added lineage diversity broke the earlier gene–lineage confound, where a
co-carried fluoroquinolone mutation topped the list). Top features per drug: parC/gyrA + qnr (cipro),
aac(3) acetyltransferases (gentamicin), sul1/dfrA (TMP-SMX), blaKPC (meropenem), ompK36 porin
(cefoxitin). **Residual co-selection:** cefoxitin still shows co-carried MDR markers high because
porin-loss resistance co-occurs with carbapenemase-carrying lineages — per-instance SHAP (the
beeswarm), not mean rank, remains the honest lens, and SHAP rank ≠ causation.

### 5.2 Cross-resistance / linked-determinant scan
(`src/interpret/novel_determinants.py`, `summary_11`) High-SHAP determinants that are not a drug's causal mechanism are, on inspection, known
AMR genes for *other* classes co-carried on the same MDR plasmids (e.g. meropenem resistance travels
with aac(6')-Ib, blaOXA, blaTEM, ble). This maps the isolate population's linked-resistance structure.
**Honest limit:** because features are AMRFinderPlus *catalogued* determinants, none is a novel gene —
genuine novel-gene discovery would require pan-genome/k-mer features (out of scope, see §7).

### 5.3 MIC regression — predicting the resistance *level* (gold-standard output)
Beyond binary R/S, the model predicts continuous **MIC** (minimum inhibitory concentration — the
clinical gold standard) per drug from the same determinant features, under the phylogeny-aware split
(`src/models/mic_regression.py`, `summary_19`). Evaluated by **Essential Agreement (EA)** — the
CLSI/FDA criterion of prediction within ±1 two-fold dilution: **ciprofloxacin 88.9% EA** (at the ~90%
clinical expectation), cefoxitin 83.6%, gentamicin 81.6%, meropenem 59.4% (harder — wide MIC range and
panel censoring; reported honestly). Pearson r 0.75–0.85. This is a richer, quantitative output scored
by the metric regulators actually use — the model predicts not just *whether* but *how strongly* a
strain resists.

#### 5.3a Protein-language-model extension (ESM-2) — allele resolution for the hard carbapenems
Presence/absence discards *which* variant of a resistance gene a strain carries, yet carbapenem MIC
depends on it (blaKPC-2 vs KPC-3; OXA-23/24/58/51 vs NDM/VIM; porin frameshifts). We extract each
resistance-gene protein per genome (AMRFinderPlus coordinates → translate), embed it with **ESM-2**
(`esm2_t30_150M_UR50D`) so different alleles get different vectors, and add these to the MIC regressor
(`src/models/esm2_mic.py`, `summary_21`). On the hardest drugs this gives a **real, significant gain**,
confirmed by 30 paired lineage-grouped folds:
- *K. pneumoniae* **meropenem**: Essential Agreement **68.4% → 73.2%** (+4.8 pts, Wilcoxon p = 0.0004),
  and fold variance nearly halved.
- *A. baumannii* **imipenem**: **53.8% → 60.5%** (+6.7 pts, p = 0.0001); meropenem +3.6 (p = 0.051,
  borderline, data-limited).
- Neutral on drugs whose determinants already saturate (gentamicin, ciprofloxacin) — exactly as the
  biology predicts. Scaling to ESM-2 **650M** gave no further gain (p = 0.68), so 150M is the honest,
  laptop-scale choice. The whole experiment dedupes 58,502 proteins to 1,945 unique alleles, embedded
  once on the Apple-Silicon GPU (MPS) in minutes.

#### 5.3b Foundation models honestly evaluated (what we tried and rejected)
We also tested an **isolate-similarity GNN** (AMR-GNN style: each genome a node, edges by genomic
similarity — the architecture whose published gains target hard drugs; `summary_20`). On identical
lineage-grouped folds it was **worse than a plain MLP** on the same features for all five drugs (robust
across graph density), because our curated-determinant features already saturate the mechanistic signal
that a population-structure graph would otherwise supply. Reported as an honest negative — "we
implemented the SOTA graph model and it did not win, for a reason we can explain." Net: ESM-2 (allele
resolution) helps where variant identity matters; GNN does not — a principled, tested model-selection
story rather than an assumed one.

### 5.4 Uncertainty & clinical abstention (conformal prediction)
Class-conditional (Mondrian) conformal prediction gives each call an **empirically-validated per-class
coverage level** (the exchangeability guarantee is deliberately broken by the lineage-held-out split,
so we measure coverage rather than assume it) and an explicit **abstain → defer to phenotypic testing**
option (`src/models/conformal.py`,
`results/figures/conformal.png`). At α=0.05 the model confidently classifies ~95% of strains for
well-powered drugs (cipro, TMP-SMX) and, for hard drugs (meropenem, cefoxitin), makes confident calls
on ~55–58% with **near-zero very-major error (0.3–0.6%)** while deferring the rest. Because we test on
*unseen lineages* (breaking the exchangeability conformal assumes), coverage is **empirically
validated** under distribution shift — stricter than standard random-split conformal, and it holds
(≥0.93–1.00) for most drug/class pairs.

### 5.5 Risk–coverage & AURC — evaluating the abstention (Summary #31)
Beyond *having* an abstain option, we *evaluate* it: ranking calls by confidence and deferring the
least-confident strains to the lab gives a **risk–coverage curve** (`results/figures/risk_coverage.png`).
Deferring the least-confident **~30%** roughly halves both error and the clinically-critical **VME**
(meropenem VME 8.7%→3.5%, cefoxitin 19.4%→9.2% at 70% coverage); AURC (area under the error curve,
lower better) is 0.013–0.024 for the well-powered drugs and 0.056 for cefoxitin — the model knows when
it doesn't know.

### 5.6 Decision-curve analysis — clinical net benefit (Summary #30)
The clinical-utility axis most AMR-ML omits: standardized **net benefit** vs the clinician's threshold
probability p_t (`results/figures/decision_curve.png`). The model yields more net benefit than treating
everyone or trusting a gene lookup across the **low-p_t** region where a clinician who fears missing
resistance operates — for meropenem/gentamicin/TMP-SMX the model beats the rules baseline over the
*entire* p_t∈[0.01,0.5] range; for ciprofloxacin only above p_t≈0.33 (where the gyrA/parC lookup is
already strong — reported honestly). This ties the model's value directly to the VME-first cost framing,
not just to AUC.

### 5.7 Imbalance-robust metrics
Because the resistant class is the minority for most drugs, the regulator table (§4.3, Summary #24) also
reports **MCC** (Matthews correlation — robust to imbalance, unlike accuracy/CA) and **PPV/NPV at a
realistic 30% local resistance prevalence** (full 10/30/50% set in `clinical_rigor.json`), since
predictive value depends on prevalence, not just sensitivity/specificity.

## 6. Demo
An **eight-organism interactive web app** (`src/app/streamlit_app.py`; 36 deployable models under
`results/models/<organism>/`). Pick an organism (K. pneumoniae / E. coli / A. baumannii / S. aureus /
P. aeruginosa / Salmonella enterica / E. faecium / S. pneumoniae) — the model set is guarded by an
MLST + intrinsic-marker **species check** that withholds predictions on a wrong-species upload —
choose a cached genome or **upload a FASTA** (`.fna/.fasta/.fa/.gz/.faa`, gzip &amp; molecule auto-detected) to run
AMRFinderPlus live, and get per-drug R/S at the VME≤3% threshold, calibrated P(resistant), and the
SHAP-ranked determinants behind each call. An optional **AI clinical narrative** renders the findings
into clinician-readable prose — an *explanation layer only* (it never changes a call), backed by
Claude, DeepSeek, or Google Gemini (whichever API key is set; deterministic templated fallback with
none). CLI equivalent: `python -m src.app.predict --organism saureus --genome <fasta>`.

## 7. Limitations
- **Scope: 8 organisms, ~3–5 drugs each** — K. pneumoniae is the deep-dive (5 drugs, full MIC +
  conformal + rigor); the other seven are validated at panel scope. Not a universal tool — organisms
  were added only where public lab data supports honest modelling (E. cloacae, C. jejuni excluded).
- **Co-selection confounding** inflates the apparent importance of co-carried determinants (§5);
  reduced but not eliminated (still visible for cefoxitin; S. pneumoniae TMP-SMX is co-selection-driven
  because the causal folA/folP mutations are not in the determinant catalog).
- **High major-error cost at VME≤3%** for hard-to-call drugs (cefoxitin ME 0.85, meropenem 0.29) —
  a clinical operating-point trade-off, not a discrimination failure.
- **Reference-database bias** — AMRFinderPlus/CARD and BV-BRC over-represent well-sequenced regions
  and pathogens; generalization beyond them is unproven.
- **Determinant features only** — k-mer/whole-genome signal deliberately out of scope (an isolate-
  similarity GNN was tested and gave no gain on these curated features, §5.3b).
- **A. baumannii meropenem MIC gain is borderline** (p = 0.051) — honestly under-powered at n=449.
- **Published comparison is band-level, not a re-run head-to-head** (§4.1) — we place our numbers
  within the published range under a stricter split rather than re-executing prior pipelines.
- **No externally-sourced validation cohort (the principal open limitation).** Training *and* testing
  draw on a single aggregator, BV-BRC. The lineage-held-out split makes the evaluation honest *within*
  that corpus, and §4.3 quantifies the population-structure inflation we avoid — but it cannot answer
  whether performance holds on isolates curated by someone else, under different laboratory AST
  practice and a different geographic/temporal sampling frame. A temporal split (`summary_12`) is a
  partial proxy only. Until a genuinely independent cohort is scored, all reported numbers should be
  read as **internal validation under a strict split**, not as external validation.
  *In progress (`summary_32`):* an independent cohort has been **scoped but not yet scored**. Note
  that BV-BRC ingests AMR phenotypes *from* NCBI BioSample/Antibiogram records, so NCBI Pathogen
  Detection is **upstream** of our corpus and a naive comparison would re-test training isolates;
  after a BioSample/assembly anti-join removed the 16.1% that overlapped, **1,038 unseen
  K. pneumoniae isolates** remain (496–780 usable R/S per drug). The Step-2 analysis plan is
  **pre-registered** in `summary_32` — frozen models, no re-tuning, all five drugs reported whatever
  the result. No external number exists in this report yet.
- **Cross-species transfer degrades outside a family** (§4.2) — zero-shot transfer is strong within
  Enterobacterales but unreliable to non-fermenters and absent across the Gram divide; the models are
  organism-specific tools, not a universal predictor.

## 8. Ethics
Research and decision-support only — **not a diagnostic device**, and no substitute for phenotypic
susceptibility testing; a missed resistant call has real clinical cost, hence the emphasis on
very-major error. Only public, de-identified genomic data are used. Interpretability is itself a
safeguard: every prediction is auditable via its determinants.

## 9. Reproducibility
Pinned `environment.yml` (incl. `torch`/`fair-esm`/`biopython`/`scipy` for the ESM-2 pipeline), fixed
seed (42), config-driven runs, `data/manifest.md` (all eight organisms + ESM-2 derived data, sources,
queries, tool versions), and `docs/decisions.md` (54 logged decisions). Tests (`pytest`) cover the
data joins, the feature builder, and — most importantly — the **lineage splitter**, asserting zero
train/test lineage overlap **for every organism** (parametrized over all 8 split CSVs). Makefile
targets regenerate each stage (`make features`, `split`, `train`, `eval`, `organisms ORG=…`,
`leakage`, `models`, `figures`). Raw data regenerates from the manifest;
`data/raw/` is git-ignored.

## References

*Every DOI below was resolved against Crossref, PubMed, doi.org, or the publisher; entries with no
DOI are standards documents that are not DOI-registered, and are marked as such rather than given a
constructed identifier.*

**Data sources and reference databases**

1. Olson RD, Assaf R, Brettin T, Conrad N, Cucinell C, Davis JJ, et al. (2023). Introducing the
   Bacterial and Viral Bioinformatics Resource Center (BV-BRC): a resource combining PATRIC, IRD and
   ViPR. *Nucleic Acids Research* 51(D1):D678–D689. DOI: 10.1093/nar/gkac1003
2. Feldgarden M, Brover V, Gonzalez-Escalona N, Frye JG, Haendiges J, Haft DH, et al. (2021).
   AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among
   antimicrobial resistance, stress response, and virulence. *Scientific Reports* 11(1):12728.
   DOI: 10.1038/s41598-021-91456-0
3. Alcock BP, Huynh W, Chalil R, Smith KW, Raphenya AR, Wlodarski MA, et al. (2023). CARD 2023:
   expanded curation, support for machine learning, and resistome prediction at the Comprehensive
   Antibiotic Resistance Database. *Nucleic Acids Research* 51(D1):D690–D699. DOI: 10.1093/nar/gkac920

**Genome quality control and typing**

4. Bowers RM, Kyrpides NC, Stepanauskas R, Harmon-Smith M, Doud D, Reddy TBK, et al. (2017). Minimum
   information about a single amplified genome (MISAG) and a metagenome-assembled genome (MIMAG) of
   bacteria and archaea. *Nature Biotechnology* 35(8):725–731. DOI: 10.1038/nbt.3893
5. Parks DH, Imelfort M, Skennerton CT, Hugenholtz P, Tyson GW (2015). CheckM: assessing the quality
   of microbial genomes recovered from isolates, single cells, and metagenomes. *Genome Research*
   25(7):1043–1055. DOI: 10.1101/gr.186072.114
6. Lam MMC, Wick RR, Watts SC, Cerdeira LT, Wyres KL, Holt KE (2021). A genomic surveillance framework
   and genotyping tool for *Klebsiella pneumoniae* and its related species complex. *Nature
   Communications* 12(1):4188. DOI: 10.1038/s41467-021-24448-3

**Prior work in genomic AMR prediction**

7. Nguyen M, Brettin T, Long SW, Musser JM, Olsen RJ, Olson R, et al. (2018). Developing an in silico
   minimum inhibitory concentration panel test for *Klebsiella pneumoniae*. *Scientific Reports*
   8(1):421. DOI: 10.1038/s41598-017-18972-w
8. Nguyen M, Long SW, McDermott PF, Olsen RJ, Olson R, Stevens RL, et al. (2019). Using machine
   learning to predict antimicrobial MICs and associated genomic features for nontyphoidal
   *Salmonella*. *Journal of Clinical Microbiology* 57(2):e01260-18. DOI: 10.1128/JCM.01260-18
9. Moradigaravand D, Palm M, Farewell A, Mustonen V, Warringer J, Parts L (2018). Prediction of
   antibiotic resistance in *Escherichia coli* from large-scale pan-genome data. *PLOS Computational
   Biology* 14(12):e1006258. DOI: 10.1371/journal.pcbi.1006258
10. Hicks AL, Wheeler N, Sánchez-Busó L, Rakeman JL, Harris SR, Grad YH (2019). Evaluation of
    parameters affecting performance and reliability of machine learning-based antibiotic
    susceptibility testing from whole genome sequencing data. *PLOS Computational Biology*
    15(9):e1007349. DOI: 10.1371/journal.pcbi.1007349 — *the population-structure leakage reference
    motivating our non-negotiable #1.*

**Methods — models, interpretability, uncertainty, clinical utility**

11. Chen T, Guestrin C (2016). XGBoost: a scalable tree boosting system. In: *Proceedings of the 22nd
    ACM SIGKDD International Conference on Knowledge Discovery and Data Mining (KDD '16)*, pp. 785–794.
    DOI: 10.1145/2939672.2939785
12. Lundberg SM, Lee S-I (2017). A unified approach to interpreting model predictions. In: *Advances
    in Neural Information Processing Systems 30 (NIPS 2017)*. arXiv:1705.07874.
    DOI: 10.48550/arXiv.1705.07874 *(NeurIPS proceedings assign no publisher DOI.)*
13. DeLong ER, DeLong DM, Clarke-Pearson DL (1988). Comparing the areas under two or more correlated
    receiver operating characteristic curves: a nonparametric approach. *Biometrics* 44(3):837–845.
    DOI: 10.2307/2531595
14. Angelopoulos AN, Bates S (2023). Conformal prediction: a gentle introduction. *Foundations and
    Trends in Machine Learning* 16(4):494–591. DOI: 10.1561/2200000101
15. Vickers AJ, Elkin EB (2006). Decision curve analysis: a novel method for evaluating prediction
    models. *Medical Decision Making* 26(6):565–574. DOI: 10.1177/0272989X06295361
16. Lin Z, Akin H, Rao R, Hie B, Zhu Z, Lu W, et al. (2023). Evolutionary-scale prediction of
    atomic-level protein structure with a language model. *Science* 379(6637):1123–1130.
    DOI: 10.1126/science.ade2574 — *ESM-2.*

**Reporting and clinical standards**

17. Collins GS, Moons KGM, Dhiman P, Riley RD, Beam AL, Van Calster B, et al. (2024). TRIPOD+AI
    statement: updated guidance for reporting clinical prediction models that use regression or
    machine learning methods. *BMJ* 385:e078378. DOI: 10.1136/bmj-2023-078378 — *see Appendix A.*
18. Clinical and Laboratory Standards Institute (2026). *Performance Standards for Antimicrobial
    Susceptibility Testing.* 36th ed. CLSI supplement M100. Wayne, PA: CLSI. **No DOI** (CLSI
    standards are not DOI-registered).
19. U.S. Food and Drug Administration (2009). *Class II Special Controls Guidance Document:
    Antimicrobial Susceptibility Test (AST) Systems.* CDRH, guidance no. 631. **No DOI.**
    *Terminology note: the FDA document says "very major discrepancy" and "major discrepancy"; the
    "very-major error / major error" phrasing used throughout this report is the CLSI/field
    convention for the same quantities.*
20. Ellington MJ, Ekelund O, Aarestrup FM, Cantón R, Doumith M, Giske C, et al. (2017). The role of
    whole genome sequencing in antimicrobial susceptibility testing of bacteria: report from the
    EUCAST Subcommittee. *Clinical Microbiology and Infection* 23(1):2–22.
    DOI: 10.1016/j.cmi.2016.11.012

---

## Appendix A — TRIPOD+AI reporting checklist

Mapped against the official **TRIPOD+AI** checklist (Collins et al., *BMJ* 2024;385:e078378 —
27 numbered items, 52 rows including sub-items; checklist version 11-Jan-2024).

**Read this appendix as a disclosure, not a compliance badge.** TRIPOD+AI was written for *clinical*
prediction models built on human participants. This study predicts a laboratory phenotype for a
**bacterial isolate**, so a number of items (informed consent, treatments received, blinding of
outcome assessors, patient and public involvement) have no meaningful analogue and are marked ➖ **N/A**
with the reason given. Several others are honestly marked ❌ or ⚠️ — most importantly the entire
**evaluation (E) track**, because no externally-sourced validation cohort has been scored (§7).

**Key:** ✅ reported · ⚠️ partially reported / caveated · ❌ not done · ➖ not applicable

| Item | Requirement (abbreviated) | Status | Where / note |
|---|---|---|---|
| 1 | Title identifies model, population, outcome | ✅ | Title + subtitle |
| 2 | Abstract | ⚠️ | Abstract present; the separate *TRIPOD+AI for Abstracts* checklist was not formally applied |
| 3a | Healthcare context & rationale; existing models | ✅ | §1; §4.1 (published methods) |
| 3b | Target population, intended purpose & users | ✅ | §1, §6, §8 — research / decision-support, explicitly **not a diagnostic device** |
| 3c | Known health inequalities between sociodemographic groups | ⚠️ | §7 — reframed for this domain: AMRFinderPlus/CARD and BV-BRC **over-represent well-sequenced (high-income) regions**, so performance in under-sequenced settings is unproven. This is the genuine equity axis here; no patient-level sociodemographic data exist |
| 4 | Objectives; development vs validation | ✅ | §1 — **development + internal validation only** |
| 5a | Data sources & representativeness | ✅ | §2, `data/manifest.md` — BV-BRC, laboratory-evidence phenotypes only |
| 5b | Dates of participant data | ⚠️ | Query dates in `data/manifest.md`; isolate collection years used for the temporal split (`summary_12`). Collection dates are missing for a subset of isolates |
| 6a | Study setting, number and location of centres | ⚠️ | Not enumerable — BV-BRC aggregates isolates from many laboratories and countries; per-centre provenance is not consistently available upstream |
| 6b | Eligibility criteria | ✅ | §2 — organism taxon, `evidence = "Laboratory Method"`, MIMAG-aligned QC gate |
| 6c | Treatments received | ➖ | N/A — bacterial isolates, not treated patients |
| 7 | Data pre-processing & quality checking | ✅ | §2 (QC gate), §3 (feature construction) |
| 8a | Outcome definition & rationale | ⚠️ | §2–§3 — binary R/S from laboratory AST. **Caveat:** `Intermediate` isolates are excluded; this is an open decision pending microbiology sign-off (`docs/decisions.md`) |
| 8b | Outcome assessor qualifications | ➖ | N/A — outcome is a standardised laboratory AST result, not a subjective human assessment |
| 8c | Blinding of outcome assessment | ➖ | N/A — phenotypes were generated independently and long before modelling |
| 9a | Choice of initial predictors; pre-selection | ✅ | §3 — the **entire** AMRFinderPlus determinant catalogue for the organism; no pre-selection or supervised filtering |
| 9b | Predictor definition & measurement | ✅ | §3 — AMRFinderPlus, tool and database version pinned (DB 2026-05-15.1) |
| 9c | Predictor assessor subjectivity | ➖ | N/A — fully automated annotation |
| 10 | Study size justification / sample size calculation | ❌ | **No formal sample size calculation was performed.** Achieved sizes and per-drug R/S counts are reported (§2, §4); drugs with too few labels were **dropped rather than modelled** (e.g. P. aeruginosa amikacin, E. faecium cipro/erythromycin) |
| 11 | Missing data handling | ✅ | §2 — per drug, genomes without that drug's phenotype are excluded from that drug's model; QC failures dropped. **No imputation of genomes, labels, or features** |
| 12a | How data were partitioned | ✅ | §3 — lineage-aware `StratifiedGroupKFold`; zero-overlap enforced by `tests/test_split.py` |
| 12b | Predictor handling (form, scaling) | ✅ | §3 — binary presence/absence; no transformation or rescaling |
| 12c | Model type, rationale, model-building, hyperparameter tuning, internal validation | ✅ | §3 — includes the explicit **"untuned defaults, tuned threshold"** disclosure |
| 12d | Heterogeneity across clusters | ✅ | §4.3 — **lineage-clustered** bootstrap (effective n = lineages); §4.2 — per-organism results |
| 12e | Performance measures & rationale | ✅ | §3, §4 — VME/ME first, plus ROC-AUC, PR-AUC, MCC, Brier, Essential Agreement, net benefit, AURC |
| 12f | Model updating from evaluation | ➖ | N/A — no external evaluation performed |
| 12g | How predictions are computed | ✅ | §6; `src/app/predict.py`; 36 serialised `.joblib` models in `results/models/` |
| 13 | Class imbalance methods | ✅ | §3 (`scale_pos_weight`), §5.7 (imbalance-robust metrics). Synthetic oversampling (SMOTE) was **considered and rejected** — interpolating binary determinant vectors invents biologically meaningless genomes (`summary_23`) |
| 14 | Fairness approaches | ⚠️ | No sociodemographic fairness analysis is possible (no patient-level data). The domain analogue — geographic/database representation bias — is disclosed in §7 but **not quantified**; a genuine gap |
| 15 | Model output & threshold derivation | ✅ | §3 — isotonic-calibrated P(resistant); threshold set on **out-of-fold training** predictions to cap VME ≤3%, never on test |
| 16 | Differences between development and evaluation data | ➖ | N/A — no separate evaluation dataset exists. **This is the study's principal reporting gap** (§7) |
| 17 | Ethical approval / informed consent | ➖ | N/A — public, de-identified bacterial genomic data; no human participants. See §8 |
| 18a | Funding | ⚠️ | To be completed at submission |
| 18b | Conflicts of interest | ⚠️ | To be completed at submission |
| 18c | Protocol availability | ✅ | `proposal/`, `PHASES.md`, `CLAUDE.md` — the scope and phase plan were fixed **before** data acquisition |
| 18d | Registration | ➖ | **Not registered.** Stated plainly; no prospective registration was made |
| 18e | Data availability | ✅ | All source data public (BV-BRC); every query, version and date in `data/manifest.md`; `data/raw/` regenerates from the manifest |
| 18f | Code availability | ✅ | Full analytical code public on GitHub; pinned `environment.yml`; `Makefile` regenerates every stage |
| 19 | Patient & public involvement | ➖ | N/A — no patient or public involvement; no human participants |
| 20a | Flow of participants through the study | ⚠️ | Counts reported at each stage (acquired → QC-passed → per-drug labelled) in §2 and the per-organism summaries; **no flow diagram** is drawn |
| 20b | Characteristics of the data | ⚠️ | Per-drug R/S counts and lineage counts reported (§4); richer isolate metadata (source, country, year) is incomplete upstream |
| 20c | Comparison of predictor distributions, development vs evaluation | ➖ | N/A — no external evaluation cohort |
| 21 | Number of participants and events per analysis | ✅ | §4 table; `results/metrics/clinical_rigor.json` reports n and n_R per drug |
| 22 | Full model specification | ✅ | 36 trained models shipped in the repository + complete training code — third parties can reproduce or re-evaluate directly |
| 23a | Performance with confidence intervals | ✅ | §4.3 — lineage-clustered bootstrap 95% CIs on VME/ME/CA/AUC, including the intervals that **cross** the clinical bar |
| 23b | Heterogeneity in performance across clusters | ✅ | §4.2 (8 organisms × 18 drugs), §4.3, `generalization_heatmap.png` |
| 24 | Results of model updating | ➖ | N/A — no external evaluation performed |
| 25 | Overall interpretation | ✅ | §4 "Honest reading", §5, §4.1 |
| 26 | Limitations | ✅ | §7 — including the absence of external validation as the principal limitation |
| 27a | Handling poor-quality or unavailable input data | ✅ | §6 — MLST + intrinsic-marker species check withholds predictions on wrong-species uploads; zero-determinant guard; conformal **abstention** defers low-confidence isolates to phenotypic testing (§5.4–5.5) |
| 27b | User interaction and expertise required | ✅ | §6, `README.md` — genome in, per-drug call + calibrated probability + determinants out; intended for users who can interpret an AST result |
| 27c | Next steps, applicability, generalizability | ✅ | §7 and **Appendix B** |

**Self-assessed tally (52 rows):** ✅ 30 reported · ⚠️ 10 partial · ❌ 1 not done · ➖ 11 not applicable.

The single ❌ (item 10, sample size justification) and the concentration of ➖ in the evaluation track
are both real and both point the same direction: this is a **development-and-internal-validation**
study. Presenting it as anything else would be the overclaim the project exists to avoid.

## Appendix B — Next steps

1. **External validation (highest priority — Step 1 complete).** Score the frozen models, unchanged,
   on a cohort assembled independently of BV-BRC, and report the result whatever it is. This is the
   one claim the study currently cannot make (§7, item 16 above). **Step 1 is done** (`summary_32`):
   1,038 genuinely unseen *K. pneumoniae* isolates identified from NCBI Pathogen Detection after
   anti-joining out the 16.1% that overlapped our training set, with the analysis plan pre-registered.
   **Step 2** — download ~5.1 GB of assemblies, annotate (~14 h), score the frozen models — is the
   next action.
2. **Resolve the `Intermediate` decision** (⚕) and, if `I` is folded into `R`, re-run the full panel.
3. **Quantify representation bias** (item 14) — stratify performance by isolate collection region and
   year where metadata permit, rather than only disclosing the bias narratively.
4. **Flow diagram and richer cohort characteristics** (items 20a/20b).
5. **Extend ESM-2 allele embeddings beyond the carbapenems** to other variant-driven resistance
   (blaCTX-M/SHV/TEM, aminoglycoside AACs, PBP mosaics), including explicit flagging of divergent
   alleles a curated catalogue would miss (`summary_23`, idea #9).
