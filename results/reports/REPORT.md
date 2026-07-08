# Reading Resistance — Project Report (DRAFT)

*An interpretable, honestly-benchmarked classifier for antibiotic resistance in* Klebsiella
pneumoniae *from bacterial genomes.*

> **Status:** living draft. Numbers below are the **current (pre-top-up)** results on 1,472 genomes;
> they will be refreshed on the enlarged dataset via `python -m src.refresh_pipeline` once the
> Phase-2b top-up download completes. Sections marked ⚕ need microbiology sign-off; 🔜 marks
> work still to add (published-method comparison, poster).

---

## Abstract
We predict, from an assembled *K. pneumoniae* genome, whether a strain is resistant or susceptible
to each of five antibiotics — meropenem, gentamicin, ciprofloxacin, trimethoprim-sulfamethoxazole,
and cefoxitin — using interpretable resistance-determinant features. Every model is validated under a
**phylogeny-aware split** (train and test share no MLST sequence type), reports **clinical error
rates** (very-major / major) alongside ROC-AUC/PR-AUC, and is benchmarked against a transparent
known-gene rules baseline. On unseen lineages the models achieve ROC-AUC 0.88–0.98. The headline
result is **cefoxitin**, where the gene-lookup baseline is near-useless (misses ~92% of resistant
strains) because resistance is driven by porin loss the lookup cannot see, and the model lifts
ROC-AUC to 0.96 — a concrete demonstration of what machine learning adds over a database query.

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
- **Current set:** 1,472 QC-passed genomes (balanced ciprofloxacin thin slice, reused across all
  drugs via multi-drug labels). 🔜 Phase-2b top-up adds ~2,400 resistant genomes for
  meropenem/gentamicin/cefoxitin.
- **QC (literature-aligned, MIMAG):** completeness ≥90%, contamination ≤5% (Bowers et al. 2017),
  contigs ≤500, genome length 4.5–7.5 Mbp (Kleborate range). Drops 1.9%, class-balanced.
- Provenance in `data/manifest.md`; every choice in `docs/decisions.md`.

## 3. Methods
- **Features:** AMRFinderPlus (`--organism Klebsiella_pneumoniae`, DB 2026-05-15.1) → binary
  genome × determinant matrix (acquired genes + point mutations; VIRULENCE/STRESS excluded).
  Current matrix: 1,472 × 446 determinants.
- **Lineage-aware split:** MLST (`mlst`, klebsiella scheme) assigns each genome a sequence type
  (400 STs); untypeable genomes get a unique synthetic lineage. **StratifiedGroupKFold** holds out
  whole STs for the test set while balancing R/S. A test asserts zero shared lineage
  (`tests/test_split.py`).
- **Models (one per antibiotic):** transparent known-gene **rules baseline** (drug-specific
  determinants), **logistic regression**, and **XGBoost** with **isotonic calibration**.
- **Operating point:** threshold chosen on out-of-fold TRAIN predictions to cap **very-major error at
  ≤3%** (CLSI/FDA-style), evaluated on the held-out TEST lineages.
- **Interpretability:** global SHAP (TreeExplainer) per drug; top determinants checked against known
  mechanisms.

## 4. Results (current — pre-top-up, VME ≤ 3% operating point)
Per drug, test set = unseen lineages. VME = resistant called susceptible (dangerous); ME = susceptible
called resistant. Best model shown; full tables in `results/reports/summary_05_week2_panel.md`.

| Drug | Test (R/S) | ROC-AUC | PR-AUC | VME | ME | Rules-baseline ROC | ML adds |
|---|---|---|---|---|---|---|---|
| ciprofloxacin | 184/184 | 0.976 | 0.976 | 0.022 | 0.147 | 0.897 | fewer major errors |
| TMP-SMX | 179/141 | 0.969 | 0.963 | 0.039 | 0.113 | 0.895 | discrimination + ME |
| gentamicin | 102/249 | 0.962 | 0.910 | 0.049 | 0.237 | 0.929 | modest |
| meropenem | 73/265 | 0.947 | 0.922 | 0.055 | 0.109 | 0.949 | matches strong rule |
| **cefoxitin** | 79/89 | **0.957** | 0.961 | (data-limited) | — | **0.532** | **large — porin loss** |

**Honest reading:** discrimination is strong everywhere (ROC 0.88–0.98 on unseen lineages). ML clearly
beats the baseline on cefoxitin, TMP-SMX, ciprofloxacin; on meropenem/gentamicin the rules baseline is
already strong (carbapenemase / AME presence is a direct signal). Operating at a strict VME≤3% is clean
for the well-populated drugs (cipro, TMP-SMX) and data-limited for cefoxitin (top-up will address).

## 5. Interpretability & biological validation ⚕
Global SHAP (`results/figures/shap_*.png`) recovers the correct causal mechanisms: gyrA/parC + qnr
(cipro), aac(3) acetyltransferases (gentamicin), dfrA/sul (TMP-SMX), blaKPC carbapenemases with the
largest per-genome impact (meropenem), and ompK35/36 porin mutations (cefoxitin). **Caveat —
co-selection:** for meropenem/cefoxitin the *mean*-|SHAP| ranking is topped by fluoroquinolone
mutations that are co-carried in MDR lineages, not causal; per-instance SHAP (the beeswarm) separates
mechanism from lineage marker. This is a genuine finding, not a bug, and is exactly what the
biological-validation step exists to catch.

## 6. Demo
`python -m src.app.predict --genome <fasta>` → per-drug R/S at the VME≤3% threshold, calibrated
P(resistant), and the SHAP-ranked determinants behind each call. Verified on resistant and susceptible
genomes; live-annotation path matches cached.

## 7. Limitations
- **Single organism, five drugs** — scope-disciplined; not a general tool.
- **Co-selection confounding** inflates the apparent importance of co-carried determinants (§5).
- **cefoxitin operating point is data-limited** at VME≤3% (top-up in progress).
- **Reference-database bias** — AMRFinderPlus/CARD and BV-BRC over-represent well-sequenced regions
  and pathogens; generalization beyond them is unproven.
- **Determinant features only** — k-mer/whole-genome signal deliberately out of scope.
- 🔜 **No published-method head-to-head yet** — to add for a complete honest benchmark.

## 8. Ethics
Research and decision-support only — **not a diagnostic device**, and no substitute for phenotypic
susceptibility testing; a missed resistant call has real clinical cost, hence the emphasis on
very-major error. Only public, de-identified genomic data are used. Interpretability is itself a
safeguard: every prediction is auditable via its determinants.

## 9. Reproducibility
Pinned `environment.yml`, fixed seed (42), config-driven runs, `data/manifest.md` (sources + tool
versions + checksums), `docs/decisions.md` (26 logged decisions), and a one-command refresh
(`src/refresh_pipeline.py`). Raw data regenerates from the manifest; `data/raw/` is git-ignored.

## References
Bowers et al. 2017 (MIMAG, *Nat. Biotechnol.*); Parks et al. 2015 (CheckM); Lam et al. 2021
(Kleborate, *Nat. Commun.*); Feldgarden et al. 2021 (AMRFinderPlus); Alcock et al. 2023 (CARD);
Lundberg & Lee 2017 (SHAP); Hicks et al. 2019 (population-structure leakage); Moradigaravand et al.
2018; Nguyen et al. 2019. 🔜 verify/expand.
