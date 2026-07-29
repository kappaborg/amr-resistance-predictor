# Summary #23 — Deep Research: Improvement Ideas to Boost the Project

**Date:** 2026-07-10 · Two parallel literature sweeps (2022–2026) — (A) SOTA genomic-AMR ML methods,
(B) evaluation / clinical-translation rigor. Below, only ideas that fit our scope (curated
determinants + ESM-2, laptop, honest science) — each tagged **ADOPT / EXPERIMENT / REJECT** with the
evidence and the leakage/honesty check.

---

## A. Cheap, high-value additions we can make from our *existing* outputs (ADOPT)

These need little or no new modelling — they reframe/strengthen what we already produce and are exactly
what reviewers and competition judges reward.

1. **Regulator-grade error thresholds in every metrics table.** State the FDA/CLSI pass bars next to
   our numbers: **VME ≤ 1.5%** (≤3% tolerated with CI caveats), **ME ≤ 3%**, **Minor ≤ 10%**,
   **Categorical Agreement ≥ 90%**, **Essential Agreement ≥ 90%** (our MIC metric *is* EA). Add a
   PASS/FAIL flag per drug. *(FDA AST guidance; CLSI M100/M52.)* — near-zero code, big credibility.
2. **Quantify the leakage we already avoid.** Report random-split vs lineage-held-out vs temporal
   side-by-side and headline the **ΔAUC inflation**. This turns our CLAUDE.md non-negotiable #1 into a
   *result*. *(PLOS Biology 2024, "Biased sampling… confounds ML prediction of AMR".)* We have both
   splits — just tabulate the delta.
3. **Bootstrap 95% CIs + DeLong tests.** Put CIs on VME/ME/AUC/EA (stratified paired bootstrap, ≥2000
   resamples) and a DeLong p-value for model-vs-rules-baseline AUC per drug. Converts "beats baseline"
   into a defensible claim. FDA wants VME with 95% CI anyway.
4. **MCC + prevalence-adjusted PPV/NPV.** One-line metrics that preempt the "accuracy misleads on
   imbalance" critique; recompute PPV/NPV at 2–3 realistic local resistance prevalences.
5. **TRIPOD+AI reporting checklist (appendix).** Map the write-up to the 27-item standard for
   prediction-model reporting *(Collins et al., BMJ 2024)*; cite EUCAST WGS-AST framing and the
   ~91–95% genotype–phenotype concordance benchmark for the external-comparison table.

## B. Distinctive additions worth a targeted experiment (EXPERIMENT)

6. **Decision-Curve Analysis (net benefit).** Plot standardized net benefit vs threshold probability
   for model / rules-baseline / treat-all / treat-none per drug. This is the *clinical-utility*
   dimension most AMR-ML submissions lack, and it ties directly to our VME-heavy cost framing.
   *(Vickers & Elkin 2006; BMJ 2016 reporting guidance.)* Medium effort, high distinctiveness.
7. **Risk–coverage curves + AURC for our conformal layer.** We already do class-conditional (Mondrian)
   conformal; add the standard selective-prediction evaluation — empirical coverage vs nominal per
   class, mean set size / singleton rate, and **VME among non-deferred vs % deferred**. Turns the
   abstention feature into a rigorously-evaluated result. *(Angelopoulos & Bates 2023; AMR precedent
   reached 93.75% coverage flagging 12.5% of cases.)*
8. **Multi-task / multi-drug joint model for low-data drugs.** A shared-trunk MLP on our
   determinant+ESM-2 features with per-drug heads (or a cross-drug meta-learner stacking per-drug
   XGBoost probabilities). Evidence that hard parameter sharing helps rare-resistance drugs
   *(MDPI Antibiotics 2024 PMC11504406; Sci Rep 2026 s41598-026-41185-z)*. **Leakage guard: one global
   lineage-grouped split across all drugs** (a lineage held out for drug A must be held out for every
   head) — add a test asserting group consistency. Treat as an ablation; the strongest evidence is
   from EHR/antibiogram data, so the transfer to determinant features is plausible-but-unproven.
9. **Extend the ESM-2 divergent-allele win beyond carbapenems.** Our biggest differentiator. Apply the
   allele-embedding features to every drug whose resistance is *variant-driven* (blaCTX-M/SHV/TEM,
   aminoglycoside AACs, PBP mosaics) and explicitly use ESM-2 to flag **novel/divergent alleles at
   <30% homology** that a curated catalog misses — the exact false-negative gap of our rules baseline.
   *(Frontiers Microbiology 2025, fmicb.2025.1628952.)*
10. **Multi-scale determinant encoding (presence + copy-number + domain roll-up).** AMRFinderPlus
    already gives coordinates/counts; add count columns and a protein-domain roll-up (interpretable,
    catches gene-dosage effects like duplicated blaOXA / efflux overexpression relevant to MIC).
    *(bioRxiv 2025 PMC12633328, median nMCC 0.89 across ESKAPE.)*

## C. Considered and deliberately rejected (REJECT — good defense material)

11. **Pan-genome / unitig / DBGWAS k-mer features.** Explicitly out of our sprint scope *and* honest-
    science risk: reported models drop from AUC 0.93 (train) to ~0.77 (independent validation) — the
    population-structure overfitting we are specifically defending against — and unitigs are far less
    interpretable. The one real benefit (finding determinants outside CARD) we already capture via the
    ESM-2 divergent-allele route (#9). *(CSBJ 2024 PMC11067008.)*
12. **Whole-genome self-supervised "genome language models" / semi-supervised on unlabelled genomes.**
    On *curated determinants* these mostly re-learn lineage — i.e. reintroduce the leakage we exclude —
    and pretraining at scale doesn't fit a laptop/4-week sprint. The affordable, on-target slice
    (protein-level pretraining) we already exploit via ESM-2. *(bioRxiv 2025.04.01.646674.)*
13. **SMOTE / synthetic oversampling.** Interpolating between binary determinant vectors invents
    biologically meaningless genomes and risks cross-lineage leakage — contradicts "never synthesize."
    Use class weighting / `scale_pos_weight` + threshold tuning instead.
14. **Attention/transformer over determinants as a *replacement*.** On curated presence/absence,
    gradient-boosted trees + SHAP already match or beat attention and stay interpretable — consistent
    with our own finding that a GNN gave no gain here (Summary #20). Only worth it as a narrative device
    for per-determinant attention weights, not for accuracy.

---

## Recommended sequencing (impact-per-effort)
1. **Now, near-zero code:** #1 (regulatory thresholds), #2 (leakage delta), #3 (bootstrap/DeLong CIs),
   #4 (MCC/PPV-NPV) — fold straight into the results tables and poster.
2. **Two distinctive experiments:** #6 (decision-curve analysis) and #7 (risk-coverage/AURC) — the two
   dimensions that separate a strong submission from a typical one (clinical utility + rigorous
   abstention).
3. **If time permits, modelling ablations:** #8 (multi-task, honestly as an ablation) and #9 (extend
   ESM-2 alleles) — both build on what we have and help exactly where the honest baseline fails.
4. **Defense section:** write up #11–14 as "considered and excluded, because they reintroduce
   population-structure leakage or break interpretability."

**Key numbers for the defense:** VME ≤ 1.5%, ME ≤ 3%, Minor ≤ 10%, CA ≥ 90%, EA ≥ 90% (FDA/CLSI);
conformal must hit nominal coverage 1−α *per class*; published genotype–phenotype concordance ≈ 91–95%.

*Full source lists retained in the decision log (#31).*
