# Summary #33 — External validation: the result

**Date:** 2026-08-22 · **Status:** complete · **Models:** frozen, scored once
**Executed exactly as pre-registered in `summary_32`, before any genome was downloaded.**

---

## Headline

**Performance degrades substantially on externally-curated isolates.** On the *K. pneumoniae*
five-drug panel, unseen-lineage ROC-AUC falls from **0.905–0.983 internally** to **0.596–0.835** on
1,143 isolates the models had never seen, assembled and phenotyped by other groups.

**The core claim survives; the headline claim does not.** "The model adds something over a known-gene
lookup" holds externally on 4 of 5 drugs in both cohorts, but with much smaller margins. "Cefoxitin is
where machine learning dramatically beats the lookup" **does not survive** — externally the gap is
+0.027, and on lineages absent from training the rules baseline is *better* than the model.

## What was run

| | |
|---|---|
| Cohort | 1,143 of 1,144 assemblies (one, `GCA_031046115.1`, was withdrawn from NCBI) |
| Annotation | AMRFinderPlus 4.2.7 / DB 2026-05-15.1 — **re-run by us**, not reused from either source |
| Typing | `mlst` 2.33.1, klebsiella scheme; 1,043 of 1,143 typeable |
| Models | `results/models/kpneu/*.joblib` as shipped — **no retraining, no re-tuning, no threshold re-selection** |
| Baseline | the same organism-aware rules baseline, corrected AmpC list (decision #57) |

## Results

ROC-AUC with lineage-clustered bootstrap 95% CIs.

| Drug | Internal (pooled CV) | External — all | External — ST-novel |
|---|---|---|---|
| meropenem | 0.968 | **0.795** [0.74–0.84] | 0.685 [0.49–0.88] |
| gentamicin | 0.981 | **0.826** [0.77–0.88] | 0.869 [0.77–0.95] |
| ciprofloxacin | 0.983 | **0.835** [0.79–0.88] | 0.821 [0.70–0.92] |
| TMP-SMX | 0.977 | **0.806** [0.75–0.86] | 0.774 [0.61–0.92] |
| **cefoxitin** | 0.905 | **0.596** [0.53–0.66] | 0.592 [0.40–0.78] |

Clinical error rates (external, all isolates): VME **5.2–15.9%**, ME **26.2–87.9%**. The VME ≤ 3%
operating point **does not hold externally** on any drug — it was selected on training out-of-fold
predictions at a different resistance prevalence (see below), and it does not transfer.

### Does ML still beat the gene lookup?

| Drug | External — all | ST-novel |
|---|---|---|
| meropenem | model 0.795 vs rules 0.768 (**+0.027**) | 0.685 vs 0.624 (+0.060) |
| gentamicin | 0.826 vs 0.795 (**+0.030**) | 0.869 vs 0.728 (+0.141) |
| ciprofloxacin | 0.835 vs 0.770 (**+0.065**) | 0.821 vs 0.758 (+0.063) |
| TMP-SMX | 0.806 vs 0.714 (**+0.092**) | 0.774 vs 0.694 (+0.079) |
| **cefoxitin** | 0.596 vs 0.569 (**+0.027**) | **0.592 vs 0.620 (−0.028 — rules win)** |

## Why — three contributing causes, in order of confidence

**1. It is not a feature-extraction artifact.** Checked first, because a broken pipeline would produce
exactly this pattern. External genomes carry **14.4 determinants on average against 15.1 in
training**, 14.2 of them inside the 688-feature vocabulary, and **not one external genome has zero
known determinants**. Out-of-vocabulary determinants are rare and biologically plausible
(`aph(3')-VIa`, `qnrB12`, `blaOXA-204`). The features are sound; the drop is real.

**2. A large resistance-prevalence shift, partly self-inflicted.** Our training set was deliberately
**enriched for resistance** — the Phase-2b top-up targeted resistant strains for
meropenem/gentamicin/cefoxitin (REPORT §2). The external cohort is not enriched:

| Drug | %R internal | %R external | shift |
|---|---|---|---|
| ciprofloxacin | 71.7% | 47.7% | **−24.0 pts** |
| TMP-SMX | 69.4% | 50.4% | −18.9 pts |
| **cefoxitin** | 60.1% | 32.9% | **−27.2 pts** |
| gentamicin | 43.1% | 46.2% | +3.1 pts |
| meropenem | 40.4% | 39.8% | −0.6 pts |

A calibrated threshold chosen at 60% prevalence over-calls badly at 33%, which is most of the ME
inflation. Note this does **not** explain the AUC drop — AUC is threshold-free — so prevalence is
part of the story, not all of it.

**3. Co-selection does not transfer, and cefoxitin is the proof.** Cefoxitin was the drug most
dependent on population-specific structure: REPORT §5.1 already flagged that its high-SHAP features
include co-carried MDR markers, because porin-loss resistance co-occurs with carbapenemase-carrying
lineages *in our corpus*. In a differently-sampled population that co-occurrence breaks, and the
model loses most of its signal (0.905 → 0.596) while the mechanistic AmpC lookup degrades far less
(0.518 → 0.569, i.e. it actually *improves* slightly). **This is the cleanest evidence in the whole
project that co-selection was inflating an internal number** — and it was predicted, in writing, as a
limitation before this test was run.

## What this changes in the write-up

1. **The cefoxitin headline must be scoped to internal validation.** It remains a true and
   well-evidenced *internal* result, and the mechanism (porin loss invisible to a gene lookup) is real
   biology. But it is no longer safe to present it as the project's transferable headline.
2. **The honest headline becomes the generalization gap itself** — measured, not asserted. Very few
   comparable studies report one at all.
3. **The VME ≤ 3% operating point is internal-only** and must be labelled as such everywhere.
4. **The ML-over-lookup claim stands**, externally, on 4 of 5 drugs, at +0.03 to +0.09.

## Honest caveats (all pre-registered in `summary_32`)

- **Independent-source, not independent-population.** Submitter ecosystems overlap.
- **NCBI does not vet submitted AST methods**, so external label noise is plausibly higher than in our
  curated training labels. This biases the external result **downward** by an unknown amount, and we
  cannot separate that from genuine generalization failure.
- **129 training genomes (3.4%) carry no accession** and could not be anti-joined; undetected overlap
  is bounded above by ≤11.3% of the cohort — which would, if anything, make the external number
  *optimistic*.
- **ST-novel intervals are wide** (n = 85–155 per drug). Read them as directional, not precise.
- 100 of 1,143 isolates were untypeable and are therefore absent from the ST-novel cohort.

## Verdict

The models are **not** ready to be pointed at isolates from an arbitrary new source, and this report
now says so with numbers. That is a more useful contribution than another internal AUC: it quantifies
how far a determinant-based AMR classifier travels beyond the corpus it was fitted on, on a
lineage-disjoint external cohort, with the analysis plan fixed in advance.
