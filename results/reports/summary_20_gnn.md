# Summary #20 — Isolate-Similarity GNN (AMR-GNN-style)

**Date:** 2026-07-09 · device `mps` · 5-fold lineage-grouped CV · kNN graph k=10 on cosine feature similarity.

Each **isolate is a node**; edges join genomically similar isolates. A 2-layer GCN does node classification (R/S), modelling population structure directly — the AMR-GNN idea (Nat. Commun. 2026), whose largest gains were on the hardest drugs. Compared, on identical phylogeny-aware folds, against XGBoost (our baseline) and an MLP on the same features with **no graph** (the control).

| Drug | n | %R | XGBoost | MLP (no graph) | **GCN (+graph)** | GCN − MLP |
|---|---|---|---|---|---|---|
| meropenem | 3529 | 40% | 0.968 | 0.967 | **0.933** | -0.034 |
| gentamicin | 3713 | 43% | 0.981 | 0.977 | **0.950** | -0.027 |
| ciprofloxacin | 3176 | 72% | 0.985 | 0.982 | **0.975** | -0.007 |
| trimethoprim_sulfamethoxazole | 2966 | 69% | 0.975 | 0.974 | **0.937** | -0.036 |
| cefoxitin | 2537 | 60% | 0.910 | 0.909 | **0.865** | -0.044 |

**Reading (mean unseen-lineage ROC-AUC).** *GCN − MLP* isolates the value of the graph: both use the identical determinant features and folds, so a positive gap means the isolate-similarity structure adds signal the features alone don't carry. *GCN vs XGBoost* asks whether the graph model beats our strong tree baseline.

**Honesty.** The GCN is **transductive** — it sees all nodes' *features* and the graph during training, but only training-lineage *labels* enter the loss; XGBoost and the MLP are inductive. Lineage grouping puts every member of a lineage in one fold, so no trivial same-lineage label leaks across the graph. This is the same evaluation contract as the rest of the project (unseen lineages), applied fairly to all three models.

## Verdict: the graph does not help on our features (and we can say why)

On every drug the GCN is **worse** than both the plain MLP and XGBoost. The gap is not a badly-tuned
graph: sweeping the neighbour count on meropenem holds it flat and negative throughout —

| kNN k | 3 | 5 | 15 | 30 | 50 |
|---|---|---|---|---|---|
| GCN − MLP (meropenem AUC) | −0.029 | −0.029 | −0.023 | −0.024 | −0.029 |

**Why our result differs from the AMR-GNN paper — the mechanistic reason.** AMR-GNN's large gains
(cefepime +28.8%, aztreonam +18.9%) came on *Pseudomonas aeruginosa* using **unitig** node features —
raw assembly k-mers, which are noisy and only *indirectly* tied to mechanism. There, smoothing labels
across genomically-similar isolates denoises a weak signal, so the graph pays off. **Our node features
are curated AMRFinderPlus determinants** — the carbapenemase, *mecA*, *gyrA*/*parC* calls whose
gene→phenotype map is already near-deterministic. XGBoost/MLP already reach 0.91–0.99 from the features
alone; there is little residual population-structure signal for a graph to add. Worse, kNN smoothing
*blurs* the sharp determinant→phenotype boundary across lineages that carry the same gene at different
prevalences, which is exactly why the GCN loses a couple of points.

**Bottom line.** Testing a GNN here was the right call — it is the on-thesis architecture (it models
clonal structure, our core confounder) and its literature wins are on hard drugs like ours. But the
honest empirical answer on *curated-determinant* features is **no**: the mechanistic signal is already
saturated, so the graph cannot beat gradient boosting. This is a genuine, defensible negative result —
"we implemented the SOTA graph model and it did not win, for a reason we can explain" — which is worth
more to the project than an unverified claim that it would.
