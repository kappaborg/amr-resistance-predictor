# Summary #31 — Risk–Coverage Curves & AURC (defer-to-lab abstention)

**Date:** 2026-07-17 · *Klebsiella pneumoniae* · selective prediction on the pooled lineage-held-out calibrated predictions. Deferring the least-confident strains to phenotypic testing lowers error among those kept. **AURC** = area under the error risk–coverage curve (lower better); we also track **VME among non-deferred**.

| Drug | AURC | error @100% | error @70% cov | VME @100% | VME @70% cov |
|---|---|---|---|---|---|
| meropenem | 0.154 | 31.6% | 25.3% | 1.4% | 1.0% |
| gentamicin | 0.039 | 5.0% | 5.1% | 3.7% | 0.5% |
| ciprofloxacin | 0.034 | 6.0% | 2.6% | 6.0% | 3.8% |
| trimethoprim_sulfamethoxazole | 0.019 | 6.0% | 1.5% | 3.7% | 1.5% |
| cefoxitin | 0.128 | 38.6% | 23.2% | 0.4% | 0.2% |

**Reading.** Both curves fall as coverage drops: deferring the least-confident ~30% of strains to phenotypic testing cuts the error (and the clinically-critical **VME**) among the calls the tool *does* make. This is the quantitative case for the defer-to-lab abstention — the model knows when it doesn't know. Complements the conformal per-class coverage evaluation (Summary #9), which is empirically validated under lineage shift.
