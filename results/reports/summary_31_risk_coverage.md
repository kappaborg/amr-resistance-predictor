# Summary #31 — Risk–Coverage Curves & AURC (defer-to-lab abstention)

**Date:** 2026-07-17 · *Klebsiella pneumoniae* · selective prediction on the pooled lineage-held-out calibrated predictions. Deferring the least-confident strains to phenotypic testing lowers error among those kept. **AURC** = area under the error risk–coverage curve (lower better); we also track **VME among non-deferred**.

| Drug | AURC | error @100% | error @70% cov | VME @100% | VME @70% cov |
|---|---|---|---|---|---|
| meropenem | 0.024 | 5.6% | 2.4% | 8.7% | 3.5% |
| gentamicin | 0.017 | 3.8% | 1.7% | 4.9% | 1.2% |
| ciprofloxacin | 0.013 | 3.6% | 0.8% | 1.8% | 0.8% |
| trimethoprim_sulfamethoxazole | 0.013 | 6.2% | 1.4% | 3.3% | 0.7% |
| cefoxitin | 0.056 | 15.0% | 8.6% | 19.4% | 9.2% |

**Reading.** Both curves fall as coverage drops: deferring the least-confident ~30% of strains to phenotypic testing cuts the error (and the clinically-critical **VME**) among the calls the tool *does* make. This is the quantitative case for the defer-to-lab abstention — the model knows when it doesn't know. Complements the conformal per-class coverage evaluation (Summary #9), which is empirically validated under lineage shift.
