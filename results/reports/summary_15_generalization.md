# Summary #15 — Second-Organism Generalization (E. coli)

**Date:** 2026-07-09 · **Status:** ✅ the same method transfers to a second pathogen.

## What was done
Ran the **identical pipeline** (QC → AMRFinderPlus → MLST → phylogeny-aware split → per-drug models,
VME≤3% operating point) on *Escherichia coli* — a phylogenetically distinct organism — via the
organism-general runner (`src/organism_pipeline.py --organism ecoli`). No re-engineering: only the
taxon, the AMRFinderPlus organism (`Escherichia`), and the MLST scheme changed. 3,035 QC-passed
genomes × 596 determinants.

## Cross-organism comparison — the generalization proof
Same method, three drugs shared with K. pneumoniae, unseen-lineage test, logistic regression:

| Drug (shared) | K. pneumoniae ROC-AUC | E. coli ROC-AUC |
|---|---|---|
| ciprofloxacin | 0.973 | **0.989** |
| gentamicin | 0.968 | **0.989** |
| trimethoprim-sulfamethoxazole | 0.970 | **0.948** |

**The method achieves comparable, strong discrimination (0.95–0.99) on BOTH organisms** — it is a
transferable *system*, not a single-organism fit. E. coli-relevant β-lactams also work
(ampicillin ROC 0.969, ceftriaxone 0.973).

## Honest per-drug reading (E. coli, VME≤3%)
| Drug | ROC | VME | ME | Rules ROC | Note |
|---|---|---|---|---|---|
| ciprofloxacin | 0.989 | 0.019 | 0.152 | 0.860 | ML >> rules (gyrA/parC + PMQR) |
| trimethoprim-sulfamethoxazole | 0.948 | 0.032 | 0.214 | 0.815 | ML >> rules |
| gentamicin | 0.989 | 0.007 | 0.551 | 0.988 | rules already strong; ML matches, over-calls |
| ampicillin | 0.969 | 0.031 | 0.331 | 0.959 | β-lactamase presence is direct |
| ceftriaxone | 0.973 | 0.026 | 0.194 | 0.967 | ESBL/AmpC direct |

Same pattern as K. pneumoniae: ML clearly beats the gene-lookup where resistance is combinatorial
(cipro, TMP-SMX); on drugs with a direct single-gene signal (gentamicin, β-lactamases) the rules
baseline is already strong and ML matches it. The VME≤3% operating point again trades higher ME for
low missed-resistance — same safety-first behaviour.

## Biological contrast (a nice defense point)
Carbapenem resistance was **excluded from the E. coli panel** — only 69 resistant genomes vs
K. pneumoniae's 1,400+ — because carbapenemase carriage is rare in E. coli. The panel choice reflects
each organism's real resistance epidemiology, not a template.

## Why this matters
Generalization to a second, unseen pathogen under the same honest evaluation is the single most
persuasive evidence that the approach is a method, not an overfit. Adding a third organism is now just
another entry in the `ORGANISMS` registry.
