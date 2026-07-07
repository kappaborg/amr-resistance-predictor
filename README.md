# Reading Resistance

An interpretable machine-learning classifier that predicts, from a bacterial genome, whether a
strain is **resistant or susceptible** to each of several antibiotics — and names the genes and
mutations behind every call. Its defining commitment is an **honest benchmark**: validated under a
**phylogeny-aware split** (generalize to unseen lineages, not memorize clones) and compared
head-to-head against a transparent known-gene baseline and published methods.

Fully computational · public data only · laptop / free Colab · four-week sprint · two-person team.

## Non-negotiables
1. Split by **lineage (MLST/cluster), never randomly** — enforced by `tests/test_split.py`.
2. Known-gene rules classifier is the **honest baseline**; report what the model *adds*.
3. **One model per antibiotic.**
4. Report **very-major & major error rates** first, alongside ROC-AUC / PR-AUC.
5. **Determinant features only** (known genes + point mutations); k-mers out of scope.
6. Always **calibrate** probabilities.
7. **Never fabricate** — report weak/negative/missing results plainly.

## Layout
```
CLAUDE.md            engineering brief (source of process truth)
PHASES.md            phase plan (13 phases, summary after each)
proposal/            proposal docx + 3 diagrams (source of scope truth)
config/config.yaml   central run configuration
data/{raw,interim,processed}/   raw is git-ignored; manifest.md tracks sources
src/{data,features,split,models,evaluation,interpret,app}/
tests/               data joins, feature builder, lineage-leakage test
results/{figures,metrics,models,reports}/   reports/ holds per-phase summaries
docs/decisions.md    decision log / lab notebook
environment.yml · Makefile
```

## Quick start
```bash
conda env create -f environment.yml
conda activate amr-resistance-predictor
make help          # list pipeline targets
make test          # run the test suite (incl. the leakage test)
```
Pipeline: `make data → features → split → train → eval → figures` (or `make all`).

## Status
**Phase 0 complete** — scaffold, environment, config, decision log in place.
Next: **Phase 1** — BV-BRC organism/drug selection (needs data plan + microbiology sign-off; no
downloads until approved).
