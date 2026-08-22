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

## Demo — eight organisms
```bash
# CLI (--organism: kpneu | ecoli | abaumannii | saureus | paeruginosa | senterica | efaecium | spneumoniae)
python -m src.app.predict --organism saureus --genome-id 1280.10000        # instant (cached)
python -m src.app.predict --organism kpneu   --genome path/to/genome.fna   # live AMRFinderPlus (~1-3 min)
python -m src.app.report  --organism abaumannii --genome-id 470.11118      # + AI clinical narrative

# Interactive web app — use `python -m streamlit` from the project env
# (a bare `streamlit` may resolve to a different/broken install on your PATH)
python -m streamlit run src/app/streamlit_app.py       # pick organism -> calls + confidence + drivers + narrative
```
Outputs per-drug resistant/susceptible + calibrated P(resistant) + the determinants behind each call.

**Genome upload** (web app): assembled nucleotide FASTA — `.fna`, `.fasta`, `.fa`, or gzipped `.gz` —
**or a protein FASTA `.faa`** (predicted proteins). Molecule type is auto-detected by content (nucleotide
→ AMRFinderPlus `-n`, protein → `-p`), and gzip is detected by magic bytes, so a mis-named file still
works. One or many contigs; typical size 2–6 MB. A non-FASTA upload gets a clear error. AMRFinderPlus
runs locally — nothing leaves the machine.

**AI narrative (optional).** The clinical narrative is an *explanation layer only* — it never changes
the model's calls — so the provider is interchangeable. Set **one** key (preference order
Claude → DeepSeek → Gemini; force one with `AMR_LLM_PROVIDER`). **Easiest: use a `.env` file** — the
app loads it automatically from the project root, so you don't have to `export` in every terminal:
```bash
cp .env.example .env        # then edit .env and paste ONE key, e.g. DEEPSEEK_API_KEY=sk-...
```
`.env` is git-ignored. A real shell `export` still takes precedence over `.env` if you prefer that.
Keys: `ANTHROPIC_API_KEY` (claude-opus-4-8) · `DEEPSEEK_API_KEY` (deepseek-chat) · `GEMINI_API_KEY`
(gemini-flash-latest, free tier at aistudio.google.com/apikey; override models with `DEEPSEEK_MODEL` /
`GEMINI_MODEL`). Gemini uses Google's native API; DeepSeek uses its OpenAI-compatible endpoint.
With **no key set the demo still works fully** — it prints a deterministic templated report instead.

## Status — externally validated, and the result is honest

**External validation is complete and the models degraded.** On **1,143 *K. pneumoniae* isolates
curated by other groups** — frozen models, no retraining, analysis plan registered *before* any genome
was downloaded — ROC-AUC falls from 0.905–0.983 internally to **0.795–0.835**, and cefoxitin to
**0.596**. The model still beats the known-gene lookup on **4 of 5** drugs (+0.027 to +0.092), but the
VME ≤ 3% operating point does not transfer. A feature-extraction artifact was ruled out first
(14.4 determinants/genome external vs 15.1 training). See `results/reports/summary_33_external_validation.md`
and REPORT §4.4. **Every internal number below should be read as internal validation under a strict split.**

## Internal results — 8 WHO-priority pathogens
- **Delivered:** end-to-end pipeline across **8 organisms** (K. pneumoniae, E. coli, A. baumannii,
  S. aureus, P. aeruginosa, Salmonella enterica, E. faecium, S. pneumoniae), per-drug calibrated
  models (VME≤3% operating point), SHAP + biological validation, interactive demo, ESM-2 allele MIC
  extension, honest GNN rejection. See `results/reports/REPORT.md` and `summary_*.md`.
- Discrimination ROC-AUC **0.84–0.998** on **unseen lineages** (single held-out split; **0.751–0.998** under pooled lineage-grouped CV); ML beats the gene-lookup where
  resistance is combinatorial/regulatory (cefoxitin porin loss, P. aeruginosa ceftazidime, E. faecium
  pbp5, S. pneumoniae pbp mosaics) and matches it on direct single-gene calls (mecA, van, carbapenemases).
- Regulator-grade evaluation (VME/ME/CA + Brier calibration + lineage-clustered CIs + DeLong vs an
  **organism-aware** rules baseline) for all 8 in `results/reports/summary_24_clinical_rigor.md`.
- E. cloacae and C. jejuni evaluated and **excluded** for insufficient public lab data (honest
  data-driven selection). Full audit trail in `docs/decisions.md`.
