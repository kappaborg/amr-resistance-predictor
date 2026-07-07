# Phase 1 — BV-BRC Organism & Drug Selection: Query Plan & Shortlist Approach

**Status:** DRAFT for review · **Gate:** organism + drug choice needs ⚕ microbiology sign-off
**Rule:** this phase pulls **counts only** (tiny metadata queries). **No genome downloads** until the
data plan is approved.

---

## 1. Goal
Pick, *from the data*, **one pathogen + a handful of antibiotics** that have **abundant and
class-balanced** resistant/susceptible genomes — enough to train *and* honestly test under a
phylogeny-aware split. Choosing blind is the single biggest scheduling risk in a 4-week sprint
(proposal §6.1); this removes it.

## 2. Data source & how we count
- **BV-BRC Data API** (REST), core `genome_amr` — one row per (genome, antibiotic) lab phenotype.
  Base URL: `https://www.bv-brc.org/api/genome_amr/`
- Fields we rely on: `genome_id`, `antibiotic`, `resistant_phenotype`
  (`Resistant` / `Susceptible` / `Intermediate`), `laboratory_typing_method`, `measurement`,
  and a taxon filter (`taxon_lineage_ids` preferred so sub-species/strains are included;
  fall back to `taxon_id` — **exact field verified on the first live call**, not assumed).
- **Counting method (reliable, no bulk pull):**
  1. Facet `genome_amr` by `antibiotic` for the organism → the list of drugs that have any data.
  2. For each drug, issue `limit(1)` count queries filtered by phenotype and read the total from
     the `Content-Range` response header → exact #Resistant, #Susceptible, #Intermediate.
- Volume: ~9 organisms × ~20 drugs × 3 phenotypes ≈ a few hundred tiny requests. Rate-limited and
  polite; responses are counts, not sequences.

## 3. Candidate organisms (starting set — the microbiologist prunes/extends)
Brief's four + common AMR-ML benchmark organisms. `taxon_id` to be confirmed on first call.

| Organism | taxon_id | Why a candidate |
|---|---|---|
| *Escherichia coli* | 562 | Huge BV-BRC AMR coverage; ref. Moradigaravand 2018 |
| *Klebsiella pneumoniae* | 573 | Carbapenem-resistance interest; well populated |
| *Salmonella enterica* | 28901 | Ref. Nguyen 2019 (NTS MIC prediction) |
| *Mycobacterium tuberculosis* | 1773 | Strong genotype→phenotype link (rpoB/katG/gyrA) |
| *Staphylococcus aureus* | 1280 | MRSA; clear determinants |
| *Acinetobacter baumannii* | 470 | MDR nosocomial |
| *Pseudomonas aeruginosa* | 287 | MDR; intrinsic + acquired resistance |
| *Neisseria gonorrhoeae* | 485 | Focused drug panel, good typing |
| *Streptococcus pneumoniae* | 1313 | Well-typed; balanced panels |

## 4. Selection criteria (how a drug qualifies for the shortlist)
Applied per (organism, drug), computed purely from the counts above:

| Criterion | Threshold (config-driven) | Rationale |
|---|---|---|
| **Abundance** | `min(#R, #S) ≥ 100` (`data.min_genomes_per_class`) | Enough minority-class genomes to survive a lineage-held-out test split. |
| **Balance** | minority fraction `min(R,S)/(R+S) ≥ 0.20`; prefer ≥ 0.30 | Avoid degenerate class imbalance; PR-AUC still reportable. |
| **Total labels** | `#R + #S ≥ 300` | Headroom for train/test after dropping Intermediate. |
| **Lineage headroom** *(proxy now)* | high raw genome count | True lineage diversity is only known after MLST (Phase 5); high N is the Phase-1 proxy. Flagged, not guaranteed. |
| **Label method** *(⚕ audit)* | prefer standardized MIC / broth microdilution | Cleaner phenotypes; microbiologist reviews `laboratory_typing_method` mix. |

`Intermediate` is **excluded from counts** for shortlisting (config `label_mapping.intermediate: drop`);
its handling is a ⚕ decision revisited in Phase 3.

## 5. Output — the shortlist we bring to the decision
For every candidate organism, a table of qualifying drugs:

```
Organism: <name> (taxon_id, total AMR rows)
| Antibiotic | #R | #S | #Intermediate | balance | total(R+S) | qualifies? |
```

Then a **ranked recommendation**: the 1–2 organisms with the most qualifying, balanced drugs, and a
proposed panel of **3–6 drugs** for the chosen organism. Written to
`results/reports/phase1_shortlist.md` + machine-readable `results/metrics/phase1_counts.csv`.

## 6. The decision (GATE)
The microbiologist + engineer jointly pick **one organism + the drug panel** from the shortlist.
Logged in `docs/decisions.md`; written into `config/config.yaml` (`organism`, `drugs`). Only then
does Phase 2 (actual genome download, size flagged first) begin.

## 7. What this phase explicitly does NOT do
- No genome/sequence downloads. No AMRFinderPlus DB download. No model training.
- No committing to an organism before the microbiologist reviews the counts.
- No fabricated numbers — if a query fails or a field name differs, we fix the query and report the
  real counts, never a placeholder.

## 8. How to run (once approved to execute the counts)
```bash
conda run -n amr-resistance-predictor python -m src.data.survey --config config/config.yaml
# add --dry-run to print the planned queries without hitting the network
```
Produces `results/metrics/phase1_counts.csv` and `results/reports/phase1_shortlist.md`.
