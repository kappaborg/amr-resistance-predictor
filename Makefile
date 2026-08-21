# Reading Resistance — one-command build.
# Processed artifacts regenerate from raw + code; data/raw is never a build target.
# Each target maps to a pipeline phase. Stubs raise until the phase is implemented.

.PHONY: all data features split train eval rigor figures leakage dca riskcov extscope extfetch extscore organisms models test clean help
.DEFAULT_GOAL := help

CONFIG := config/config.yaml
PY     := python

help:  ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN{FS=":.*?## "}{printf "  \033[36m%-10s\033[0m %s\n", $$1, $$2}'

data:  ## Phase 2: acquire genomes + phenotypes from BV-BRC (flags size first)
	$(PY) -m src.data.acquire --config $(CONFIG)

features:  ## Phase 4: annotate genomes -> determinant feature matrix
	$(PY) -m src.features.build --config $(CONFIG)

split:  ## Phase 5: assign lineages + phylogeny-aware train/test split
	$(PY) -m src.split.make_split --config $(CONFIG)

train:  ## Phase 6/7: train per-drug models (baseline + gradient boosting)
	$(PY) -m src.models.train --config $(CONFIG)

eval:  ## Phase 8/9: regulator-grade metrics (VME/ME/CA), calibration+Brier, CIs, DeLong — all 8 orgs
	$(PY) -m src.evaluation.clinical_rigor --config $(CONFIG)

rigor: eval  ## Alias for eval (regulator-grade evaluation)

figures:  ## Phase 8/10: reliability curves, SHAP plots, benchmark table, 8-organism heatmap
	$(PY) -m src.evaluation.figures --config $(CONFIG)

leakage:  ## Quantify population-structure leakage: random vs lineage vs temporal AUC
	$(PY) -m src.evaluation.leakage_delta --config $(CONFIG)

dca:  ## Decision-curve analysis (clinical net benefit) — needs eval first (pooled_predictions.pkl)
	$(PY) -m src.evaluation.decision_curve

riskcov:  ## Risk-coverage curves + AURC for the defer-to-lab abstention — needs eval first
	$(PY) -m src.evaluation.risk_coverage

extscope:  ## External-validation Step 1: scope an independent NCBI PD cohort (metadata only, ~179 MB)
	$(PY) -m src.evaluation.external_validation

extfetch:  ## External-validation Step 2a: download + AMRFinderPlus + mlst the external cohort (resumable)
	$(PY) -m src.evaluation.external_fetch_annotate

extscore:  ## External-validation Step 2b: score FROZEN models on the external cohort (ALL + ST-novel)
	$(PY) -m src.evaluation.external_score

organisms:  ## Annotate + MLST + train one organism (ORG=paeruginosa|senterica|efaecium|...)
	$(PY) -m src.organism_pipeline --organism $(ORG)

models:  ## Fit + save deployable per-drug models for all eight organisms (demo)
	$(PY) -m src.models.save_models --config $(CONFIG)

test:  ## Run the test suite (data joins, feature builder, per-organism LEAKAGE test)
	pytest -q

all: data features split train eval figures  ## Full pipeline end-to-end (K. pneumoniae panel)

clean:  ## Remove regenerable artifacts (keeps data/raw and results/reports)
	rm -rf data/interim/* data/processed/* results/figures/* results/metrics/* results/models/*
