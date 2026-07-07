# Reading Resistance — one-command build.
# Processed artifacts regenerate from raw + code; data/raw is never a build target.
# Each target maps to a pipeline phase. Stubs raise until the phase is implemented.

.PHONY: all data features split train eval figures test clean help
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

eval:  ## Phase 8/9: metrics, calibration, honest benchmark
	$(PY) -m src.evaluation.evaluate --config $(CONFIG)

figures:  ## Phase 8/10: reliability curves, SHAP plots, benchmark table
	$(PY) -m src.evaluation.figures --config $(CONFIG)

test:  ## Run the test suite (data joins, feature builder, LEAKAGE test)
	pytest -q

all: data features split train eval figures  ## Full pipeline end-to-end

clean:  ## Remove regenerable artifacts (keeps data/raw and results/reports)
	rm -rf data/interim/* data/processed/* results/figures/* results/metrics/* results/models/*
