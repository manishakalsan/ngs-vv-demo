.PHONY: help data run test unit trace demo-change giab clean

help:            ## list targets
	@grep -E '^[a-z-]+:.*##' $(MAKEFILE_LIST) | sed 's/:.*##/\t/'

data:            ## regenerate the synthetic dataset (deterministic)
	python scripts/simulate_test_data.py

run:             ## run the pipeline on synthetic data
	snakemake --cores 2

unit:            ## unit tests with coverage
	pytest tests/unit --cov=vvdemo --cov-report=term-missing

test:            ## full verification suite
	pytest --cov=vvdemo --junitxml=reports/junit.xml

trace: test      ## rebuild traceability matrix (release gate)
	python scripts/build_traceability.py

demo-change:     ## show regression suite catching CR-001
	bash scripts/demo_regression_catch.sh

giab:            ## fetch GIAB slice and run benchmark
	bash scripts/fetch_giab.sh
	snakemake --cores 4 --configfile config/giab_hg002_chr20.yaml

clean:           ## remove run outputs
	rm -rf results results_* .snakemake reports .pytest_cache .coverage
