.PHONY: validate python-check bash-check r-parse r-tests functional-tests gene-environment-tests portability-tests regional-tests manifest

validate: python-check bash-check r-parse r-tests functional-tests gene-environment-tests regional-tests portability-tests
	python tools/validate_repository.py

python-check:
	python -m compileall -q scripts workflow tools tests

bash-check:
	@find workflow scripts -name '*.sh' -print0 | xargs -0 -n 1 bash -n

r-parse:
	Rscript -e 'files <- list.files("scripts", pattern="[.]R$$", recursive=TRUE, full.names=TRUE); invisible(lapply(files, parse)); cat("R_PARSE_PASS:", length(files), "files\n")'

r-tests:
	Rscript tests/test_taxonomic_invariants.R

functional-tests:
	OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m unittest discover -s tests -p 'test_functional_*.py'

gene-environment-tests:
	OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m unittest discover -s tests -p 'test_gene_environment_*.py'

regional-tests:
	Rscript --vanilla tests/test_regional_functional_statistics.R

portability-tests:
	python -m unittest discover -s tests -p 'test_workflow_*.py'

manifest:
	python tools/build_manifest.py

