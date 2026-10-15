.PHONY: validate python-check bash-check r-parse r-tests functional-tests manifest

validate: python-check bash-check r-parse r-tests functional-tests
	python tools/validate_repository.py

python-check:
	python -m compileall -q scripts workflow

bash-check:
	@for script in workflow/*.sh; do bash -n "$$script" || exit $$?; done

r-parse:
	Rscript -e 'files <- list.files("scripts", pattern="[.]R$$", recursive=TRUE, full.names=TRUE); invisible(lapply(files, parse)); cat("R_PARSE_PASS:", length(files), "files\n")'

r-tests:
	Rscript tests/test_taxonomic_invariants.R

functional-tests:
	OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m unittest discover -s tests -p 'test_functional_*.py'

manifest:
	python tools/build_manifest.py

