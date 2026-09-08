.PHONY: validate python-check bash-check r-parse manifest

validate: python-check bash-check r-parse
	python tools/validate_repository.py

python-check:
	python -m compileall -q scripts

bash-check:
	bash -n workflow/*.sh

r-parse:
	Rscript -e 'files <- list.files("scripts", pattern="[.]R$$", recursive=TRUE, full.names=TRUE); invisible(lapply(files, parse)); cat("R_PARSE_PASS:", length(files), "files\n")'

manifest:
	python tools/build_manifest.py

