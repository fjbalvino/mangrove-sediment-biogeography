# Contributing

This repository is the frozen analytical companion to a scientific study.
Changes must therefore preserve both reproducibility and provenance.

## Before proposing a change

1. Open an issue describing the scientific or technical reason for the change.
2. Create a short branch from `main`.
3. Keep canonical filenames. Record changes in Git history and the change log;
   preserve released data archives instead of adding numbered source duplicates.
4. Do not commit raw reads, large matrices, RDS files, databases, credentials or
   generated result directories.
5. Regenerate source hashes with `make manifest`, then run `make validate` before opening a pull request.

## Script requirements

- Keep the workflow prefix and use an objective-based filename.
- Record input paths, parameters, package versions, random seeds and output
  paths.
- Write a terminal `PASS` or `FAIL` state where the script implements run
  directories.
- Use `run_accession` as the sample key and validate one-to-one joins.
- Preserve the distinction among the 341-profile taxonomic reference, 400-profile iCAMP analysis,
  191-profile functional atlas and the separately documented regional/network
  analyses. Functional geographic comparisons use 184 profiles from 57 sites.
- Update `docs/provenance/SCRIPT_INDEX.tsv` when a canonical script changes.

## Pull requests

Describe the affected manuscript statement, inputs, expected output and
validation performed. A result-changing modification must include a comparison
against the prior canonical output.


Keep datasets, numerical outputs, maps and figures outside Git. The source export and validation tools enforce this boundary. Text report sources and small accession/rule configurations may remain.
