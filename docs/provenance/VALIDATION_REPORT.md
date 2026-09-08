# Repository validation report

Validated on **2026-09-08 UTC**.

## Package integrity

- 43 executable/helper files in the canonical `scripts/` tree, including the
  optional Bracken sensitivity script.
- Historical scripts and exact duplicates are excluded from the public
  repository; their selection audit is preserved in `DUPLICATE_AUDIT.tsv`.
- No byte-identical files remain duplicated inside the canonical tree.
- All canonical analytical copies are byte-identical to their selected uploaded
  source files; only their filesystem names/locations changed.
- `SCRIPT_INDEX.tsv`, `METHODS_CODE_CROSSWALK.tsv` and
  `DUPLICATE_AUDIT.tsv` have consistent tabular structure.
- `CITATION.cff` parses as valid YAML and contains seven authors.

## Syntax checks

- Python: all canonical `.py` files passed `compileall`.
- Bash: all three workflow scripts passed `bash -n`.
- R: the local packaging runtime does not contain `Rscript`, so a new local
  `parse()` sweep was not possible. The selected R scripts are unchanged from
  the supplied/previously executed versions. Run the following on Landa Lab
  before creating a public release:

```bash
find scripts -type f -name '*.R' -print0 | while IFS= read -r -d '' f; do
  /home/fjbalvino/.mamba/envs/rtool/bin/Rscript \
    -e 'parse(file=commandArgs(TRUE)[1])' "$f" || exit 1
done
```

## Authoritative checksums

| Component | SHA-256 |
|---|---|
| Primary iCAMP script 2.2.0 | `423595c69897b51ab027613c5fdda8d09bf30d62d21b54368056d159e5700d1c` |
| Functional regional script 180.2.3 | `e915927f23170c2631ebebabf1b441adf48a8ece7eca2a0b4d4a10d783d8f026` |
| Corrected KO annotation snapshot builder | `d28c94b6d4a1f9e28c17a5e48de8206d34fac35af0924c4d507c7509e36d4b4b` |
| Corrected KO CLR/network preparation | `7d27d9489df6125db6ede65d9672447a1fab5d0996073c11ae288ba4d367d33a` |

The complete repository inventory is in `MANIFEST.sha256`.

## Evidence from completed server runs

- Primary iCAMP: `PASS`, 400 samples, 1,892 genera, 79,800 total pairs,
  24,692 within-realm pairs and 999 randomizations.
- Functional regional analysis 180.2.3: `PASS` with normal PID closure,
  `Metadata_2026_completa.csv`, three explicit unmatched-accession exclusions,
  191 retained profiles and 17,482 non-zero KOs.

These completed-run statements are based on the status/log outputs supplied
during the project audit. Raw result directories are not embedded in this code
archive and should be deposited separately with checksums.

## Release blockers

1. Add and validate the missing MEGAHIT-to-KO scripts in block 07.
2. Freeze external datasets/caches and add their checksums and retrieval dates.
3. Deposit the final accession table and derived-data bundle in Zenodo.
4. Run the configured GitHub Actions R parse check and obtain a green status.
