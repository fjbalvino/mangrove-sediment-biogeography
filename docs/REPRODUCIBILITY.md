# Reproducibility guide

## 1. Obtain the code

Use a tagged release for manuscript reproduction. The moving `main` branch is
for development and should not be cited as an immutable version.

## 2. Obtain the data

Download the matching Zenodo data release and verify its SHA-256 manifest.
Keep the archive outside the Git repository. Raw reads remain in SRA/ENA and
can be reacquired from the frozen accession table.

## 3. Configure the server

Copy `config/landalab.env.example` to `config/landalab.env`, then set paths to
the downloaded data, software environments and frozen reference assets. Never
commit this local configuration.

The repository records Python requirements and audited R package versions in
`config/`. Exact external command-line tool and database versions must also be
reported in the Zenodo provenance bundle.

## 4. Validate before computing

Run:

```bash
make validate

python scripts/00_PREFLIGHT_LANDALAB.py \
  --config config/landalab.env \
  --scope all \
  --out preflight.json
```

The first command checks code syntax and repository integrity. The second
checks the Landa Lab paths and inputs without changing data.

## 5. Execute in scientific order

Blocks 01 and 02 have dry-run orchestrators under `workflow/`. Later blocks are
executed through their individual scripts because their tested command-line
contracts differ. The definitive command, parameters, seed and environment for
each manuscript result must be retained in its result directory and Zenodo
provenance bundle.

Block 07 is currently pending and must be completed before the first archival
release claims end-to-end reconstruction of the KO matrix.

## 6. Preserve each run

Use a time-stamped, non-overwriting output directory. Preserve, where
implemented:

- `STATUS.txt`;
- logs;
- input and parameter tables;
- software/session information;
- validation tables;
- principal result tables and figures; and
- checksums for frozen inputs and outputs.

Only a documented `PASS` run should support manuscript values. Failed and
sensitivity runs should remain distinguishable from the primary analysis.

## 7. Build the public data bundle

`workflow/99_BUILD_EXPORT_BUNDLE.sh` inventories the principal Landa Lab inputs
and is dry-run by default. Review its file list, add the final functional
annotation products and permitted environmental caches, then create a versioned
bundle for Zenodo. Do not include raw FASTQ files or databases.

