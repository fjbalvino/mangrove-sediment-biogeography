# Gene biogeography configuration

These small files contain the frozen operational rules, curated process descriptions, nitrogen reaction labels and 368 queried KO display labels. They are configuration, not a downloaded functional annotation database and not sample-level measurements. `KO_labels.tsv` explicitly excludes all historical counts, prevalence and sample fields.

See `docs/GENE_BIOGEOGRAPHY.md` for the input contract and execution stages. New matrices and metadata must be supplied externally. The provenance record identifies the reviewed source package and original source-file hashes before portability changes. Requirements pin the Python analysis environment. Never add matrices, full metadata or generated HTML containing embedded data to this directory.
