# Data directory

Only small, publication-safe tables belong in this directory. Large or binary
data products are released through Zenodo and raw reads through SRA/ENA.

- `accessions/`: final small accession and provenance tables.
- `zenodo/`: inventories and checksums describing the matching Zenodo release.

Do not commit FASTQ files, count matrices, RDS objects, environmental rasters,
Kraken2 databases, cached service responses or generated analysis directories.

