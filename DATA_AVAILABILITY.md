# Data availability

## Raw sequencing data

Public source studies are accessed through their NCBI SRA/ENA run accessions.
Raw shotgun metagenomic reads generated for the Mexican Atlantic component of
the study were submitted to NCBI SRA under BioProject
[PRJNA1502944](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1502944).

The final accession table will map every `run_accession` to BioSample,
BioProject, study, country, locality, marine realm, source publication and
whether the sequence data were newly generated or reused.

## Frozen derived data

The companion Zenodo record will contain the immutable inputs needed to rerun
the analyses without placing large files in GitHub:

- curated and environmental metadata;
- the accession/provenance table;
- taxonomic count matrices and the canonical phyloseq objects;
- the functional KO matrix;
- the exact GTDB R226 genus tree used by iCAMP;
- frozen environmental assets or caches where redistribution is permitted;
- primary and supplementary result tables required to reproduce reported
  values; and
- a SHA-256 manifest.

The Zenodo DOI and version will be added here, to the README and to
`CITATION.cff` before the public manuscript release.

## Licensing and attribution

Repository code is MIT-licensed. Newly authored tabular documentation and
metadata will be released under the license stated in the Zenodo record.
Third-party sequence records, occurrence data and environmental products retain
their original attribution and licensing requirements; their sources and
versions must remain explicit in the archived provenance table.

