# Global biogeography of mangrove sediment microbiomes

[![Repository validation](https://github.com/fjbalvino/mangrove-sediment-biogeography/actions/workflows/validate.yml/badge.svg)](https://github.com/fjbalvino/mangrove-sediment-biogeography/actions/workflows/validate.yml)
[![License: MIT](https://img.shields.io/badge/Code%20license-MIT-blue.svg)](LICENSE)
[![Status: work in progress](https://img.shields.io/badge/status-work%20in%20progress-orange.svg)](#release-status)

Reproducible code and documentation for the global analysis of mangrove
sediment microbiomes developed within the MAGENTA project. The repository is
organized by scientific objective, not by historical figure number.

The code supports taxonomic profiling, environmental covariates, alpha and beta
diversity, taxonomic biogeography, iCAMP community-assembly inference,
functional biogeography and KO co-variation network analyses.

## Workflow

| Block | Objective | Scripts |
|---|---|---|
| 01 | Metadata, reads, FastQC, fastp, Kraken2, taxonomic matrix and phyloseq | `01A`–`01H` |
| 02 | Human Footprint, WorldClim, GBIF, SoilGrids and MEOW | `02A`–`02G` |
| 03 | Hill diversity and environmental GAMs | `03A`–`03D` |
| 04 | NMDS, distance decay, PERMANOVA and taxonomic dbRDA/VPA | `04A`–`04D` |
| 05 | Taxonomic core, realm restriction and core/non-core turnover | `05A`–`05D` |
| 06 | Genus-level iCAMP and robustness to sampling structure | `06A`–`06B` |
| 07 | MEGAHIT → Prodigal → MMseqs2 → eggNOG → KO matrix | **Pending** |
| 08 | Functional persistence, PCA and taxonomic–functional coupling | `08A`–`08C` |
| 09 | Functional prevalence and regional differentiation | `09A` |
| 10 | KO annotation and consensus co-variation network | `10A`–`10G` |

## Analysis cohorts

These cohorts support different questions and must not be interchanged.

| Analysis | Audited cohort |
|---|---|
| Taxonomic Set401 | 401 metagenomes |
| Primary iCAMP | 400 metagenomes, 1,892 genera, 24,692 within-realm pairs and 999 null randomizations |
| Functional atlas | 194 assembly-quality-filtered profiles and 17,508 KOs |
| Regional functional analysis | 191 metadata-matched profiles and 17,482 non-zero KOs after three declared exclusions |
| Functional network | 194 profiles and 17,508 KOs before network feature selection |

## Using the workflow

Users must download or clone the repository and adapt paths, software
environments, database locations, thread and memory settings, and job-submission
procedures to their own computer or computing cluster before execution.

## Reproducibility conventions

- `run_accession` is the canonical sample identifier.
- Each analysis writes to a time-stamped directory and records status, logs,
  parameters, session information and checksums where implemented.
- The primary iCAMP inference uses within-realm comparisons; between-realm
  pairs are retained for audit but are not pooled into the primary estimate.
- Live GBIF and SoilGrids queries and gridded covariates must be represented by
  frozen caches/assets in the Zenodo release for byte-identical reproduction.
- The KO network describes robust co-variation among functional features; it
  is not evidence of direct molecular interaction.

## Release status

The canonical scripts are organized and repository validation is configured.
Block 07, which constructs the KO matrix from assemblies, is deliberately
marked pending until its scripts are supplied and audited. The first public
archival release should be tagged only after Block 07, Zenodo checksums and the
final accession table are complete.

## Relationship to MAGENTA

This is the manuscript-specific, auditable workflow derived from the broader
[MAGENTA project repository](https://github.com/fjbalvino/magenta). It freezes
the scripts and cohort contracts used for the global biogeography study rather
than importing changing files from MAGENTA.

## Citation and license

Citation metadata are provided in [`CITATION.cff`](CITATION.cff). Code is
released under the [MIT License](LICENSE). Data licenses and third-party terms
are documented per archived object in the Zenodo record; source datasets retain
their original licenses and accession-level attribution.
