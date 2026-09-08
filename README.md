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

## Where the public research objects live

| Resource | Public home | Contents |
|---|---|---|
| Source code | This GitHub repository | Scripts, configuration templates, workflow launchers and provenance |
| Frozen derived data | Zenodo release (DOI pending) | Curated metadata, accession table, count/KO matrices, phyloseq objects, the GTDB R226 genus tree, principal result tables and checksums |
| Raw reads | NCBI SRA | Public source studies and newly generated reads under BioProject [PRJNA1502944](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1502944) |

Raw FASTQ files, large databases and complete analysis directories are not
stored in GitHub. This separation keeps the software reviewable while Zenodo
provides an immutable, DOI-bearing data release and SRA preserves raw reads.

## Workflow

| Block | Objective | Scripts |
|---|---|---|
| 00 | Landa Lab preflight | `scripts/00_PREFLIGHT_LANDALAB.py` |
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

The authoritative file list and the Methods-to-code mapping are in
[`docs/provenance/SCRIPT_INDEX.tsv`](docs/provenance/SCRIPT_INDEX.tsv) and
[`docs/provenance/METHODS_CODE_CROSSWALK.tsv`](docs/provenance/METHODS_CODE_CROSSWALK.tsv).

## Analysis cohorts

These cohorts support different questions and must not be interchanged.

| Analysis | Audited cohort |
|---|---|
| Taxonomic Set401 | 401 metagenomes after the four declared Giralia exclusions |
| Primary iCAMP | 400 metagenomes, 1,892 genera, 24,692 within-realm pairs and 999 null randomizations |
| Functional atlas | 194 assembly-quality-filtered profiles and 17,508 KOs |
| Regional functional analysis | 191 metadata-matched profiles and 17,482 non-zero KOs after three declared exclusions |
| Functional network | 194 profiles and 17,508 KOs before network feature selection |

## Quick start on Landa Lab

1. Clone the repository and enter it.

   ```bash
   git clone https://github.com/fjbalvino/mangrove-sediment-biogeography.git
   cd mangrove-sediment-biogeography
   ```

2. Create a private local configuration from the template and verify every
   path. The local file is ignored by Git.

   ```bash
   cp config/landalab.env.example config/landalab.env
   ```

3. Run the read-only preflight.

   ```bash
   python scripts/00_PREFLIGHT_LANDALAB.py \
     --config config/landalab.env \
     --scope all \
     --out preflight.json
   ```

4. Inspect the dry-run workflow before executing long jobs.

   ```bash
   bash workflow/01_RUN_TAXONOMIC_PREPROCESSING.sh
   bash workflow/02_RUN_ENVIRONMENTAL_COVARIATES.sh
   ```

5. Add `--execute` only after the printed commands, paths and checksums have
   been reviewed.

Detailed instructions are in
[`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md). A beginner-oriented
publication checklist is in
[`docs/GITHUB_ZENODO_GUIDE_ES.md`](docs/GITHUB_ZENODO_GUIDE_ES.md).

## Reproducibility conventions

- `run_accession` is the canonical sample identifier.
- Each analysis writes to a time-stamped directory and records status, logs,
  parameters, session information and checksums where implemented.
- The primary iCAMP inference uses within-realm comparisons; between-realm
  pairs are retained for audit but are not pooled into the primary estimate.
- The exact GTDB R226 genus tree used by iCAMP has SHA-256
  `15761713cbfbe4edcf0600bb65afbea0f39cfe6a01a7b111a92aeeba59ffd3fd`.
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
than importing changing files from MAGENTA. See
[`docs/RELATIONSHIP_TO_MAGENTA.md`](docs/RELATIONSHIP_TO_MAGENTA.md).

## Citation and license

Citation metadata are provided in [`CITATION.cff`](CITATION.cff). Code is
released under the [MIT License](LICENSE). Data licenses and third-party terms
are documented per archived object in the Zenodo record; source datasets retain
their original licenses and accession-level attribution.

