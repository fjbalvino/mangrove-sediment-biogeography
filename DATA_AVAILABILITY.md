# Data availability

## Raw sequencing data

Public source studies are accessed through their NCBI SRA/ENA run accessions.
Raw shotgun metagenomic reads generated for the Mexican Atlantic component of
the study were submitted to NCBI SRA under BioProject
[PRJNA1502944](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1502944).

The final accession table will map every `run_accession` to BioSample,
BioProject, study, country, locality, marine realm, source publication and
whether the sequence data were newly generated or reused.

## Inputs for the taxonomic analyses

Place the external input files under `data/`, or set the environment variables
in `config/taxonomic.env.example` to their locations. The sample manifests are included in this repository and define cohort
membership explicitly. Each identifier matches the `run_accession` field;
newly generated samples may retain their study sample names until their public
accession mapping is included in the final provenance table.

| File | Contents and use | Configuration |
| --- | --- | --- |
| `data/taxonomic_counts.rds` | A phyloseq object with counts, taxonomy and sample metadata for 405 metagenomes. The supplied counts are already rarefied to 18,145 reads per metagenome. | `MAGENTA_PHYLOSEQ_405` |
| `data/taxonomic_core_counts.rds` | The 401-metagenome subset used for core and community-assembly analyses. Create it from the 405-metagenome object with `scripts/01_taxonomic_preprocessing/01H_SELECT_CORE_ANALYSIS_SAMPLES.R` and `config/samples/taxonomic_core_accessions.txt`; taxonomic features and counts within retained samples are preserved. | `MAGENTA_PHYLOSEQ_401`; `CORE_SAMPLE_MANIFEST` selects the accession manifest. |
| `data/environmental_metadata.csv` | Sample identifiers, coordinates, regional assignments, study identifiers and environmental covariates, joined to counts by `run_accession`. | `MAGENTA_METADATA_CSV` |
| `data/genus_phylogeny_GTDB_R226_bac120.newick` | The external bacterial genus phylogeny derived from the GTDB R226 bac120 reference tree, with 1,905 named tips and branch lengths. It is supplied independently of read classification. | `TREE_FILE` or the `--tree` option of `06A_ICAMP_GENUS_COMMUNITY_ASSEMBLY.R`. |

The diversity cohort is defined by
`config/samples/taxonomic_diversity_accessions.txt` (345 accessions), configurable
with `ALPHA_SAMPLE_MANIFEST`. Script `03A` retains those samples and standardizes
counts to the minimum input library size without replacement. With the supplied
object, every selected profile already contains 18,145 reads, so this operation
does not introduce a further reduction in depth.

The count object must retain `run_accession`, `REALM`, `Locality` and
`EcologicalClass` in its sample metadata and genus annotations in its taxonomy
table. Environmental metadata must provide `run_accession`, latitude/longitude
(`lat`, `lon`), locality, country, marine realm, province and ecoregion assignments,
and study or BioProject identifiers for the corresponding regional checks.

The Gaussian GAM uses these 14 environmental columns:

| Column | Variable |
| --- | --- |
| `human_footprint_2020` | Human footprint index (2020) |
| `bio1` | Annual mean temperature |
| `bio2` | Mean diurnal temperature range |
| `bio3` | Isothermality |
| `bio13` | Precipitation of the wettest month |
| `bio15` | Precipitation seasonality |
| `gbif_shannon_H_2km` | GBIF plant Shannon diversity within 2 km |
| `gbif_richness_species_2km` | GBIF plant species richness within 2 km |
| `soc_mean` | Soil organic carbon |
| `ph_mean` | Soil pH |
| `sand_mean` | Soil sand content |
| `silt_mean` | Soil silt content |
| `clay_mean` | Soil clay content |
| `bdod_mean` | Soil bulk density |

Use the archived metadata values and scales to reproduce the fitted models.
Missing environmental values are handled by the analysis-specific complete-case
filters recorded in each output audit.

## Reference input checksums

SHA-256 identifies the exact file bytes. Renaming a file does not change its
checksum; reserializing an R object or rewriting a CSV may do so.

| File | SHA-256 |
| --- | --- |
| `data/taxonomic_counts.rds` | `2c3123e46e8b35f240c86ad7b1a2ac01a660dce26248b328baeadce22bc72181` |
| `data/environmental_metadata.csv` | `b04e97e92b24cb8a892cf47548ad2027a780cf1b5583aa302fa0ff4de5420bd7` |
| `data/genus_phylogeny_GTDB_R226_bac120.newick` | `15761713cbfbe4edcf0600bb65afbea0f39cfe6a01a7b111a92aeeba59ffd3fd` |

The derived 401-metagenome object is identified by its source input and accession
manifest. The reference tree checksum is checked when the default tree file is
used. Alternative trees supplied with `--tree` are recorded in the output audit;
their taxonomic coverage and genus-name matches must be inspected separately.
Name matching is not a taxonomic synonym crosswalk or a test of phylogenetic
niche signal.

## Companion data archive

Large derived inputs are stored separately from GitHub. The companion archive
is intended to include the inputs above, the accession/provenance table, the
functional KO matrix, environmental assets or caches where redistribution is
permitted, primary and supplementary result tables, and a SHA-256 manifest.

A public archive DOI and version have not yet been added to this repository.
They will be recorded here, in the README and in `CITATION.cff` when available.

## Licensing and attribution

Repository code is MIT-licensed. Newly authored tabular documentation and
metadata will be released under the license stated in the Zenodo record.
Third-party sequence records, occurrence data and environmental products retain
their original attribution and licensing requirements; their sources and
versions must remain explicit in the archived provenance table.

