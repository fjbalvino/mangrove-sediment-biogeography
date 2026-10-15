# Data availability

The Figure 6 rerun uses the same 191-profile functional cohort documented
below. The companion `MAGENTA_functional_network.zip` additionally supplies the
frozen KO dictionary, restricted historical KEGG MODULE mapping, corrected
network scripts, final tables and execution logs. Input checksums are in
`config/network_input_checksums.json`; mapping limitations and reproduction
instructions are in `docs/FUNCTIONAL_NETWORK_ANALYSES.md`.

## Taxonomic reference

The current reference contains 341 profiles. The companion data archive supplies
these files, used directly without another classification or rarefaction step:

| File | Contents | Configuration |
| --- | --- | --- |
| `data/phyloseq_taxonomic.rds` | Counts, taxonomy and sample metadata; 341 profiles, 24,045 input features, 18,145 reads per profile | `MAGENTA_PHYLOSEQ` |
| `data/metadata_taxonomic.csv` | Environmental covariates, regional assignments, sampling units and source-material annotations | `MAGENTA_METADATA_CSV` |
| `config/samples/taxonomic_accessions.txt` | Reference accession manifest, also present in this repository | `ALPHA_SAMPLE_MANIFEST`; `CORE_SAMPLE_MANIFEST` |

SHA-256 of the supplied files:

| File | SHA-256 |
| --- | --- |
| `phyloseq_taxonomic.rds` | `bc28139721edac3eaf9e4cb60f4333d43943d15ee6afef1a46668310c882ac84` |
| `metadata_taxonomic.csv` | `4179f08b18c133cd2da7355c101125a501a9c03a2e6243fca55186d5745f7ce1` |

`run_accession` is the join key. `sampling_site_id` groups identical coordinates
rounded to five decimal places, with the published Mud Island sampling site
represented as a single unit across campaigns and depths. `core_locality`
identifies named localities for core pooling; unresolved exact localities remain
missing. `study_id` combines the eight sample-level Texas BioProjects into their
common published study. Other public studies use their archive study accession.

Missing coordinates and covariates remain missing. Each analysis records its
required complete-case subset. Source-material review fields distinguish
supported sediment/soil records, explicitly different material and unresolved
records; they are not silently converted into further exclusions.

## Environmental variables

The Shannon GAM includes human footprint (2020); annual mean temperature (BIO1);
mean diurnal temperature range (BIO2); isothermality (BIO3); precipitation of the
wettest month (BIO13); precipitation seasonality (BIO15); GBIF plant Shannon
diversity and species richness within 2 km; and SoilGrids organic carbon, pH,
sand, silt, clay and bulk density. The metadata column `phh2o_mean` is accepted as
`ph_mean` by model code. SoilGrids pH is stored as pH × 10 and organic carbon in
dg/kg; figure axes state those scales. Use the supplied covariate values to
reproduce the models.

## Separate iCAMP analysis

The retained iCAMP estimates use their original 400-profile, 1,892-genus analysis,
24,692 within-realm pairs and 999 randomizations. They were not recalculated from
the 341-profile reference. The verified process table is included under
`results/icamp/`; the companion archive also includes the original supplementary
workbook and audit text supporting it.

Independent reproduction requires the original count object, embedded metadata
and GTDB R226 bacterial genus reference phylogeny. The count object is supplied
to script 06A through `ICAMP_PHYLOSEQ` (default
`data/icamp_taxonomic_counts.rds`); the source accession selection is recorded in
`config/samples/icamp_accessions.txt`. The original phylogeny has 1,905 named tips
and SHA-256 `15761713cbfbe4edcf0600bb65afbea0f39cfe6a01a7b111a92aeeba59ffd3fd`.
The original tree and fitted pairwise RDS are not included in the current archive.
Do not substitute the 341-profile object when reporting these retained estimates.

## Functional atlas inputs

The compact companion package `Functional_analysis_package.zip` contains the
three files below under `functional_analysis/inputs/`. Copy that directory to
`data/functional/`, or set `FUNCTIONAL_INPUT_DIR` to its extracted location.
The input files remain outside the code repository.

| File | Contents |
| --- | --- |
| `functional_KO_counts_original.tsv.gz` | Source KO matrix: 17,508 KOs and 194 profiles, including three air controls |
| `metadata_functional.csv` | Corrected realm, locality, site and study metadata for the 191 retained profiles |
| `taxonomic_matched_counts.npz` | Paired taxonomic counts, identifiers and features; 18,145 reads per profile |

The analysis selects the 191 accessions in
`config/samples/functional_accessions.txt`; only the documented air controls
may be omitted from the source KO matrix. The local sample identifier `TX1_3`
is retained. Twenty-six KOs absent from all retained profiles are removed.
The paired taxonomic matrix contains 17,906 nonzero features.

SHA-256 checksums, filenames and byte sizes are recorded in
`config/functional_input_checksums.json`. The decompressed source KO matrix has
SHA-256 `57455728bd3d1e7b56a7d83471e8eed2e9f0b547a5e50331d314ef270f629de1`.
Counts represent protein-cluster members assigned to KO terms; each listed KO
receives the full protein count. They are not read-mapping abundances.

The source protein matrix contained 258 profiles; its 64-entry exclusion list
produced the 194-profile KO matrix. Of the 341 taxonomic profiles, 195 were in
the protein matrix before sample filtering and four were on its exclusion list,
leaving 191 matched profiles. The numeric rule used to construct that exclusion
list is not yet independently verified: `metagenome_totals.tsv` and
`process_matrix.py` are required to establish the proposed 2,000-protein cutoff.
The current Results therefore do not assert that cutoff.

## Functional prevalence and regional differentiation (Figure 5)

The companion `MAGENTA_Figure_5.zip` provides the original KO matrix and corrected
functional metadata above, plus `KO_persistence_classes.tsv` (the current atlas
classification used for independent reconciliation) and `KO_function_dictionary.tsv`
(a frozen KO description dictionary). Their hashes and sizes are in
`config/functional_regional_input_checksums.json`. The same 191-accession manifest
governs selection; no taxonomic object is required. The package also includes full
per-selection results, RDS and NPZ analysis objects, assembled and separate figure panels,
two Word reports with images, logs, exact commands for landalab and the proposed
repository changes. Inputs and full analysis objects remain outside GitHub.

Classification uses 191 profiles and 26 named realm-localities. Balanced comparisons
use 188 eligible profiles in Central Indo-Pacific, Tropical Atlantic and Western
Indo-Pacific; Tropical Eastern Pacific has only one locality and remains descriptive.
The primary and alternative selections contain 18 and 12 profiles, respectively.
The selected-process catalog, its inclusion audit, process-score effect sizes, locality-turnover summaries and base-R verification are included. The map retains all 26 locality centroids but displays only the 23 with coordinates (184 profiles); the three unlocated units remain in the global reference.
Do not substitute the 184-profile coordinate-complete functional subset or the
341-profile taxonomic and 400-profile iCAMP cohorts for these units.

## Raw sequences and archive

Public source studies are identified by NCBI SRA/ENA run, BioSample and BioProject
accessions. Newly generated Mexican Atlantic reads were submitted under
[PRJNA1502944](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1502944); local study
identifiers are retained until the complete public crosswalk is supplied.

Large input objects and full results are distributed in the companion archive,
not in GitHub. No public archive DOI has yet been assigned. Third-party sequence,
occurrence, environmental and reference-tree data retain their original licenses
and accession-level attribution.
