# Data availability

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

## Raw sequences and archive

Public source studies are identified by NCBI SRA/ENA run, BioSample and BioProject
accessions. Newly generated Mexican Atlantic reads were submitted under
[PRJNA1502944](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1502944); local study
identifiers are retained until the complete public crosswalk is supplied.

Large input objects and full results are distributed in the companion archive,
not in GitHub. No public archive DOI has yet been assigned. Third-party sequence,
occurrence, environmental and reference-tree data retain their original licenses
and accession-level attribution.
