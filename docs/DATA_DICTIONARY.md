# Public data dictionary

## Accession table

The publication accession table should contain one row per metagenomic run and
use `run_accession` as its unique key.

| Field | Meaning |
|---|---|
| `run_accession` | Canonical SRA/ENA run accession |
| `sample_accession` | BioSample or equivalent public sample accession |
| `bioproject_accession` | NCBI BioProject accession |
| `study_id` | Harmonized study/BioProject grouping used in sensitivity analyses |
| `new_sequence_data` | Whether the reads were generated for this study |
| `country` | Harmonized country |
| `locality` | Harmonized sampling locality |
| `latitude`, `longitude` | Public sampling coordinates and coordinate precision |
| `marine_realm` | Harmonized MEOW marine realm |
| `taxonomic_included` | Inclusion in the 401-sample taxonomic cohort |
| `icamp_included` | Inclusion in the 400-sample iCAMP cohort |
| `functional_atlas_included` | Inclusion in the 194-profile functional cohort |
| `functional_regional_included` | Inclusion in the 191-profile regional functional cohort |
| `exclusion_reason` | Explicit reason when any analysis excludes the run |
| `source_citation` | Dataset or article attribution |

Every derived matrix must state whether samples are rows or columns, its
feature identifier namespace, units, filtering rules and transformation state.

## Functional terminology

Functional tables represent KEGG orthologues (KOs), modules or derived
functional features. They should not be labeled as counts of “unique genes.”
Network edges represent reproducible co-variation under the documented
procedure, not direct biochemical interactions.

