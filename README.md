# Mangrove microbiome biogeography — manuscript analysis code

Code and documentation accompanying **A global metagenomic atlas reveals shared functional potential across geographically structured mangrove microbiomes**.

## Overview

Mangrove microorganisms contribute to carbon and nutrient transformations across contrasting coastal environments. MAGENTA (**MAngrove GENe caTAlogue**) brings together read-based taxonomic profiles, gene-centred functional profiles and environmental metadata to ask how geographically different communities retain shared functional potential.

This repository accompanies the manuscript and contains its analysis scripts, configuration templates, accession manifests and documentation. **Input datasets, reference databases, numerical result tables, fitted objects and figures are distributed separately.** The companion Zenodo data and figure deposit is being prepared; no DOI is available yet. See [data availability](DATA_AVAILABILITY.md).

## Manuscript summary

The atlas combines **341 taxonomic profiles across five marine realms** with functional profiles from **191 of those metagenomes**, comprising **17,482 KEGG orthologues (KOs)**. Taxonomic diversity and composition varied geographically and with environmental setting. Core-associated taxa comprised 27.4% of eligible features and accounted for 62.6% of counts assigned to those features. Absolute-core, near-global and realm-core KOs together represented 14.5% of detected KOs and 86.2% of KO annotation counts. Taxonomic and functional dissimilarities were strongly correlated, while functional similarity declined more slowly with geographic distance. Co-abundance analyses identified recurring functional modules, with membership sensitive to sample inclusion. These findings describe shared encoded potential alongside regional differentiation; study differences limit attribution to ecological mechanisms.

Functional counts represent protein-cluster members assigned to KOs. They do not measure expression, metabolic activity or flux. A separate **400-profile iCAMP analysis** must not be confused with the taxonomic or functional reference cohorts.

## Pipeline workflow

| Stage | Purpose | Source scripts | Manuscript context |
| --- | --- | --- | --- |
| 01 | Acquire metadata and reads; QC, classify and select profiles | [01_taxonomic_preprocessing](scripts/01_taxonomic_preprocessing) | Atlas construction; Figure 1 |
| 02 | Reconstruct environmental and geographic covariates | [02_environmental_covariates](scripts/02_environmental_covariates) | Environmental setting; Figure 2 |
| 03 | Diversity, GAMs, composition, distance decay and taxonomic core | [03](scripts/03_alpha_diversity), [04](scripts/04_beta_diversity), [05](scripts/05_taxonomic_core_biogeography) | Figures 2–3 |
| Separate cohort | Community assembly and robustness | [06_community_assembly](scripts/06_community_assembly) | Figure 3; Table S8 |
| Pending upstream block | Assembly → proteins → clustering → KO matrix | [07_functional_annotation_PENDING](scripts/07_functional_annotation_PENDING) | Methods; original validated code still needed |
| 04 | Functional atlas, regional tests and taxonomic–functional coupling | [08_functional_atlas](scripts/08_functional_atlas) | Figure 4 |
| 05 | KO co-abundance network and sensitivity analyses | [10_functional_network](scripts/10_functional_network) | Figure 5 |
| 06 | Prevalence, persistence and balanced regional differentiation | [09_functional_regional_differentiation](scripts/09_functional_regional_differentiation) | Figure 6 |

The [manuscript guide](docs/MANUSCRIPT_GUIDE.md) maps every result section, supplementary table group and analysis cohort to the code. Existing script names are retained to preserve provenance; the execution guide follows the current manuscript order. This is a read-based/gene-centred workflow, with no MAG-binning stage.

## Prerequisites and installation

Use Bash, Python and R on Linux or a compatible Unix environment. The [installation guide](docs/INSTALLATION.md) separates dependencies for statistical analyses from tools and external databases needed for raw-read processing. Resource requirements depend on the chosen stage; raw-read classification and assembly need substantially more storage and memory than documentation checks.

```bash
git clone https://github.com/fjbalvino/mangrove-sediment-biogeography.git
cd mangrove-sediment-biogeography
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r config/functional_requirements.txt
cp config/project.env.example config/project.env
```

Edit `config/project.env` for your local inputs, outputs and executables. Obtain the required external inputs before executing analyses. Cloning this repository alone does not supply them.

## Usage

Preview the main entry points without running analyses or creating output directories:

```bash
bash workflow/01_RUN_TAXONOMIC_PREPROCESSING.sh
bash workflow/02_RUN_ENVIRONMENTAL_COVARIATES.sh
bash workflow/03_RUN_TAXONOMIC_ANALYSES.sh --core
bash workflow/04_RUN_FUNCTIONAL_ANALYSES.sh
bash workflow/06_RUN_FUNCTIONAL_REGIONAL.sh
```

Add `--execute` to a shell workflow after checking its configuration. The network driver has its own `--dry-run` option. See [00_README_HOW_TO_RUN](workflow/00_README_HOW_TO_RUN.md) for the complete commands, inputs and independent iCAMP execution.

Paths are configurable and do not depend on an institutional server. The supplied recipes retain the manuscript's cohort manifests, scientific thresholds and validation checks. Adapting them to another study requires reviewing those definitions as well as replacing the inputs.

## Repository and output directories

| Directory | Contents |
| --- | --- |
| `config/` | Portable settings, dependencies, cohort IDs, input hashes and analysis rules |
| `scripts/` | Numbered scientific analysis modules |
| `workflow/` | Execution guide and orchestration entry points |
| `docs/` | Methods, manuscript context and text sources for report builders |
| `tests/`, `tools/` | Statistical regression checks, validation and code export |
| `data/` | Local inputs only; excluded from Git except its README |
| `results/`, `figures/` | Locally generated outputs; excluded from Git except the results README |

Tiny accession manifests and curated process definitions remain as configuration, not biological count databases. Figure-generation code is included; image files are not.

## Validation and code export

```bash
make validate
bash workflow/99_BUILD_EXPORT_BUNDLE.sh --output ../manuscript_code.zip
bash workflow/99_BUILD_EXPORT_BUNDLE.sh --output ../manuscript_code.zip --execute
```

Validation checks source syntax, statistical contracts, workflow portability, the code-only boundary and SHA-256 integrity. It does not replace rerunning the manuscript with its archived inputs. After source edits, run `make manifest` before validation.

## Citation, license and acknowledgements

Use [CITATION.cff](CITATION.cff) for software attribution and cite the manuscript when its publication details are available. Code is released under the [MIT license](LICENSE); external datasets retain their own licenses.

The presentation of the workflow, installation and usage sections was inspired by [MetaGenome-MAG-Analysis](https://github.com/SongzeCHEN/MetaGenome-MAG-Analysis). Its MAG assembly/binning methods were not imported into this study. [REA](https://github.com/morluto/rea) was used for static inspection of the reference source archive; [redesign provenance](docs/provenance/REPOSITORY_REDESIGN.md) records its scope.
