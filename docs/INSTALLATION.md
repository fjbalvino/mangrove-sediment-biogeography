# Installation and environments

The workflows use Bash, Python and R. Work from a clone or a verified code release. Paths are supplied through `config/project.env` or command-line arguments; no institution-specific directories are required.

## Python

Python 3.12 is used by repository CI. Create a virtual environment, then install the functional dependencies:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r config/functional_requirements.txt
```

For acquisition and environmental reconstruction, also install:

```bash
python -m pip install -r config/python_requirements.txt -c config/functional_requirements.txt
```

The pinned functional requirements record the supplied analytical environment. Environmental requirements are compatibility ranges, not a complete historical lockfile. Geospatial Python packages may require system GDAL/PROJ libraries if binary wheels are unavailable.

## R

[config/R_packages.tsv](../config/R_packages.tsv) inventories analysis packages. Install CRAN dependencies and the two Bioconductor packages into the R environment you intend to use:

```r
packages <- read.delim("config/R_packages.tsv")$package
bioc <- c("phyloseq", "KEGGREST")
install.packages(setdiff(packages, bioc), repos = "https://cloud.r-project.org")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
BiocManager::install(bioc, ask = FALSE, update = FALSE)
```

This installs available compatible versions; it is not a substitute for the archived analytical environment. The regional analysis records its validated versions in [functional_regional_R_packages.tsv](../config/functional_regional_R_packages.tsv). The network driver deliberately requires **R 4.3.3, igraph 2.2.1 and data.table 1.17.8** and checks them before inference. Use an environment with those versions to reproduce the network. Do not relax its version or checksum checks to make a different environment appear identical. Compressed R input also requires `R.utils`.

For lightweight validation, base R runs syntax and taxonomic invariant checks; the regional statistical checks additionally need `vegan`. Plotting packages and real biological inputs are unnecessary for these tests. CI tests source contracts; it does not reproduce the complete analysis.

## Raw-read tools and external references

| Block | Additional dependencies | External assets |
| --- | --- | --- |
| 01: acquisition/QC/classification | FastQC, fastp, Kraken2; SRA Toolkit if the ENA download route is unavailable | Raw-read accessions, curated metadata, a versioned Kraken2 database and matching NCBI taxdump |
| 02: covariates | Geospatial Python packages; access to documented GBIF and SoilGrids services | Exact Human Footprint, WorldClim and MEOW assets, cached responses |
| 06: iCAMP | iCAMP (reference analysis: 1.8.6) and its R dependencies | Independent count object and GTDB R226 bacterial genus tree |
| 07: upstream functional annotation | Manuscript methods: MEGAHIT 1.2.9, Prodigal 2.6.3, MMseqs2 18-8cc5c, eggNOG-mapper 2.1.13 / database 5.0.2 | Original validated scripts and reference snapshots still need recovery |

The manuscript reports MEGAHIT meta-sensitive assembly and MMseqs2 easy-linclust at 90% identity and 80% coverage. This description does not constitute a validated executable reconstruction of the missing block.

Resource needs depend on input size and reference databases. Set worker/thread counts in `config/project.env` to fit your machine; these settings do not provision a server. The code uses no remote execution endpoint.

## Configure and check

```bash
cp config/project.env.example config/project.env
# Edit the copied file, especially MAGENTA_DATA_DIR and output paths.
make validate
```

The shell wrappers resolve their source directory independently of the current working directory. When invoking an individual R/Python script, use its documented flags or source the configuration first. A custom `--config` file is Bash code and must be a trusted local file.
