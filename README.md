# Global biogeography of mangrove microbiomes

Reproducible code for taxonomic diversity, composition, core-associated taxa,
community assembly and functional biogeography in the MAGENTA project.

## Reference data and analyses

The taxonomic reference contains **341 metagenomes**, each with **18,145 retained
reads**, and one accession manifest: `config/samples/taxonomic_accessions.txt`.
The metadata and count object are joined by `run_accession`.

| Analysis | Analytical units |
| --- | --- |
| Hill diversity q0, q1 and q2 | 341 profiles across five marine realms |
| Primary diversity comparisons | 120 sites in three realms with at least five sites |
| Environmental Shannon GAM | 274 profiles from 33 localities; 14 environmental predictors |
| Ordination and geographic distance decay | 323 profiles from 122 sampling sites |
| Compositional variation partitioning | 99 sites with complete environmental covariates |
| Core-associated taxa | 339 profiles from 47 realm-locality units |
| Functional atlas and taxonomic–functional coupling | 191 matched profiles; 17,482 KOs |
| Functional geographic analyses | 184 profiles from 57 sites; primary realm tests use 56 sites |
| iCAMP, separate completed analysis | 400 profiles, 1,892 genera, 24,692 within-realm pairs; 999 randomizations |

The eight profiles from Mud Island represent one sampling site. The Temperate
Northern Atlantic and Tropical Eastern Pacific each have one site and are shown
descriptively, without contributing to primary site-level realm comparisons.
Two profiles from Fujian have supported realm assignments but unresolved exact
localities; they remain in the reference and are omitted only where locality or
coordinates are required.

Sampling material varies across the retained reference: 242 profiles have
record-level support for mangrove sediment or soil, 46 represent enriched,
mixed, cultured or plant-tissue material, and 53 require further source
clarification. These annotations are retained in the metadata. Results describe
this specified reference; they do not establish that every profile represents
an unmanipulated natural sediment community.

## Reproduce the taxonomic analyses

Place the files specified in [DATA_AVAILABILITY.md](DATA_AVAILABILITY.md) under
`data/` and review `config/taxonomic.env.example`. Preview the analysis sequence:

```bash
bash workflow/03_RUN_TAXONOMIC_ANALYSES.sh --config config/taxonomic.env.example --core
```

Add `--execute` to run it. Each analysis writes to its current output directory;
filenames have no version suffixes. iCAMP retains its separately completed result
and is not rerun by this workflow. Its input, metadata and estimates are kept
separate from the 341-profile reference.

[The analysis guide](docs/TAXONOMIC_ANALYSES.md) defines the statistical units,
models, sample handling and interpretation. Current numerical summaries are in
`results/`; full tables, fitted objects and figures accompany the data archive.

## Reproduce the functional atlas

The functional atlas uses the 191 taxonomic profiles represented in the filtered
protein-annotation matrix. Counts quantify assembled proteins assigned to KOs.
The primary TA–CIP comparison uses one mean CLR profile per sampling site.

Install `config/functional_requirements.txt` and place the three compact inputs
listed in [DATA_AVAILABILITY.md](DATA_AVAILABILITY.md) under `data/functional/`.
Preview or run the workflow:

```bash
bash workflow/04_RUN_FUNCTIONAL_ANALYSES.sh --config config/functional.env.example
bash workflow/04_RUN_FUNCTIONAL_ANALYSES.sh --config config/functional.env.example --execute
```

[Functional methods and interpretation](docs/FUNCTIONAL_ANALYSES.md) describe the
class definitions, PERMANOVA and dispersion tests, count-adjusted analyses,
taxonomic–functional coupling, and geographic sampling sensitivity.
[Current results and Figure 4 caption](docs/FUNCTIONAL_RESULTS.md) use the same outputs.
Compact numerical results are in `results/functional/`.

The [functional network workflow](docs/FUNCTIONAL_NETWORK_ANALYSES.md) now uses
the same 191-profile cohort. [Figure 6 results and caption](docs/FUNCTIONAL_NETWORK_RESULTS.md)
report 1,060 associations, 327 connected KOs and five major modules.
[Supplementary results and methods](docs/FUNCTIONAL_NETWORK_SUPPLEMENTARY_RESULTS.md)
provide the complete enrichment tests and sensitivity analyses. Source statistics are in
`results/functional_network/`. The standalone driver records all commands and
performs independent output checks; images are generated locally.

## Workflow blocks

| Block | Objective |
| --- | --- |
| 01 | Metadata, sequence QC, Kraken2, taxonomic counts and reference selection |
| 02 | Human Footprint, WorldClim, GBIF, SoilGrids and MEOW covariates |
| 03 | Hill diversity, environmental GAM, dominance and Figure 2 |
| 04 | NMDS, distance decay, regional PERMANOVA and dbRDA/VPA |
| 05 | Core membership, realm restriction, core composition and Figure 3 |
| 06 | Independent genus-level iCAMP analysis and sampling-structure assessment |
| 07 | Assembly-to-KO annotation workflow; pending completion |
| 08 | Functional atlas, regional tests, taxonomic–functional coupling and Figure 4 |
| 09–10 | Separate regional functional and KO co-variation network analyses |

Functional analyses retain their separately documented cohorts. Repository checks
validate syntax, file integrity and statistical implementation invariants. They
do not replace execution of the analyses with their input data.

## Availability and citation

A public archive DOI has not yet been assigned. Input checksums and requirements
are in [DATA_AVAILABILITY.md](DATA_AVAILABILITY.md). Code is MIT licensed;
third-party data retain their source attribution and terms. Citation metadata
are in `CITATION.cff`.
