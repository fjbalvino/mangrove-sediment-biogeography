# How to run the manuscript analyses

Use the [installation guide](../docs/INSTALLATION.md) and [input inventory](../DATA_AVAILABILITY.md) first. Copy `config/project.env.example` to `config/project.env` and edit it. The default local configuration is loaded automatically; `--config FILE` or `MAGENTA_CONFIG` selects another file.

Shell entry points **preview by default**. Add `--execute` only when ready to run. Previews do not read biological inputs, start downloads, create output directories or verify that the scientific inputs are correct.

## Choose a starting point

- With the archived 341-profile phyloseq and metadata, start at workflow 03.
- To rebuild from raw reads, run 01 and 02 before 03, using the exact upstream metadata and reference assets.
- With the archived functional matrices, start at 04; run 05 and 06 as separate branches.
- iCAMP uses its own count object/tree and is independent of workflow 03. Never pass the 341-profile object as a substitute.

## 01–03: taxonomy and environmental context

```bash
bash workflow/01_RUN_TAXONOMIC_PREPROCESSING.sh --config config/project.env
bash workflow/02_RUN_ENVIRONMENTAL_COVARIATES.sh --config config/project.env
bash workflow/03_RUN_TAXONOMIC_ANALYSES.sh --config config/project.env --core
```

Workflow 01 supports `--local-fastq-map FILE` for local reads. Workflow 02 accepts `--base-metadata FILE`, and refuses to overwrite its final metadata destination. Its example coastal anchor limit is the scripts' 2,000 m default; reproducing environmental extraction also requires the exact source assets and cached service responses.

Workflow 03 covers Hill diversity, the Shannon GAM, regional composition, geographic decay, variation partitioning and optional taxonomic core. It uses the external iCAMP process summary only for the legacy figure assembler; it does not rerun iCAMP.

For an independent iCAMP run, review the source script's arguments and its separate cohort before executing:

```bash
source config/project.env
"$R_BIN" scripts/06_community_assembly/06A_ICAMP_GENUS_COMMUNITY_ASSEMBLY.R --help
```

[Taxonomic methods](../docs/TAXONOMIC_ANALYSES.md) explain sample filtering, units, inferential guards and the separate assembly analysis.

## 04: functional atlas

Required under `FUNCTIONAL_INPUT_DIR`: the original KO count matrix, corrected metadata and matched taxonomic counts.

```bash
bash workflow/04_RUN_FUNCTIONAL_ANALYSES.sh --config config/project.env
# Add --execute to run analyses and render Figure 4 locally.
```

The downstream network/regional workflows also require the frozen persistence classification and annotation snapshots. Do not substitute an annotation database downloaded today for a historical snapshot. See [functional methods](../docs/FUNCTIONAL_ANALYSES.md).

## 05: functional co-abundance network — manuscript Figure 5

The Python driver runs when invoked normally. Use `--dry-run` to inspect its plan without execution:

```bash
source config/project.env
"$PYTHON_BIN" workflow/05_RUN_FUNCTIONAL_NETWORK.py \
  --counts "$FUNCTIONAL_INPUT_DIR/functional_KO_counts_original.tsv.gz" \
  --metadata "$FUNCTIONAL_INPUT_DIR/metadata_functional.csv" \
  --atlas-classes "$FUNCTIONAL_INPUT_DIR/KO_persistence_classes.tsv" \
  --annotations "$FUNCTIONAL_INPUT_DIR/network_annotations" \
  --accessions "$FUNCTIONAL_SAMPLE_MANIFEST" \
  --outdir "$BIOGEOGRAPHY_RESULTS/functional_network_full" \
  --rscript "$R_BIN" --workers 4 --stop-after G --dry-run
```

Remove `--dry-run` to execute. `network_annotations/` must contain `137_KO_function_dictionary_used.tsv` and `144_KO_to_KEGG_MODULE_long.tsv`. The driver verifies the source hashes and pinned R environment before inference. Stopping at F omits study-omission sensitivity and the main figure; the complete manuscript workflow runs through G.

The [network methods](../docs/FUNCTIONAL_NETWORK_ANALYSES.md) describe restart controls, frozen annotation coverage, subsampling and independent checks. Changing the accession/checksum file does not by itself adapt every manuscript-specific validation to a new study.

## 06: regional functional differentiation — manuscript Figure 6

Four inputs are required: the original KO matrix, corrected metadata, `KO_persistence_classes.tsv` and `KO_function_dictionary.tsv`. The Natural Earth basemap is an **external** file, configured with `WORLD_MAP` or `--world-map`.

```bash
bash workflow/06_RUN_FUNCTIONAL_REGIONAL.sh --config config/project.env
# Add --execute after checking all inputs and the basemap.
```

Use `--inputs DIR` and `--outdir DIR` to override configured locations. `--complete-panels` reuses and verifies existing balanced results before completing process/map analyses. Numerical verification runs before rendering. The [regional methods](../docs/FUNCTIONAL_REGIONAL_ANALYSES.md) describe the 3 × 6 and 3 × 4 locality designs. Historical output filenames beginning `Figure_5` refer to the current Figure 6.

## Outputs, validation and release

Generated matrices, figures, fitted objects and reports stay in local output directories. None is uploaded by these workflows.

```bash
make validate
bash workflow/99_BUILD_EXPORT_BUNDLE.sh --output ../manuscript_code.zip
bash workflow/99_BUILD_EXPORT_BUNDLE.sh --output ../manuscript_code.zip --execute
```

The exporter packages only verified source files, even when local data/results are present. It does not include Git history. A recovery backup is a separate artifact.

## Adapting the code

Change local paths and executables in the configuration. For another dataset, also supply appropriate accession manifests, metadata columns, biological definitions and input checksums. Review expected cohort/class counts and interpretation in the individual scripts; the distributed wrappers are the manuscript's reproducibility recipes. Do not bypass their scientific guards and report a new analysis as a reproduction of MAGENTA.
