# Re-running biogeochemical gene biogeography

This stage reproduces the eight-cycle atlas, operational rule evidence, environmental design diagnostics and exploratory gene–environment slopes. It integrates nitrogen with carbon, methane/C1, sulfur, phosphorus, iron acquisition, hydrogen/electron transfer, and respiration/phototrophy. The last two panels provide metabolic context; iron redox remains unresolved.

The source derives from version 7 of the reviewed 191-profile package. That reference contains 97 rules, 65 curated processes, 69 displayed reaction groups and 368 unique queried KOs; 328 have matrix rows. These are historical reference counts, **not imposed dimensions for the next matrix**. A new run derives sample, coordinate, study, available-KO and calibration-family sizes from its inputs. The rule catalogue remains fixed until an explicit scientific curation revision.

## Scope and measurement

- Input measurements are **predicted-protein annotation records** associated with KOs, not mapped reads, TPM, coverage or absolute gene abundance. Integer values alone do not establish this provenance.
- Normalize each KO against the sum across **all rows of the complete input KO matrix**, before selecting a panel, then multiply by one million.
- A missing KO row is unknown; a zero in an available row is an observed non-detection. Neither establishes biological absence.
- Runs are averaged within BioSample labels, then labels equally within study–coordinate groups. Missing BioSample IDs use separate run labels. Study labels and coordinate bins do not certify biological independence.
- Operational AND/OR rules describe community co-detection. They establish neither a complete pathway in an organism nor activity, metabolic direction or flux. Subunits are not summed as pathway abundance.
- Observed-gene associations remain exploratory. No observed P/q values or confirmatory intervals are generated, even if a future calibration passes; interpretation still requires scientific review.

## Inputs and dependencies

Use Python 3.12 and the pinned versions in `config/gene_biogeography/requirements.txt`:

```bash
python3.12 -m venv .venv-gene
source .venv-gene/bin/activate
python -m pip install -r config/gene_biogeography/requirements.txt
```

The calibration implementation was developed with R 4.3.3 and `clubSandwich` 0.5.10. Verify these versions locally; the calibration check compares the Python implementation with the installed official R implementation. Node with `jsdom` 26.1.0 is optional for viewer interaction tests. No installation downloads are performed by the analysis driver.

| Input | Contract |
| --- | --- |
| KO matrix | TSV or TSV.gz; first column unique `K#####`; remaining columns unique run IDs; finite nonnegative integer record counts; positive total per profile |
| Metadata | CSV; one unique `run_accession` per row; a larger universal metadata table is allowed |
| Selected cohort | Explicit accession text file, or explicit `--all-matrix-samples`; no automatic restriction to the old 191 |
| Required metadata | `run_accession`, `study_id`, `sample_accession`, `lat`, `lon`, `Country`, `Locality`, `REALM`, `bio1`, `soc_mean`, `phh2o_mean` |
| Environmental units | BIO1 in degrees C. Default SOC/pH storage is SoilGrids encoded values ×10; `--soilgrids-storage physical` accepts SOC in g/kg and pH units |
| Visual assets | External `atlas/assets` directory containing `world_110m.geojson`, `natural_earth_relief.jpg`, `mangroves_USGS.jpg`, `nitrogen_reference_template.svg`; use the accepted package's assets |

Rows are matched **only by `run_accession`**. Duplicate or unmatched IDs fail explicitly. Extra metadata rows do not enter the cohort. Original metadata and matrix files are never edited. Unit conversions and grouping fields exist only in derived run files. Environmental missing values remain missing and are handled by complete-case diagnostics; they are not imputed.

Every selected profile needs valid coordinates for this geographic pipeline. If some coordinates remain unresolved, use the existing `workflow/07_RUN_BIOGEOCHEMICAL_CYCLES.py` for coordinate-free rule evidence while completing them; do not invent coordinates or silently drop profiles. No habitat decisions are made automatically: the approved cohort is supplied by the researcher.

## 1. Reproduce the reference or prepare a new cohort

Run from the repository root. Commands preview by default; add `--execute` after reviewing paths. `--run-dir` must be new/empty for preparation, outside the repository or under its ignored `results/` directory.

Reference example, using the previously delivered archive after extraction:

```bash
python workflow/08_RUN_GENE_BIOGEOGRAPHY.py --stage prepare \
  --matrix /path/to/biogeochemical_atlas_191/inputs/functional_KO_counts_original.tsv.gz \
  --metadata /path/to/biogeochemical_atlas_191/inputs/metadata_functional.csv \
  --samples config/samples/functional_accessions.txt \
  --unit predicted_protein_records --expected-profiles 191 \
  --run-dir /path/to/runs/reference-191 --execute
```

For Dani's complete matrix, first confirm measurement and annotation compatibility. Then use the approved new accession list:

```bash
python workflow/08_RUN_GENE_BIOGEOGRAPHY.py --stage prepare \
  --matrix /path/to/KEGG_ABUNDANCE_MATRIX.tsv \
  --metadata /path/to/current_universal_metadata.csv \
  --samples /path/to/approved_functional_accessions.txt \
  --unit predicted_protein_records \
  --run-dir /path/to/runs/functional-complete --execute
```

Alternatively replace `--samples ...` with `--all-matrix-samples` when every matrix column is explicitly approved. Do not retain `--expected-profiles 191` for the expanded dataset. New sample counts are determined from the selected IDs, not from metadata row count or file names. If the measurement changes to mapped counts/coverage/TPM, stop and revise the quantification/statistical contract instead of relabelling it.

Preparation stages: validate inputs → stage immutable analysis-local copies → nitrogen summaries/slopes → design support → remaining-cycle summaries/slopes → curate labels → integrate shared KOs by ID → verify counts, full-matrix denominators, missingness, hierarchical weights and rule states. Shared KOs must agree before integration.

## 2. Generate all eight viewers and local figures

```bash
python workflow/08_RUN_GENE_BIOGEOGRAPHY.py --stage atlas \
  --run-dir /path/to/runs/functional-complete \
  --assets /path/to/biogeochemical_atlas_191/atlas/assets --execute
```

Open the run's `atlas/index.html`. Every panel includes map zoom and overlap selection, a biochemical reference diagram, an individual-gene heatmap, environmental sensitivity, rules and methods. Plasma is shared across all panels within a release; the upper scale is recalculated for that release and recorded. Compare colour scales explicitly between releases. Heatmap columns decrease in signed longitude, east to west. Figures and derived data stay local.

The map retains the corrected Leaflet layer order: offline overview below terrain; sample markers above terrain. OpenTopoMap tiles load on demand, never in bulk. Terrain availability depends on internet/provider access and does not improve coordinate precision. Asset provenance and attribution are visible in the viewers. The landing-page manuscript links refer to the historical supplements, not automatically to the new cohort.

Local verification:

```bash
python /path/to/runs/functional-complete/verify_outputs.py --viewers
# Optional: jsdom must be installed in the Node environment.
ATLAS_TEST_ROOT=/path/to/runs/functional-complete/atlas \
  node /path/to/runs/functional-complete/atlas/tests/test_viewers.cjs
```

The Node checks cover selection, marker colour changes, thresholds, map zoom, overlap lists, CSS sizing and layer order. They do not replace visual inspection in a real browser. Local supplementary panels are generated under `atlas/supplement/`; figure numbers are assigned editorially rather than carried into every future dataset.

## 3. Prepare and verify a fresh calibration

```bash
python workflow/08_RUN_GENE_BIOGEOGRAPHY.py --stage calibration-plan \
  --run-dir /path/to/runs/functional-complete --execute
python workflow/08_RUN_GENE_BIOGEOGRAPHY.py --stage calibration-check \
  --run-dir /path/to/runs/functional-complete --execute
```

Two families remain separate: available nitrogen KOs and available remaining-panel KOs, three covariates and two response components. Overlapping KOs occur in both historical families; presentation integration does not invent a pooled inferential family. Plans report the current family size and retain available members of the prespecified biochemical signal-marker lists. Unsupported all-zero positive margins or undefined covariate/recovery scales stop calibration explicitly while descriptive results remain available.

The dominant-heterogeneity scenario freezes the study contributing the largest within-study temperature sum of squares before responses are simulated; it is not permanently tied to a historical BioProject. With the reference data this retains the original study. Simulation mechanisms, signal magnitudes, seeds, 1,000 null plus 500 mixed datasets per scenario, multiplicity rules and acceptance criteria are retained. Five scenarios therefore generate 7,500 datasets per family. No criterion is relaxed to obtain a favourable result.

## 4. Execute calibration explicitly

```bash
python workflow/08_RUN_GENE_BIOGEOGRAPHY.py --stage calibrate \
  --run-dir /path/to/runs/functional-complete --workers 2 --execute
```

This is computationally expensive and is not part of the ordinary atlas build. Run on Linux/Unix with sufficient resources. Resuming requires unchanged input, source and plan fingerprints. Never reuse the old 191-profile checkpoints for the new matrix. Failed or non-estimable cells remain reported; invalid tests retain P=1 within fixed simulation families. The historical reference did not pass all combined validity, coverage, type-I/FDR criteria under either candidate uncertainty procedure. This result must not be copied into a future run without recalibration.

Rebuild the atlas after calibration to refresh the displayed status. Model summaries are not automatically translated into observed-gene significance claims.

## Reproducibility and relationship to the rest of the manuscript

`run_manifest.json` records source input paths/hashes, measurement, units, cohort size, staged source/configuration hashes and Python version. Each stage writes a log and pass/fail record. Model outputs record NumPy/pandas versions; reference package hashes are in `config/gene_biogeography/source_provenance.json`. Preserve the run directory and original inputs together in the project's data storage.

The existing functional backbone and regional workflows remain under `workflow/04_RUN_FUNCTIONAL_ANALYSES.sh` and `workflow/06_RUN_FUNCTIONAL_REGIONAL.sh`; see their input requirements and explicit cohort settings before rerunning. The old network code is retained for provenance but **network redesign is pending** and is not called by this new driver. Upstream assembly/protein prediction/clustering/annotation is not implemented by the gene-matrix workflow. The historical 191 profiles remain the accepted reference cohort.

For public hosting, see [Atlas publication](ATLAS_PUBLICATION.md). Only source code, small curated rule/label configuration, tests and documentation belong in this update. Input matrices, metadata, generated viewers containing data, images, figures, simulations and archives remain outside Git.
