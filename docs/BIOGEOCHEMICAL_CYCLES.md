# Biogeochemical-cycle evidence from the functional KO matrix

This section evaluates community-level genetic evidence for carbon, methane,
nitrogen, sulfur and phosphorus transformations, with iron and energy-related
context. It starts from the **191 accepted mangrove-sediment metagenomes** in
`config/samples/functional_accessions.txt`. It does not reassess habitat eligibility.

The current implementation contains **97 operational rules** covering **56 of 65
curated categories**, plus **32 KEGG reference modules**. Nine categories explicitly
require evidence unavailable from this matrix. This is a defined analysis panel,
not an exhaustive inventory of all genes involved in biogeochemical cycling.

## What the numbers mean

The input contains counts of predicted protein records assigned to KOs, recovered
from assemblies and the protein-clustering workflow. These are not mapped-read
counts, coverage, TPM, absolute gene abundance or expression. A rule describes
components recovered together in a community; it does not establish their presence
in one organism, pathway activity, reaction direction or biogeochemical flux.

For each KO, the evaluator compares its record count against two thresholds:
`>=1` (primary description) and `>=2` (sensitivity). The second threshold is an
exploratory record-count sensitivity, not an experimentally calibrated detection
limit or independent confirmation of an annotation.

| Logical result | Meaning |
| --- | --- |
| `TRUE` / matrix value `1` | The selected rule is satisfied at that threshold. |
| `FALSE` / matrix value `0` | Observed values make the selected rule fail; this is not biological absence. |
| `UNKNOWN` / matrix value `NA` | Missing KO rows prevent a decision for that rule. |

An absent KO row remains unknown. An existing row with zero is an observed zero in
this annotation table. Neither is silently replaced by the other. In three-valued
logic, a false required component makes an AND false even when another is unknown;
a true alternative makes an OR true even when another is unknown. Other unresolved
combinations remain unknown. Optional components do not gate rule satisfaction.

The parser implements the subset of [KEGG module syntax](https://www.kegg.jp/kegg/module.html)
used by the frozen definitions: spaces and `+` are AND, commas are OR, parentheses
group terms, and `-` marks optional components. Dependencies M00597 and M00598 are
expanded before evaluation. Unsupported tokens and cyclic dependencies fail.

## Decision record and current scope

| Decision | Current choice and consequence |
| --- | --- |
| D0: available data | An existing WGS-derived KO matrix; no sequence reannotation in this stage. |
| D1: objective | Describe recovered components and their detection sensitivity. Geographic hypothesis tests are a subsequent stage. |
| D2: sampling design | Preserve the selected accessions; audit study/realm structure, supplied site labels and repeated BioSample identifiers. Coordinates alone do not establish independent sampling units. |
| D3: measurement | Audit count units, IDs and annotation-yield dependence. Record totals are not a sequencing-depth calibration. |
| D4: analytical route | Apply versioned KO rules directly. Do not claim this is a METABOLIC, FeGenie or dbCAN run. |
| D5: functional evidence | Separate catalytic sets, ambiguous markers, contextual components and module definitions. Retain alternatives and limitations. |
| D6: statistics | Descriptive frequencies and a Spearman detection diagnostic; no regional P-values or network inference in this stage. |
| D7: robustness | Compare thresholds, missing-row semantics, sample-order invariance and independently calculated catalytic sets. |
| D8: delivery | External numerical outputs, input/configuration/code hashes, versions, execution log and PASS/FAIL record. |

**Viable:** panel evaluation and descriptive comparisons of recovered components.
**Conditional:** regional or environmental inference after selecting a defensible
estimand, independent units, study structure and measurement model.
**Not estimable from this matrix alone:** expression, flux, organism-specific
pathway completeness, taxonomic carriers and interaction mechanisms.

The design diagnostic reports the rank increment from adding realm to study fixed
effects. A zero increment means that particular conditional fixed-effect contrast
is not identifiable. It does not prohibit every between-study regional comparison:
replicated studies within regions can support an appropriately specified
between-study model, under explicit comparability assumptions. A region represented
by one study has especially limited generalizability. Adding random effects does
not recover information absent from the sampling design.

## Inputs and installation

Use Python 3.12 and the pinned dependencies in
[`functional_requirements.txt`](../config/functional_requirements.txt). No online
service or reference-database download is needed during evaluation.

| Argument | Contract |
| --- | --- |
| `--matrix` | External TSV or TSV.gz; first column contains unique `K#####` identifiers, subsequent columns unique sample IDs. Selected values must be finite, nonnegative integers with the declared provenance. Extra sample columns are ignored. |
| `--metadata` | External CSV with unique `run_accession`; required columns: `study_id`, `sampling_site_id`, `REALM`, `lat`, `lon`. Study and realm must be populated for selected profiles. Site/coordinates may be blank. IDs are read as strings. |
| `--samples` | Unique accessions, one per line; all must exist in both inputs. File order is preserved. Defaults to the existing 191-accession manifest. |
| `--unit` | Explicitly declare `predicted_protein_records`. Integer values do not establish this unit by themselves. Normalized or mapped-read inputs require a separately reviewed measurement contract. |
| `--ko-labels` | Optional external TSV with unique `KO` and `name` columns; display names only, with no effect on rule evaluation. |
| `--outdir` | A new, empty local output directory. Within this repository it must be below `results/`. |

Optional metadata fields include `sample_accession`, `study_id_basis`,
`sampling_site_basis`, `Country` and `Locality`. A repeated BioSample is reported
for review; profiles are not automatically merged or deleted. Record the origin of
study/site labels in the supplied metadata. Missing coordinates are reported;
the workflow neither geocodes nor imputes them.

## Run and reproduce

From the repository root, preview without creating files:

```bash
python workflow/07_RUN_BIOGEOCHEMICAL_CYCLES.py \
  --matrix /path/to/functional_KO_counts_original.tsv.gz \
  --metadata /path/to/metadata_functional.csv \
  --unit predicted_protein_records \
  --outdir results/biogeochemical_cycles_current
```

Add `--execute` to evaluate and independently verify the output tables. Add
`--plots` as well to generate four descriptive figures **locally**, after numerical
verification. The workflow does not upload anything. It refuses an occupied output
directory so that input releases do not overwrite each other. Paths can contain
spaces; commands are executed as argument lists, without shell interpolation.

The individual scripts accept `--help`:

| Script | Responsibility |
| --- | --- |
| [11A evaluator](../scripts/11_biogeochemical_cycles/11A_EVALUATE_BIOGEOCHEMICAL_RULES.py) | Input checks, ternary logic, rule/block tables, design and detection diagnostics. |
| [11C verifier](../scripts/11_biogeochemical_cycles/11C_VERIFY_BIOGEOCHEMICAL_OUTPUTS.py) | Reconcile long tables with matrices and summaries, retain unknown values and check threshold monotonicity. |
| [11B plotter](../scripts/11_biogeochemical_cycles/11B_PLOT_BIOGEOCHEMICAL_EVIDENCE.py) | Render the threshold comparison, detection diagnostic and three category atlases. |

Run Python without `-O`; verification assertions must remain enabled. Successful
execution is a computational check, not validation of a biological mechanism.

## Local outputs

| File/group under `tables/` | Interpretation |
| --- | --- |
| `matriz_97_reglas_ge1.tsv`, `...ge2.tsv` | Rule-by-sample logical matrices: 1, 0 or NA. |
| `evidencia_por_regla_muestra.tsv.gz` | Every rule, sample and threshold; observed, subthreshold and missing-row KOs. Includes the 32 reference modules. |
| `evidencia_por_bloque_muestra.tsv.gz` | Logical evidence for each top-level expression block. |
| `resumen_reglas.tsv` | TRUE/FALSE/UNKNOWN frequencies with the full sample denominator. Unknown profiles are retained. |
| `matriz_65_reglas_representativas_ge*.tsv` | Display states for one representative rule per category; not an aggregate score over all alternatives. |
| `alcance_65_procesos.tsv` | The exact representative scope, alternatives count and category-specific limitation. |
| `inventario_KO_panel.tsv`, `conteos_KO_panel.tsv.gz` | Recoverability and external matrix values for the evaluated KO panel. |
| `muestras_diseno_y_deteccion.tsv` and design tables | Supplied study/site structure, repeated sample accessions and annotation-yield diagnostics. |
| `diagnostico_y_validacion.json` | Dimensions, numerical checks, versions and input/configuration/evaluator SHA-256 hashes. |
| `verificacion_salidas.json` | Independent table-consistency verification. |

The output root contains `execution.log` and `RUN_STATUS.json`, including command
arguments and output hashes. `figures/` is optional. The 65-row atlas has a separate
display state for unresolved categories; its colors do not replace the ternary
logical matrices. File suffixes `97` and `65` identify this frozen configuration;
sample count is taken from the selected manifest rather than hard-coded to 191.

## Interpretation safeguards and the next matrix

NifHDK is a catalytic-set rule; NifHDKENB is an additional conservative criterion,
not a universal definition of all nitrogen fixers. McrABG alone does not establish
methanogenesis rather than methane/alkane oxidation. Amo/Pmo and Nar/Nxr require
disambiguation. Sat–AprAB–DsrAB does not by itself identify sulfur reaction
direction. PhnJ and PhnGHIJKLM evaluate different scopes. Transporters and generic
cytochromes must not be relabeled as complete biogeochemical processes. These
limitations are retained per rule in the configuration; see the
[supporting literature](BIOGEOCHEMICAL_REFERENCES.md).

The number of satisfied rules includes overlapping sets and is used only as a
detection diagnostic, not a functional-diversity index. The plotted Spearman
coefficient describes dependence on total annotation yield; it is not a causal
effect, regional test or correction for detection bias.

For the updated matrix, retain the same frozen rules for a controlled comparison,
provide new input paths and a new output directory, and review the provenance/unit
contract. Supply a different accession manifest only when intentionally changing
the selected profiles. Preserve all eligible KO rows, including documented zeros;
record annotation/database versions and upstream gene/contig/cluster identifiers
outside this source repository. An upstream database update and a change to the
analytical rules are distinct changes and should be assessed separately.

Sequence-aware iron, CAZyme, hydrogenase and other specialist analyses require the
corresponding proteins, annotation confidence and genomic context. Replacing this
KO matrix alone does not resolve all nine currently blocked categories. Likewise,
a co-abundance edge cannot demonstrate cooperation or biochemical coupling.

## Verification and contribution boundary

The repository adaptation was rerun on the current 191 profiles. Both threshold
matrices, representative-category tables, summaries and long-form rule/block
evidence were compared with the preceding analysis. Synthetic tests additionally
exercise missing catalytic rows, duplicate headers, invalid counts, reordered
samples, preview behavior, output overwrite protection and deliberate corruption
of a result summary.

```bash
python -m unittest discover -s tests -p 'test_functional_cycles.py' -v
make manifest
make validate
```

Commit scripts, tests, documentation and small curated rule definitions only.
**Do not commit matrices, biological databases, raw KEGG records, numerical
outputs, plots, reports or archives.** Existing repository validation enforces
this boundary. Keep the broad literature catalogue and manuscript annex material
separate until explicitly incorporated into the manuscript workflow.
