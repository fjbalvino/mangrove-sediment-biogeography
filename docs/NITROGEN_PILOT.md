# Nitrogen gene biogeography pilot

Executed 2 November 2026 (Mexico City date). The analysis retains all 191 accepted
mangrove-sediment profiles and joins metadata exclusively by `run_accession`.
Neither metadata source is edited. The 401-row master is not a functional matrix.

## Current result and limits

The updated functional metadata contains coordinates for all 191 profiles.
Recomputing rounded coordinate groups in analysis outputs gives 63 coordinate
groups and 85 study-coordinate groups across 44 study labels. Seven stored group
labels were missing even after coordinate recovery. `12A --coordinate-grouping
current` refreshes those labels locally and records the differences without
changing metadata. These groups are cartographic summaries, not verified
independent biological replicates.

BIO1 is available for 185 profiles; SOC and pH for 179. The three covariates jointly
support 81 study-coordinate groups from 41 studies. Across all profiles, BIO1
varies within eight studies and SOC/pH within nine. Support is lower for some
positive-only KO analyses. Adding realm to all study fixed effects still adds no
rank. Two BioSample identifiers each have three runs; the input records remain
separate. The relationship among their libraries still needs resolution before
run-level biological inference.

The frozen nitrogen rules consult 46 KOs. Thirty-four have matrix rows and twelve
do not. This is the nitrogen pilot panel, not a claim of exhaustive coverage of
all nitrogen genes or architectures. Missing rows remain unavailable, not zero.
The reaction catalogue separates 16 graphical cards, including substeps of the
same process; these are not 16 statistically independent processes.

## Measurement and descriptive summaries

For each selected profile, relative representation is `KO assignments / total
assignments across the complete input KO matrix * 1e6`. The denominator is fixed
before selecting nitrogen genes. Counts represent predicted-protein records
assigned to KOs, not mapped reads, TPM, copies per genome, enzyme concentrations,
or activity. Normalizing these records does not eliminate assembly or annotation
bias. Genes and required subunits are never summed as pathway abundance.

The per-run tables retain all 191 profiles. For a descriptive study-coordinate
summary, first average runs sharing a BioSample label, then average those labels
equally. If BioSample is unavailable, retain a distinct run-specific label. This
avoids treating the six known repeated-BioSample runs as six equally weighted
biological samples in the group view, without claiming that the label resolves
the sampling design. Never use this operation as a substitute for mapping reads
or combining count libraries. The geographic point layer shows the median of
study-coordinate summaries sharing a coordinate. The overall cycle displays the
median across 85 study-coordinate summaries. Neither is a reconstructed organism.

All views use the same `log10(1 + value)` color scale. Observed zeros are white;
unavailable KO rows are grey. The viewer includes group summaries and all
individual profiles. Within-group detection fractions at thresholds one and two
are retained separately in result tables, without claiming calibrated detection
probabilities. Study-point sample counts remain visible.

## Reaction interpretation

The catalogue links each displayed gene to its KO, reference transformation,
source and limitation. NifHDK supports nitrogenase components, NapAB a nitrate to
nitrite step, NirK/S alternative NO-forming nitrite reductases, NorBC a cNOR
architecture, NosZ N2O reduction, and NrfA nitrite ammonification. Community
co-detection does not assign a complete pathway or complex to an organism.

The historical Amo/Pmo annotation cannot be reinterpreted as current AmoA solely
because KEGG now uses K28504. Nar/Nxr direction remains unresolved. NirBD does not
distinguish assimilation from dissimilation. Hao is a candidate family, not a
sufficient nitrification marker. HzsABC/Hdh have no rows, so anammox is not
evaluable. Alternative nitrogenase, assimilation, urea and cyanate components
remain in the complete panel even when not shown in the compact cycle diagram.
Arrows are reference biochemical transformations, never measured net flux.

## Environmental diagnostic models

`12B` calculates exploratory within-study OLS slopes for two descriptive group
outcomes: log2 of positive group mean representation, and the group mean of
observed detection indicators. Each model includes study fixed effects and
standardized BIO1, SOC and pH jointly. The scale is the standard deviation across
all complete study-coordinate groups, frozen before per-KO positive filtering.
No pseudocount is used in positive-only log2 outcomes. Linear detection-fraction
slopes are diagnostic, not a fitted binomial occupancy model.

The estimates describe within-study differences among cartographic summaries.
They do not implement the proposed MaAsLin 3 run-level two-part model and are not
substitutes for its primary inference. There are no p-values, q-values, confidence
intervals or discovery claims. Constant-within-study detection outcomes and rank
failures are explicitly retained as non-estimable diagnostic states.

For each estimable KO, omit whole studies that have any within-study predictor
variation. Removing studies without contrasts cannot count as stability
evidence. The recorded minimum/maximum is a sensitivity range, not a confidence
interval. The original coefficient remains the reported coefficient; no model
is selected based on sign. A separate sensitivity adds log2 annotation yield and
is not treated as a definitive correction for the response denominator.

Among the 34 positive-representation KO models, at least one study omission
reverses the slope sign for 22 temperature associations, 13 SOC associations and
19 pH associations. These are diagnostic counts, not independent hypothesis
tests. This instability makes premature ecological attribution unsafe. There
are 186 estimated diagnostic coefficients in total (102 positive representation,
84 detection-fraction), plus explicit unavailable/constant/rank-status records.

SOC and pH are converted only in the derived analysis objects: stored SOC / 10
to g/kg and stored pH / 10. BIO1 remains degrees C. Environmental reconciliation
is documented in [the recovery record](ENVIRONMENTAL_METADATA_RECOVERY_2026-11.md):
WorldClim 2.1, 1970–2000, 30 arcsec; SoilGrids 2.0 mean, 0–5 cm. These are climatic
and modelled soil proxies, not sediment temperature at collection or bioavailable
carbon. Historical extraction-anchor discrepancies and approximate recovered
coordinates remain sensitivity issues. Human Footprint conflicts and GBIF
incompleteness are outside the three-variable pilot model.

## Reproduce

Use the existing external functional input bundle. Keep data and generated
figures outside Git. Python dependencies are numpy, pandas and matplotlib. The
plotter also requires a local Natural Earth 1:110 million countries GeoJSON. The
executed file came from `nvkelso/natural-earth-vector`,
`geojson/ne_110m_admin_0_countries.geojson`; its exact digest is in the evidence
package. Natural Earth data is public domain.

```bash
python scripts/12_gene_environment/12A_AUDIT_GENE_ENVIRONMENT_DESIGN.py \
  --metadata /path/metadata_functional.csv \
  --matrix /path/functional_KO_counts_original.tsv.gz \
  --samples config/samples/functional_accessions.txt \
  --unit predicted_protein_records --coordinate-grouping current \
  --outdir /path/new_audit

python scripts/12_gene_environment/12B_BUILD_NITROGEN_PILOT.py \
  --metadata /path/metadata_functional.csv \
  --matrix /path/functional_KO_counts_original.tsv.gz \
  --samples config/samples/functional_accessions.txt \
  --rules config/biogeochemical_cycles/reglas_operativas.tsv \
  --reactions config/gene_environment/nitrogen_reactions.json \
  --outdir /path/new_nitrogen_results

python scripts/12_gene_environment/12C_PLOT_NITROGEN_PILOT.py \
  --results /path/new_nitrogen_results \
  --reactions config/gene_environment/nitrogen_reactions.json \
  --world /path/ne_110m_admin_0_countries.geojson \
  --outdir /path/nitrogen_figures

python -m unittest discover -s tests -p 'test_gene_environment_*.py'
```

The HTML template must remain next to `12C`. The generated viewer embeds its own
data and geography and requires no remote map requests. Its interactive behavior
was checked across all 277 selection options, 46 map genes and 63 coordinate
clicks using a DOM integration test. Static PNG/PDF source plots were visually
inspected. A full browser screenshot check was not completed because the browser
download failed; do not claim browser-specific layout certification.

## Next inferential gate

Resolve biological/technical relationships and the scale of independent sampling
units for each contrast. Quantify study dependence and coordinate/proxy
uncertainty. Pin and validate the primary engine on synthetic nulls matching the
actual study/site structure and zeros, including calibration of error rates and
cluster dependence, before computing primary significance. The current unit
tests validate arithmetic and data contracts only, not model error rates or
biological validity. A low-dimensional nitrogen comparison can then precede the
broader 418-KO panel. Repeat unchanged measurement-compatible steps on Dani's
matrix; if quantification units change, revise the measurement model explicitly.

Sources revisited for this pilot include [KEGG nitrogen metabolism](https://www.kegg.jp/entry/map00910),
[M00531](https://www.kegg.jp/entry/M00531),
[K28504](https://www.kegg.jp/entry/K28504),
[K10944](https://www.kegg.jp/entry/K10944),
[K00370](https://www.kegg.jp/entry/K00370),
the K22896–K22899 entries, and [Natural Earth terms](https://www.naturalearthdata.com/about/terms-of-use/).
The frozen catalogue retains its original references. The
[earlier protocol](GENE_ENVIRONMENT_BIOGEOGRAPHY.md) contains historical audit
counts; this page supersedes its coordinate and environmental availability
counts, not its stated methodological limitations.
