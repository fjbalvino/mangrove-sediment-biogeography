# Marine salinity reference covariates

This extension adds two surface-water climatological predictors to the restart
metadata without selecting a cohort or rerunning manuscript analyses. Accessions,
row order, coordinates, every other field value and relative column order are
preserved. These are regional marine reference values, not measurements of
sediment porewater or conditions at the sampling date.

| Metadata field | Bio-ORACLE statistic | Definition |
|---|---|---|
| `marine_salinity_mean` | `so_mean` | Equal mean of the two decadal mean surface salinities |
| `marine_salinity_annual_range` | `so_range` | Equal mean of the two decadal averages of annual maximum-minus-minimum salinity |

Both statistics require both complete decades. Zero is valid; missing values are
`NA`, never zero. The annual range is neither the extreme range across twenty
years nor an uncertainty interval. Ten significant digits retain computational
precision without implying equivalent measurement accuracy. The provider reports
practical salinity as PSU.

## Source and period

- Bio-ORACLE v3.0, ERDDAP dataset `so_baseline_2000_2019_depthsurf`.
- Reference period 2000–2019, with decadal coordinates `2000-01-01` and `2010-01-01`
  for 2000–2009 and 2010–2019. These are not single-day observations. The website
  also describes present conditions as 2000–2020; this code follows the dataset's
  explicit 2000–2019 identifier and title.
- WGS84 grid, 0.05°; surface source depth 0–0.49 m. Interpolated grid resolution
  does not establish ecological representativeness within small lagoons/estuaries.
- [ERDDAP dataset](https://erddap.bio-oracle.org/erddap/griddap/so_baseline_2000_2019_depthsurf.html).
- [Official statistic definitions and surface depth](https://www.bio-oracle.org/documentation.php).
- Assis et al. (2024), [doi:10.1111/geb.13813](https://doi.org/10.1111/geb.13813).

## Update all restart rows

Use external paths and a separate output directory. The independent step preserves
all rows, including samples outside historical cohorts and missing coordinates.
Do not run 02A to update the entire restart universe: 02A deliberately selects the
historical taxonomic reference.

```bash
python scripts/02_environmental_covariates/02H_EXTRACT_BIOORACLE_SALINITY.py \
  --base "$METADATA_INPUT_DIR/Metadata_reinicio_401_reconciliada.csv" \
  --out "$METADATA_OUTPUT_DIR/Metadata_reinicio_401_reconciliada.csv" \
  --audit-dir "$METADATA_OUTPUT_DIR/salinity_audit" \
  --cache-dir "$ENVIRONMENTAL_CACHE/biooracle"
```

The extractor needs only Python's standard library. Network failures and malformed
or incomplete responses stop the run with `FAIL`; they do not become missing
salinity. Rerun verified cached requests with `--offline`; missing or changed cache
entries fail. Use `--replace-salinity` explicitly when the two fields already exist.
The original snapshot is never overwritten by the script.

The existing environmental orchestrator optionally runs 02H before 02G validation
when `BIOORACLE_ENABLE=1`. Default is off, preserving historical behavior and its
existing cohort selection. `BIOORACLE_REVIEW` accepts an external review CSV.
The validator reports salinity missingness separately; the fourteen-predictor GAM
is unchanged.

## Coastal assignment

The automatic extraction uses only the native containing grid cell. Exact grid
boundaries select north/east; endpoints clamp to the grid. A 7×7 patch is downloaded
per native location. Original coordinates are not rounded together: neighboring
points may share a grid cell but remain separate coordinate records and accessions.

Coordinate-review statuses containing `conflicto`, `por_verificar`, `por_resolver`
or `pendientes_validacion` withhold assignment even when the cell has data. Other
inherited/approximate coordinate warnings remain in the audit. This is not a fresh
validation of every original GPS record.

Nearby complete cells are exported as candidates only. An explicit review may
approve a cell after examining shoreline and waterbody correspondence, or hold a
native value after detecting a spatial concern. An unresolved lagoon connection,
opposite coast, offshore mangrove coordinate or locality conflict stays `NA`.
The 10 km review bound is operational, not a universal ecological validity limit;
it can be changed explicitly up to 20 km. Candidate patches do not wrap at the
dateline, conservatively limiting fallback coverage without affecting native cells.

External review fields:

```text
run_accession,lat,lon,decision,cell_lat,cell_lon,reason,evidence
```

`decision` is `approve` or `hold`. Original coordinate strings must match exactly.
Approved coordinates must identify grid centers with both complete decades within
the explicit distance bound. Every decision needs a reason and evidence; an
approval must explicitly resolve any inherited coordinate conflict. Apply the
review with `--review`. Keep review tables and map evidence with external data.
An image can support regional coastal correspondence, not validate local salinity,
tidal exchange or hydrological equivalence to the pixel.

## Audit outputs

- `sample_salinity_audit.csv`: one row per accession, original warnings, decision,
  native/selected cells, distances, evidence, four decadal and two final values.
- `candidate_cells.csv`: available candidate cells, not approved assignments.
  Coordinates without candidates remain represented in the sample audit.
- `salinity_manifest.json`: parameters, Python version, counts, input/output/review/
  script hashes and source URLs, retrieval timestamps and payload hashes.
- `source_metadata.json`, `salinity_dictionary.json`, `STATUS.json`.
- Keep the request cache and spatial review evidence for an offline reconstruction.

The original `covariate_review_status` describes legacy covariates and is preserved;
this separate audit documents the new fields. Shared coordinates/pixels do not
create independent environmental replicates. Future models should examine
native-only versus reviewed-coastal subsets, missingness and positional uncertainty
before interpreting salinity associations. No fitted model changes in this step.

## Recover remaining NA using cardinal searches (02J)

The requested follow-up uses the radial-search pattern in
`02E_EXTRACT_SOILGRIDS_TOPSOIL_0_30CM.py`: query the original point, increase the
radius, and accept the first complete result in the direction sequence. 02E uses
250, 500, 1000 and 2000 m, with eight configurable bearings by default. For this
marine layer, 02J uses exactly **N, E, S, W** and reuses
`metadata_pipeline_common.offset_point`, including its geodesic calculation and
longitude wrapping. It does not claim to find the globally nearest marine pixel.

The marine search retains those four initial radii and extends them to 5, 10,
20, 30 and 50 km to accommodate the 0.05° marine grid and coastal NoData mask.
These are explicit operational choices, not ecological validity thresholds.
Both the nominal query radius and actual distance to the cell center must be
within the declared maximum. The audit records both: they are not interchangeable
when a displaced query falls within a relatively large grid cell. Repeated cells
are logged without downloading them again. Both statistics and both decades must
be available at the same cell.

```bash
python scripts/02_environmental_covariates/02J_RECOVER_BIOORACLE_SALINITY_NA.py \
  --base "$PREVIOUS_DIR/Metadata_reinicio_401_reconciliada.csv" \
  --previous-audit "$PREVIOUS_DIR/salinity_audit/sample_salinity_audit.csv" \
  --out "$METADATA_OUTPUT_DIR/Metadata_reinicio_401_reconciliada.csv" \
  --audit-dir "$METADATA_OUTPUT_DIR/salinity_audit" \
  --cache-dir "$ENVIRONMENTAL_CACHE/biooracle" \
  --radial-steps-m 250,500,1000,2000,5000,10000,20000,30000,50000 \
  --max-anchor-distance-m 50000
```

02J fills only rows where both salinity fields are missing. Every existing value,
other field, coordinate, accession and row position is preserved. It verifies that
the preceding audit matches the input metadata. Partially populated salinity pairs
fail rather than mixing source cells. Cached source responses from 02H are reused,
newly needed cells are fetched on demand, and `--offline` reproduces the result.
The shared environmental helper requires NumPy and pandas, already declared by the
project. The whole workflow enables this optional stage with
`BIOORACLE_FILL_MISSING_CARDINAL=1` alongside `BIOORACLE_ENABLE=1`; historical
workflow defaults remain unchanged.

Unlike the initial 02H screen, this explicitly requested regional extraction does
not block values because of inherited coordinate warnings or previous spatial
holds. It retains `previous_status`, `previous_reason` and the original evidence,
plus `recovery_note`, in the separate audit. A recovered numerical value does not
validate the coordinate, waterbody connectivity or sediment salinity. No original
coordinates are corrected. Missing/invalid coordinates remain NA, as do locations
without a complete result along the declared cardinal queries within the bound.
That last status does not assert an exhaustive absence of marine data in a circle.

Additional outputs are `cardinal_search_attempts.csv`, `recovered_accessions.csv`,
`registros_sin_salinidad.csv` and `cardinal_recovery_manifest.json`. The combined
`sample_salinity_audit.csv` retains all accessions and the earlier evidence. No
additional columns are added to the main metadata table. Preserve the preceding
metadata/audit and the updated cache for exact reconstruction of the filling step.

## Coverage figure

With Matplotlib and NumPy installed, render an external PNG/SVG from the audit:

```bash
python scripts/02_environmental_covariates/02I_PLOT_BIOORACLE_AUDIT.py \
  --audit "$METADATA_OUTPUT_DIR/salinity_audit/sample_salinity_audit.csv" \
  --land-geojson "$REFERENCE_DIR/ne_10m_land.geojson" \
  --out-prefix "$METADATA_OUTPUT_DIR/Salinidad_BioORACLE_v3_cobertura"
```

[Obtain Natural Earth land polygons](https://www.naturalearthdata.com/downloads/10m-physical-vectors/10m-land/)
separately and retain their checksum. Coordinate points and accession counts are
reported separately. The figure is a coverage summary; individual coastal reviews
require the detailed maps and evidence retained in the data companion.
