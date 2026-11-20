# Mangrove landscape cover and regional tidal range

These independent environmental stages preserve every input row, accession,
coordinate, existing value and relative column order. They do not select a cohort
or rerun manuscript models. Keep tables, maps, model atlases and credentials
outside GitHub. Use the full restart metadata directly: the earlier 02A step
deliberately selects the historical taxonomic cohort.

## 02K: mapped mangrove cover in a 2 km buffer

`mangrove_cover_2020_2km_pct` is 100 times the union area of mapped mangrove
polygons inside a 2,000 m buffer, divided by the entire buffer area, including
both land and water. It is landscape extent, not canopy closure, stem density,
biomass, plant richness, measured productivity or percentage of land alone.

The frozen source is **Global Mangrove Watch 4.1.12, year 2020**:

- DOI: [10.5281/zenodo.21346457](https://doi.org/10.5281/zenodo.21346457).
- File: `gmw_v4112_2020_mng_ext_cntry_info_vec.gpkg.gz`.
- Publisher MD5: `5d75f466c31f6a48e54d616ee1b75c91`; 167,847,900 bytes.
- License: CC BY 4.0. See the record for full author attribution and methods.
- The publisher now recommends v4.1 in place of v3.0. The extractor pins this
  particular release instead of silently following a future latest version.

2020 is a common landscape reference, also used by the project's human-footprint
layer. It does not claim to describe vegetation on each sampling date. The public
release includes an annual series, but selecting a sample-specific year would
require reliable collection dates, currently absent from the master columns.

The buffer boundary is built on the WGS84 ellipsoid using 1,440 geodesic bearings.
Intersection and union areas are measured in a local Lambert azimuthal equal-area
projection. Overlapping features and country-boundary splits cannot double-count
area. The national accuracy-adjustment factor in the publisher's country summaries
is not applied to individual local polygons. Geometry repairs are recorded.

**Zero is a valid result:** no mangrove is mapped in that fixed buffer. It is not
field-confirmed absence; omissions, the product habitat mask and inaccurate sample
coordinates remain possible. There is no offset search or movement toward mapped
vegetation. Missing/invalid coordinates remain NA. Existing coordinate warnings
are copied to the audit without changing the coordinates or validating them anew.

```bash
python -m pip install -r config/coastal_requirements.txt
python scripts/02_environmental_covariates/02K_EXTRACT_GMW_COVER.py \
  --base "$INPUT/Metadata_reinicio_401_reconciliada.csv" \
  --out "$OUTPUT/Metadata_reinicio_401_reconciliada.csv" \
  --audit-dir "$OUTPUT/mangrove_cover_audit" \
  --cache-dir "$ENVIRONMENTAL_CACHE/gmw"
```

The script checks the publisher's archive hash before extraction and records
SHA-256 for the archive and expanded GeoPackage. Reuse that cache with `--offline`.
Download/read/checksum errors fail the run; they do not become zero or NA.
The main table receives only one new field, adjacent to human footprint. The
separate audit records area, denominator, radius, source, warnings and zero cases.

## 02L: FES2022b tidal range, external atlas required

This stage is implemented and tested, but **the real FES2022b atlas was unavailable
for the current extraction**. AVISO's official file endpoint returned HTTP 401
without an authorized account. No real sample tidal ranges were computed. The
pending-input audit must not be described as a completed tidal extraction.

Obtain the **34 non-extrapolated Cartesian ocean-height constituent files** from
`fes2022b/ocean_tide_20241025` through an authorized AVISO account. Use the `.nc`
files after decompressing `.nc.xz`. Ocean loading, extrapolated coast grids,
generic test grids and alternative tide models are not substituted.

- [Official model and access information](https://www.aviso.altimetry.fr/index.php?id=1279).
- DOI: [10.24400/527896/A01-2024.004](https://doi.org/10.24400/527896/A01-2024.004).
- [Official PyFES documentation](https://cnes.github.io/aviso-fes/).
- Validated engine: `pyfes==2026.5.2`, Darwin/FES defaults, two threads,
  astronomical angles recalculated at each time step, `dynamic: [A5]` as in the
  provider's FES2022 configuration. The ocean-height files must use cm and degrees.

Run preflight first:

```bash
python scripts/02_environmental_covariates/02L_EXTRACT_FES_TIDAL_RANGE.py \
  --base "$OUTPUT/Metadata_reinicio_401_reconciliada.csv" \
  --atlas-dir "$FES_ATLAS_DIR" \
  --audit-dir "$OUTPUT/tidal_range_audit" --preflight
```

Preflight writes one status per accession. It returns exit code 2 and
`AWAITING_INPUT` if files are absent, or `READY_NOT_EXTRACTED` if all paths exist.
It does not modify the master or add an entirely empty column. A ready preflight
only confirms file presence; full extraction additionally validates files and units.

After the atlas is available:

```bash
python scripts/02_environmental_covariates/02L_EXTRACT_FES_TIDAL_RANGE.py \
  --base "$OUTPUT/Metadata_reinicio_401_reconciliada.csv" \
  --out "$NEXT_OUTPUT/Metadata_reinicio_401_reconciliada.csv" \
  --atlas-dir "$FES_ATLAS_DIR" \
  --audit-dir "$NEXT_OUTPUT/tidal_range_audit"
```

`tidal_range_mean_m` is calculated from a complete 15-minute astronomical tide
series over 2000-01-01 inclusive to 2020-01-01 exclusive (UTC). It is the mean
absolute height difference between successive high/low turning points, with
three-point parabolic refinement of each extremum. Each rising/falling half-cycle
has equal weight. Endpoints are excluded; there is no arbitrary daily boundary,
monthly maximum, single-constituent amplitude or meteorological storm component.
The audit includes the 5th and 95th percentiles of the individual ranges and the
number of half-cycles. Those percentiles are variation, not uncertainty intervals.

Only full four-point native interpolation is accepted automatically. Missing or
extrapolated predictions remain NA; partial coastal interpolation needs explicit
review. No automatic cardinal search is used. An optional `--review FILE` can
approve an anchor within 5 km after checking the same coastline/estuary or hold an
assignment. Required review columns are:

```text
run_accession,lat,lon,decision,anchor_lat,anchor_lon,reason,evidence
```

`decision` is `approve` or `hold`; original coordinate strings must match. Every
decision needs evidence and a reason. A reviewed anchor still needs a finite,
non-extrapolated prediction. Five km is an operational bound, not proof of water
connectivity. Retain the atlas hashes, exact review, original warnings and audit.

The final result is a regional astronomical tide reference. It does not measure
plot hydroperiod, sediment elevation, channel connectivity, current velocity or
redox state, and it does not reconstruct historical coastal geometry. Missing
collection dates do not prevent a common astronomical reference period.

## Workflow and verification

The existing environmental workflow enables these stages independently through
`GMW_ENABLE=1` and `FES_TIDAL_ENABLE=1`; both default to zero. An enabled tide stage
with missing inputs fails rather than silently publishing an empty tidal result.
All normal stages still pass through the metadata validator.

Tests cover geodesic-buffer area, latitude, date-line query windows, polygon
overlap, clipping, zeros versus missing coordinates, field/order preservation,
the factor-of-two difference between amplitude and range, incomplete tide-series
rejection, reviewed-anchor bounds and missing-atlas behavior. A real PyFES engine
and a synthetic NetCDF atlas test the computation path. **Synthetic tests are not
validation against the unavailable real atlas or local tide gauges.**
