# Recover sampling locations and environmental observations

Use publications and supplements to recover locations before accepting a gazetteer
point. Link `run_accession` to `experiment_accession`, then to the author's sample
label. A geographic-name match by itself does not identify the sampling station.

## Evidence hierarchy

1. Published coordinates linked to the original sample/site, including verbatim
   representation, source, datum status and uncertainty when reported.
2. Georeferenced sampling-map markers, using identified control points, projection,
   transformation, independent residual checks, marker extent and uncertainty.
3. A locality footprint from a gazetteer with documented point-radius uncertainty.
4. A broad locality name without a defensible footprint: position remains missing.

Decimal places are not positional accuracy. Unknown uncertainty must remain
unknown, never zero. Do not convert an administrative centroid or a harbor/bay
representative point into a sampled sediment position. Marine Regions entries
that cite GeoNames are not an independent second coordinate observation.

For raster extraction, assess spatial support for each variable. Preserve original
and extraction coordinates separately, raster version/date/resolution, coastal
anchor distance and uncertainty. A nearby land pixel can represent different
substrate or land use. SoilGrids is not a direct sediment measurement; climate
normals are not the temperature at sampling. Do not replace historical GBIF
metrics by querying only a subset under a changed rarefaction rule.

## Literature recovery on 2026-10-30

The seven records without coordinates were examined. Al-Amoudi et al. (2016)
provide a site coordinate for Rabigh harbor lagoon. The archived experiment
alias MNR links the relevant record to that location; the similarly named
Al Kharrar lagoon must not be substituted. Datum and quantitative uncertainty
were not established from the consulted text. This is a published site coordinate,
not proof of centimeter-level GPS accuracy.

Jin et al. (2019) link six experiment IDs to sample names. Figure 1 distinguishes
five locations: two Guangxi positions (the two 2016 collections share one depicted
position) and three Hainan positions. The figure has scale bars but no coordinate
grid; the supplement has no numeric coordinates. The five figure markers were subsequently registered to georeferenced imagery
(see the figure-registration check below). They are approximate positions, not
reported GPS coordinates.
The six samples remain in the 191-profile inventory.

Table S1 provides temperature, pH, salinity, TN, AN, TC and TOC, plus month/year.
Extract values by sample label, never row order. Preserve TN as the table label:
the methods use inconsistent nitrogen terminology. Do not silently rename TN to
nitrite or harmonize it as total nitrogen until that discrepancy is resolved.
Retain measured TOC separately from SoilGrids SOC. These observations do not fill
the existing WorldClim/SoilGrids/footprint/GBIF columns.

The same paper describes viral-particle purification and MDA. Record that
preparation as technical provenance without changing the accepted sediment
classification or deleting samples. Quantitative comparisons with bulk DNA
require an explicit assessment of preparation effects; the authors themselves
interpret their viromes qualitatively.

## Reproducible merge

The curated evidence table and generated metadata stay outside GitHub. Required
evidence columns are documented in the script's `REQUIRED` constant. In that
separate table, newly collected attributes have the `recovery_` prefix; the two
identity keys are unprefixed. None of these columns is added to the metadata.
Allowed statuses are `published_site_coordinate`,
`locality_resolved_position_pending` and `figure_georeferenced_approximate`.
The latter requires `--allow-approximate`, a positive externally documented
uncertainty radius, its interpretation and registration evidence. This is an
explicit scientific processing option, not a request for user permission. Gazetteer context coordinates must remain
separate `recovery_gazetteer_*` evidence attributes. The script does not import
these into the metadata or promote them to sampling coordinates.

```bash
python scripts/12_gene_environment/12B_APPLY_METADATA_RECOVERY.py \
  --metadata /path/to/metadata_functional.csv \
  --evidence /path/to/Recuperacion_literatura_7_muestras.csv \
  --out-dir /path/to/new_recovery_output
```

The output `metadata_functional.csv` preserves the original column names and order,
row order, sample IDs, all existing environmental values, locality labels and
site/region assignments. It fills only previously absent coordinate pairs with
the specified published-coordinate status. Evidence remains in its separate input
table; input/output hashes remain in an external audit JSON. No new variables or
explanations are inserted into the metadata. Historical site labels and basis
columns remain historical: consult the external evidence table for recovered records. Rebuild derived
spatial groups/region assignments during the next explicit spatial processing
stage, rather than silently editing them here. No environmental extraction or
gene-environment association has been run by this recovery step.

## Figure-registration check, 2026-10-30

Figure 1 panels B and C were registered separately against Esri World Imagery
exported in EPSG:3857. The actual returned extent and dimensions, not the
requested bounding box, convert image pixels to projected coordinates. Stable
river junctions, coastlines and land features were compared visually. Initial
feature-only fits that did not align the coast were rejected. SIFT/RootSIFT,
geographic initialization and RANSAC supplied candidate matches (OpenCV 4.11).
The final affine fits use 48 unique image correspondences for Guangxi and 22 for
Hainan. Inputs, matched points, source hashes, export metadata and overlays stay
with the external evidence, not in GitHub.

Internal fit RMSE is approximately 39 m and 33 m; maximum leave-one-out errors
are 99 m and 80 m. These are registration diagnostics of selected image matches,
**not independent sampling-location validation or GPS accuracy**. The final
Guangxi registration is stable across image scales 0.5, 0.75 and 1.0. The colored
markers themselves cover approximate radii of 420–620 m. Their centers were
computed from their colored pixels, then transformed to WGS84 and rounded to
three decimals. The two 2016 Guangxi samples retain one depicted position; the
2015 position is about 2.7 km away. No additional spatial replicate is invented.

A **provisional 1 km operational radius** accompanies each recovered coordinate
in the external evidence. It allows sensitivity analysis but is neither an
empirically calibrated uncertainty bound nor a 95% confidence interval. Unknown
original marker placement, field GPS datum, composite-core extent and reference
imagery accuracy prevent such a claim. Source-map marker centers can fall on
water or urban land near the described river/coast; they must not be snapped to a
convenient soil or mangrove pixel. The figures establish locality-scale positions,
not exact sediment sampling points. Preserve the accepted 191-sample cohort.

Use these positions for locality-scale mapping and evaluate positional sensitivity
before environmental extraction. SoilGrids, GBIF and human-footprint values must
not be assigned from a single uncertain center as if the point were exact.
Existing environmental values and historical spatial assignments are unchanged.
The former `coordinates_unavailable` basis fields remain historical provenance;
the external recovery table identifies the new approximate coordinates. Rebuild
site/region assignments in the explicit spatial-processing stage.

Recompute coordinates and internal diagnostics from the archived registration:

```bash
python scripts/12_gene_environment/12C_RECOMPUTE_FIGURE_GEOREFERENCES.py \
  --evidence /path/to/georeferencing_evidence.json \
  --output /path/to/recomputed_coordinates.json

python scripts/12_gene_environment/12B_APPLY_METADATA_RECOVERY.py \
  --metadata /path/to/metadata_functional.csv \
  --evidence /path/to/Recuperacion_coordenadas_6_muestras.csv \
  --allow-approximate --out-dir /path/to/new_recovery_output
```

Do not rerun the merge on the already updated base: it correctly refuses to
replace existing coordinates. The six-row evidence must be applied to the prior
version that still lacks those coordinate pairs. No environmental raster values
or gene–environment associations were calculated in this check.

## Sources

- Al-Amoudi S et al. (2016). Metagenomics as a preliminary screen for antimicrobial
  bioprospecting. *Gene* 594:248–258. https://doi.org/10.1016/j.gene.2016.09.021
- Jin M et al. (2019). Diversities and potential biogeochemical impacts of mangrove
  soil viruses. *Microbiome* 7:58. https://doi.org/10.1186/s40168-019-0675-9
  Figure 1 and Additional file 1, Table S1.
- Chapman AD, Wieczorek JR (2020; online version 1.2.1, 2025).
  *Georeferencing Best Practices*. https://doi.org/10.15468/doc-gg7h-s853
  https://docs.gbif.org/georeferencing-best-practices/1.0/en/
- Marine Regions, Yalong Bay, MRGID 33088:
  https://www.marineregions.org/gazetteer.php?p=details&id=33088
- GeoNames, Tielu Gang, ID 1887384:
  https://www.geonames.org/1887384/tielu-gang.html

- Esri World Imagery: https://services.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer
- ArcGIS export-map extent semantics: https://developers.arcgis.com/rest/services-reference/enterprise/export-map/

Sources consulted 2026-10-30. Numerical evidence and source hashes accompany
the external recovery deliverable, not this code repository.
