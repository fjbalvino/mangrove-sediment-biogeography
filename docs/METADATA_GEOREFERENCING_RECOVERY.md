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
grid; the supplement has no numeric coordinates. Sampling positions therefore
remain unresolved pending defensible map registration or author information.
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
Allowed status values are `published_site_coordinate` and
`locality_resolved_position_pending`. Gazetteer context coordinates must remain
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

Sources consulted 2026-10-30. Numerical evidence and source hashes accompany
the external recovery deliverable, not this code repository.
