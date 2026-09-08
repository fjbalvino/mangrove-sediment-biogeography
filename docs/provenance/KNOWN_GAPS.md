# Known gap and external dependencies

## Missing code

The complete **MEGAHIT → Prodigal → MMseqs2 → eggNOG-mapper → KO matrix**
chain was not present among the supplied scripts. This is the only missing
Methods analysis in the canonical sequence and the only deliberately empty
numbered block (`07`). No placeholder implementation was fabricated.

The validated GTDB R226 genus tree is treated as a frozen external input to the
iCAMP workflow, not as a second missing analysis block. Its exact file and
checksum must accompany the Zenodo release. Descriptive figure assembly and
historical analyses that are not part of the current Methods workflow are not
counted as reproducibility gaps in this canonical repository.

## External assets that must be frozen

- Kraken2 PlusPF RefSeq January 2023 database and NCBI taxonomy dump.
- GTDB R226 genus tree used for iCAMP.
- Human Footprint 2020 raster, WorldClim 2.1 BIO rasters and MEOW vector release.
- Cached GBIF and SoilGrids responses with retrieval dates and checksums.
- Canonical phyloseq objects, complete/curated metadata and the KO count matrix.

## Metadata/cohort distinction

The functional atlas/network scripts use the 194-profile assembly-quality
cohort. The final regional functional analysis uses
`Metadata_2026_completa.csv` and excludes `SRR31095011`, `SRR31095534` and
`SRR31095535`, leaving 191 profiles and 17,482 non-zero KOs. Do not overwrite
the historical 194-profile inputs to make these counts appear identical.
