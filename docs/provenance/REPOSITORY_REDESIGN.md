# Repository redesign provenance

## Sources and preservation

- Original repository: [mangrove-sediment-biogeography](https://github.com/fjbalvino/mangrove-sediment-biogeography), main commit `dc848dc2377b54f9d737df7fdc2ed8cad86eebe1`.
- Before editing, all three branches and 15 reachable commits were preserved in a verified Git bundle and a ZIP with the original main snapshot. Git object hashes, bundle restoration and ZIP integrity were checked.
- Manuscript and five supplements reviewed on 26 October 2026. The [manuscript guide](../MANUSCRIPT_GUIDE.md) records the current title, organization and cohorts.
- Structural reference: [SongzeCHEN/MetaGenome-MAG-Analysis](https://github.com/SongzeCHEN/MetaGenome-MAG-Analysis), commit `cc206cd0c9c448e96bc352492027f3f49c4b56da`.

The reference inspired an overview, ordered workflow, prerequisites, installation, usage, output inventory and citation sections. Its MAG-binning/assembly code was not copied. Existing MAGENTA analytical filenames and algorithms were retained.

## REA inspection

REA `rea-agents@6.8.0` inspected a ZIP generated from the reference commit using `inspect-artifact --detail summary`.

- Archive SHA-256: `0b6f94c53759077538c48fe08098d4859c67f466e19b11d035a59b4a540d39c4`.
- Inventory: 12 files and two directories; no integrity contradictions reported.
- Parent evidence: `ev_0b1c67f0f77dd329b81406886e3d31deaf8b0c0f471a8c817e02a2dd421b4b1c`.
- Scope: bounded static archive inspection, followed by direct source review.
- Limitation: inspection did not execute the reference pipeline, establish scientific equivalence or rerun MAGENTA.

## Changes and data separation

Portable configuration replaces machine-specific defaults. Shell previews no longer create analysis directories. The network has an explicit preview option. The regional plotter receives its basemap as an external input.

Seventy-two tracked files were preserved outside the code tree: numerical summaries and their readmes, one basemap and an execution log. A migration manifest records each original path, destination, size and SHA-256. The original bytes are unchanged. This staging package is incomplete as a final data deposit: full inputs, analysis objects and publication figures still need consolidation.

The current code tree contains no biological/geospatial databases, numerical result tables or figure files. Tiny accession manifests, curated process rules and checksums remain as reproducibility configuration. Historical Git commits remain intact.

## Scope of verification

Repository validation checks text-only source boundaries, script indexing, source checksums, Python/Bash/R syntax and existing statistical contracts. Workflow smoke tests check previews and configuration from a different working directory. None of these checks certifies a full data-backed rerun.

The final data/figure deposit and its DOI remain pending. Update archive links and the manuscript availability statement when that deposit exists.
