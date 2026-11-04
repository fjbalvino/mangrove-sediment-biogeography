# Gene-biogeography portability validation

Validation completed on 4 November 2026 against version 7 of the accepted reference package. The source archive identifier and SHA-256 are recorded in `config/gene_biogeography/source_provenance.json`. All reference data and generated outputs remained outside the repository.

| Check | Result |
| --- | --- |
| Historical 191-profile preparation | Passed: original selected cohort, ID matching and full-matrix normalization reproduced |
| Numerical reference comparison | Nine analytical tables agreed within absolute/relative tolerance `1e-10`: raw and normalized KO representation, geographic summaries, two detection thresholds, two rule-state tables, environmental diagnostics and model support |
| Cohort-size portability | Passed with 16 synthetic profiles and shuffled metadata containing an extra accession; only selected IDs entered analysis |
| Input and run integrity | Duplicate IDs, fractional normalized inputs, reused preparation directories and changed protected inputs were rejected |
| Non-estimable designs | Constant predictors were reported without manufacturing estimable coefficients |
| Fresh calibration preparation | Both family plans were built from current input dimensions and frozen signal-marker availability |
| Statistical implementation | Nitrogen's official R check passed; 240 remaining-panel fits gave 1,440 Python/R comparisons with identical finite/non-finite patterns and maximum scaled error below `1.3e-10` (acceptance tolerance `1e-5`) |
| Atlas generation | All eight panels, supplementary visualizations and output invariants passed |
| Viewer regression checks | All eight panels passed selection, overlapping locations, gene/colour updates, zoom, map-layer order, opacity and fallback checks under jsdom |
| Local website export | Nine HTML pages and supporting text exported successfully; internal links validated and original matrices/full metadata excluded |

The R comparison used R 4.3.3 and clubSandwich 0.5.10. Python dependencies match `config/gene_biogeography/requirements.txt`; viewer checks used jsdom 26.1.0. Non-numerical nitrogen labels were standardized to English.

This validates source portability and numerical implementation. **The full two-family simulation campaign was not rerun for this source update.** Each future input release must receive a fresh calibration plan, implementation check and explicit simulation run. The historical calibration outcome is not a result for Dani's future matrix. The DOM checks do not establish live terrain-provider availability or replace final visual inspection before public deployment.

The static export still embeds selected derived values, coordinates, environmental estimates and rendered images. It was inspected locally and was not uploaded or deployed. No input matrix, metadata table, image, generated viewer or simulation checkpoint is included in the source update.
