# Gene environment design audit

`12A_AUDIT_GENE_ENVIRONMENT_DESIGN.py` audits the support for the proposed
[functional gene biogeography protocol](../../docs/GENE_ENVIRONMENT_BIOGEOGRAPHY.md).
It does not fit environmental effects, infer biological independence from
coordinates, convert units, impute values or modify inputs.

Required metadata fields are `run_accession`, `study_id`, `sampling_site_id`,
`lat`, `lon` and `REALM`. `sampling_site_id` is the existing provisional
coordinate-group label. `sample_accession`, when available, identifies repeated
BioSamples. Required matrix values are finite nonnegative integer counts of
predicted protein records assigned to KOs; the first column is the unique KO ID.
The explicit sample manifest determines the selected profiles and their order.

The auditor writes external TSV summaries and an audit JSON with input hashes,
software versions, environmental coverage, within-study variation, repeated
BioSamples, coordinate-value conflicts and algebraic design ranks. Algebraic rank
does not establish model validity, biological replication or statistical power.
Output must be a new or empty directory. The unit flag is an explicit provenance
declaration; integer values alone cannot establish the measurement unit.

Use the command in the protocol and run `make gene-environment-tests` for the
four synthetic input-integrity checks. `make validate` includes these checks.
Association and atlas implementation will follow the conditions in
`config/gene_environment/analysis_plan.json`; they are not implemented by 12A.

`12B_APPLY_METADATA_RECOVERY.py` matches external literature evidence by run and
experiment ID and fills only absent coordinate pairs documented as published
sampling-site coordinates. It preserves the original columns and their order;
evidence and direct observations stay in the external evidence table. Gazetteer
representative points are not imported. See [recovery protocol](../../docs/METADATA_GEOREFERENCING_RECOVERY.md).
