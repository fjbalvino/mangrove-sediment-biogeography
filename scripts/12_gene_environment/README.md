# Gene–environment diagnostics and the eight-cycle atlas

Use `workflow/08_RUN_GENE_BIOGEOGRAPHY.py` from the repository root. The `template/` tree is executable source staged into a new local run directory, preserving the nitrogen/remaining-cycle separation and relative paths of the reviewed analytical package. It contains no measured matrices or generated figures/viewers.

The driver validates and joins external inputs by run_accession, copies the source/configuration snapshot, and runs the requested stage. It keeps the model and visualization source with each run, preventing a later checkout change from silently altering a running calibration. See `docs/GENE_BIOGEOGRAPHY.md` and `docs/ATLAS_PUBLICATION.md`.
