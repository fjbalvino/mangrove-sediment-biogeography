# Biogeochemical cycles

Start with [methods, input contracts and commands](../../docs/BIOGEOCHEMICAL_CYCLES.md).

Run through `workflow/07_RUN_BIOGEOCHEMICAL_CYCLES.py`: evaluate rules (11A),
verify tables (11C), then optionally render figures locally (11B). Evaluation is
offline, uses the existing sample manifest by default, and records missing KO rows
as unknown. All biological inputs and generated outputs remain outside Git.

The [curated configuration](../../config/biogeochemical_cycles) contains operational
analysis rules rather than an annotation database. The current stage describes
community genetic potential; it does not estimate activity, flux or interactions.
