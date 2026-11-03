#!/usr/bin/env bash
set -euo pipefail
: "${RECOVERY_WORKDIR:?Set an absolute directory outside git for private inputs, caches and outputs}"
: "${METADATA_MASTER:?Set the original master CSV path}"
: "${METADATA_FUNCTIONAL:?Set the original functional CSV path}"
repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
recovery_dir="$repo_dir/scripts/02_environmental_covariates/recovery"
python3 "$recovery_dir/prepare_recovery.py" --master "$METADATA_MASTER" --functional "$METADATA_FUNCTIONAL"
python3 "$recovery_dir/recover_worldclim30.py"
python3 "$recovery_dir/recover_rasters.py" human
python3 "$recovery_dir/recover_soil.py"
python3 "$recovery_dir/recover_gbif.py"
python3 "$recovery_dir/merge_verified.py"
# CSV identities are retained. This workflow never writes back to an input file.
