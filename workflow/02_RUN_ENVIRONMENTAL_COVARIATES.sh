#!/usr/bin/env bash

# Safe orchestrator for environmental, anthropogenic, biodiversity, and MEOW
# covariates. Default is dry-run; pass --execute only after the exact archived
# raster/vector assets have been declared in config/landalab.env.

set +e

PACKAGE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG="$PACKAGE_ROOT/config/landalab.env"
BASE_METADATA=""
EXECUTE=0

while [ "$#" -gt 0 ]; do
  case "$1" in
    --config)
      CONFIG="$2"
      shift 2
      ;;
    --base-metadata)
      BASE_METADATA="$2"
      shift 2
      ;;
    --execute)
      EXECUTE=1
      shift
      ;;
    *)
      printf "Argumento desconocido: %s\n" "$1"
      exit 2
      ;;
  esac
done

if [ ! -f "$CONFIG" ]; then
  printf "ERROR: no existe config: %s\n" "$CONFIG"
  exit 2
fi

# shellcheck disable=SC1090
source "$CONFIG"
export PATH="/data/bin:$(dirname "$PYTHON_BIN"):$(dirname "$R_BIN"):$PATH"

if [ -z "$BASE_METADATA" ]; then
  BASE_METADATA="$MAGENTA_METADATA_CURATED"
fi

STAMP="$(date -u +%Y%m%d_%H%M%S)"
RUN_DIR="$ENVIRONMENTAL_RESULTS/02_ENVIRONMENTAL_COVARIATES_${STAMP}"
TABLE_DIR="$RUN_DIR/tables"
DIAG_DIR="$RUN_DIR/diagnostics"
MANIFEST_DIR="$RUN_DIR/manifests"
LOG_DIR="$RUN_DIR/logs"

mkdir -p "$TABLE_DIR" "$DIAG_DIR" "$MANIFEST_DIR" "$LOG_DIR" "$ENVIRONMENTAL_CACHE"

SET401_METADATA="$TABLE_DIR/02A_metadata_set401.csv"
HII_METADATA="$TABLE_DIR/02B_metadata_human_footprint.csv"
CLIMATE_METADATA="$TABLE_DIR/02C_metadata_worldclim.csv"
GBIF_METADATA="$TABLE_DIR/02D_metadata_gbif.csv"
SOIL_METADATA="$TABLE_DIR/02E_metadata_soilgrids.csv"
MEOW_METADATA="$TABLE_DIR/02F_metadata_environmental_complete.csv"

run_step() {
  STEP="$1"
  shift
  printf "\n=== %s ===\n" "$STEP"
  printf "COMMAND:"
  printf " %q" "$@"
  printf "\n"
  if [ "$EXECUTE" -ne 1 ]; then
    return 0
  fi
  "$@" 2>&1 | tee "$LOG_DIR/${STEP}.log"
  RC="${PIPESTATUS[0]}"
  if [ "$RC" -ne 0 ]; then
    printf "FAILED: %s (exit=%s)\n" "$STEP" "$RC"
    return "$RC"
  fi
  printf "PASS: %s\n" "$STEP"
  return 0
}

RC=0

if [ "$RC" -eq 0 ]; then
  run_step "02A_PREPARE_SET401_METADATA" \
    "$PYTHON_BIN" "$PACKAGE_ROOT/scripts/02_environmental_covariates/02A_PREPARE_SET401_METADATA.py" \
    --base "$BASE_METADATA" \
    --out "$SET401_METADATA" \
    --audit-out "$DIAG_DIR/02A_set401_filter_audit.csv" \
    --manifest-out "$MANIFEST_DIR/02A_set401_manifest.json" \
    --expected-input-rows "$EXPECTED_PHYLOSEQ_405_SAMPLES" \
    --expected-output-rows "$EXPECTED_PHYLOSEQ_401_SAMPLES"
  RC=$?
fi

if [ "$RC" -eq 0 ]; then
  run_step "02B_HUMAN_FOOTPRINT" \
    "$PYTHON_BIN" "$PACKAGE_ROOT/scripts/02_environmental_covariates/02B_EXTRACT_HUMAN_FOOTPRINT_2020.py" \
    --base "$SET401_METADATA" \
    --human-footprint-tif "$HUMAN_FOOTPRINT_2020_TIF" \
    --out "$HII_METADATA" \
    --diag-out "$DIAG_DIR/02B_human_footprint_diagnostics.csv" \
    --manifest-out "$MANIFEST_DIR/02B_human_footprint_manifest.json" \
    --max-anchor-distance-m "$COASTAL_ANCHOR_MAX_M" \
    --overwrite
  RC=$?
fi

WORLDCLIM_COMMAND=(
  "$PYTHON_BIN" "$PACKAGE_ROOT/scripts/02_environmental_covariates/02C_EXTRACT_WORLDCLIM_BIOCLIM_1970_2000.py"
  --base "$HII_METADATA"
  --clim-dir "$WORLDCLIM_BIO_DIR"
  --out "$CLIMATE_METADATA"
  --diag-out "$DIAG_DIR/02C_worldclim_diagnostics.csv"
  --manifest-out "$MANIFEST_DIR/02C_worldclim_manifest.json"
  --max-anchor-distance-m "$COASTAL_ANCHOR_MAX_M"
  --overwrite
)
if [ -n "$WORLDCLIM_BIO_PREFIX" ]; then
  WORLDCLIM_COMMAND+=(--prefix "$WORLDCLIM_BIO_PREFIX")
fi

if [ "$RC" -eq 0 ]; then
  run_step "02C_WORLDCLIM" "${WORLDCLIM_COMMAND[@]}"
  RC=$?
fi

if [ "$RC" -eq 0 ]; then
  run_step "02D_GBIF_PLANT_DIVERSITY" \
    "$PYTHON_BIN" "$PACKAGE_ROOT/scripts/02_environmental_covariates/02D_EXTRACT_GBIF_PLANT_DIVERSITY_2KM.py" \
    --base "$CLIMATE_METADATA" \
    --out "$GBIF_METADATA" \
    --cache-dir "$ENVIRONMENTAL_CACHE/gbif" \
    --diag-out "$DIAG_DIR/02D_gbif_diagnostics.csv" \
    --manifest-out "$MANIFEST_DIR/02D_gbif_manifest.json" \
    --rarefaction-n min \
    --rarefaction-reps 100 \
    --seed "$RANDOM_SEED" \
    --max-workers "$GBIF_WORKERS" \
    --global-rps "$GBIF_GLOBAL_RPS" \
    --overwrite
  RC=$?
fi

if [ "$RC" -eq 0 ]; then
  run_step "02E_SOILGRIDS_TOPSOIL" \
    "$PYTHON_BIN" "$PACKAGE_ROOT/scripts/02_environmental_covariates/02E_EXTRACT_SOILGRIDS_TOPSOIL_0_30CM.py" \
    --base "$GBIF_METADATA" \
    --out "$SOIL_METADATA" \
    --cache-dir "$ENVIRONMENTAL_CACHE/soilgrids" \
    --diag-out "$DIAG_DIR/02E_soilgrids_diagnostics.csv" \
    --manifest-out "$MANIFEST_DIR/02E_soilgrids_manifest.json" \
    --properties clay,soc,sand,silt,bdod,cec,cfvo,phh2o \
    --max-anchor-distance-m "$COASTAL_ANCHOR_MAX_M" \
    --global-rps "$SOILGRIDS_GLOBAL_RPS" \
    --overwrite
  RC=$?
fi

if [ "$RC" -eq 0 ]; then
  run_step "02F_ASSIGN_MEOW" \
    "$PYTHON_BIN" "$PACKAGE_ROOT/scripts/02_environmental_covariates/02F_ASSIGN_MEOW_REGIONS.py" \
    --base "$SOIL_METADATA" \
    --meow "$MEOW_VECTOR" \
    --out "$MEOW_METADATA" \
    --diag-out "$DIAG_DIR/02F_meow_diagnostics.csv" \
    --manifest-out "$MANIFEST_DIR/02F_meow_manifest.json" \
    --max-anchor-distance-m "$COASTAL_ANCHOR_MAX_M" \
    --overwrite
  RC=$?
fi

if [ "$RC" -eq 0 ]; then
  run_step "02G_VALIDATE_METADATA" \
    "$PYTHON_BIN" "$PACKAGE_ROOT/scripts/02_environmental_covariates/02G_VALIDATE_ENVIRONMENTAL_METADATA.py" \
    --metadata "$MEOW_METADATA" \
    --out-dir "$RUN_DIR/02G_validation" \
    --expected-rows "$EXPECTED_PHYLOSEQ_401_SAMPLES"
  RC=$?
fi

if [ "$EXECUTE" -eq 1 ] && [ "$RC" -eq 0 ]; then
  if [ -e "$ENVIRONMENTAL_METADATA_FINAL" ]; then
    printf "ERROR: el destino final ya existe y no será sobrescrito: %s\n" "$ENVIRONMENTAL_METADATA_FINAL"
    RC=3
  else
    cp --preserve=mode,timestamps "$MEOW_METADATA" "$ENVIRONMENTAL_METADATA_FINAL"
    RC=$?
  fi
fi

if [ "$EXECUTE" -eq 1 ] && [ "$RC" -eq 0 ]; then
  printf "%s\n" "$RUN_DIR" > "$ENVIRONMENTAL_RESULTS/LATEST_02_ENVIRONMENTAL_COVARIATES.txt"
  printf "PASS\n" > "$RUN_DIR/STATUS.txt"
elif [ "$EXECUTE" -eq 1 ]; then
  printf "FAILED\nexit_code=%s\n" "$RC" > "$RUN_DIR/STATUS.txt"
else
  printf "\nDRY RUN: no se ejecutó ninguna etapa. Complete primero las rutas de los activos ambientales.\n"
fi

printf "\nRUN_DIR=%s\n" "$RUN_DIR"
printf "FINAL_METADATA=%s\n" "$ENVIRONMENTAL_METADATA_FINAL"
printf "STATUS_CODE=%s\n" "$RC"
printf "La terminal permanece abierta.\n"
exit "$RC"
