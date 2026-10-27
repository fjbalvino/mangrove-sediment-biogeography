#!/usr/bin/env bash
# Reproduce the manuscript Figure 6 cohort audit, balanced analyses and final figures.
set -euo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PACKAGE_ROOT="$REPO_ROOT"
CONFIG=""
INPUT_DIR=""
OUTPUT_DIR=""
MAP_FILE=""
EXECUTE=0
COMPLETE_PANELS=0
while (($#)); do
  case "$1" in
    --config|--inputs|--outdir|--world-map)
      [[ $# -ge 2 ]] || { echo "$1 requires a value" >&2; exit 2; }
      case "$1" in
        --config) CONFIG="$2" ;;
        --inputs) INPUT_DIR="$2" ;;
        --outdir) OUTPUT_DIR="$2" ;;
        --world-map) MAP_FILE="$2" ;;
      esac
      shift 2 ;;
    --help|-h) echo "Usage: $0 [--config FILE] [--inputs DIR] [--outdir DIR] [--world-map FILE] [--execute] [--complete-panels]"; exit 0 ;;
    --execute) EXECUTE=1; shift ;;
    --complete-panels) COMPLETE_PANELS=1; shift ;;
    *) printf 'Unknown argument: %s\n' "$1" >&2; exit 2 ;;
  esac
done
source "$PACKAGE_ROOT/workflow/lib/config.sh"
magenta_load_config "$CONFIG"
INPUT_DIR="${INPUT_DIR:-$FUNCTIONAL_INPUT_DIR}"
OUTPUT_DIR="${OUTPUT_DIR:-$FUNCTIONAL_REGIONAL_OUTPUT_DIR}"
MAP_FILE="${MAP_FILE:-$WORLD_MAP}"
RSCRIPT_BIN="${RSCRIPT:-$R_BIN}"
PYTHON_BIN="${PYTHON:-$PYTHON_BIN}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 TZ=UTC
BLOCK="$REPO_ROOT/scripts/09_functional_regional_differentiation"
MATRIX="$INPUT_DIR/functional_KO_counts_original.tsv.gz"
META="$INPUT_DIR/metadata_functional.csv"
ATLAS="$INPUT_DIR/KO_persistence_classes.tsv"
DICTIONARY="$INPUT_DIR/KO_function_dictionary.tsv"
ACCESSIONS="$FUNCTIONAL_SAMPLE_MANIFEST"
COMMAND=("$RSCRIPT_BIN" --vanilla "$BLOCK/09A_FUNCTIONAL_PREVALENCE_REGIONAL_DIFFERENTIATION.R"
  --base "$REPO_ROOT" --outdir "$OUTPUT_DIR"
  --ko_matrix "$MATRIX" --metadata "$META" --atlas_classes "$ATLAS" --accessions "$ACCESSIONS"
  --realm_col REALM --locality_col core_locality --country_col Country --study_col study_id
  --strict_expected 1 --strict_class_counts 1 --make_figures 0 --full_sample_tests 0
  --n_iter 499 --n_diagnostic_iter 499 --n_perm_balanced 499
  --n_sensitivity_iter 99 --n_equal_ko_iter 99 --min_occurrence_sensitivity 1,10 --seed 20261015)
if (( ! EXECUTE )); then
  printf 'Preview; add --execute to run.\n'
  printf '%q ' "${COMMAND[@]}"; printf '\n'
  printf 'External basemap: %s\n' "$MAP_FILE"
  printf 'Then: verify cohorts, calculate process scores and locality turnover, independently verify in R, and render manuscript Figure 6 (historical Figure_5 filenames).\n'
  if (( COMPLETE_PANELS )); then printf 'Existing balanced results will be verified and reused; only missing process/map analyses are recalculated.\n'; fi
  exit 0
fi
for f in "$MATRIX" "$META" "$ATLAS" "$DICTIONARY" "$ACCESSIONS" "$MAP_FILE"; do
  [[ -s "$f" ]] || { printf 'Missing input: %s\n' "$f" >&2; exit 1; }
done
mkdir -p "$OUTPUT_DIR"
if (( ! COMPLETE_PANELS )); then
  "${COMMAND[@]}" 2>&1 | tee "$OUTPUT_DIR/execution.log"
fi
"$PYTHON_BIN" "$BLOCK/09D_VERIFY_FUNCTIONAL_REGIONAL_OUTPUTS.py" \
  --tables "$OUTPUT_DIR/tables" --accessions "$ACCESSIONS" --metadata "$META" --counts "$MATRIX"
"$PYTHON_BIN" "$BLOCK/09C_AUDIT_FUNCTIONAL_CONTRIBUTIONS.py" \
  --tables "$OUTPUT_DIR/tables" --counts "$MATRIX" --metadata "$META" --dictionary "$DICTIONARY"
"$PYTHON_BIN" "$BLOCK/09E_FUNCTIONAL_PROCESSES_AND_MAP.py" \
  --counts "$MATRIX" --metadata "$META" --accessions "$ACCESSIONS" \
  --catalog "$REPO_ROOT/config/functional_process_catalog.tsv" \
  --tables "$OUTPUT_DIR/tables" --objects "$OUTPUT_DIR/rds"
"$RSCRIPT_BIN" --vanilla "$BLOCK/09F_VERIFY_FUNCTIONAL_PROCESS_MAP.R" \
  --counts "$MATRIX" --metadata "$META" --catalog "$REPO_ROOT/config/functional_process_catalog.tsv" \
  --tables "$OUTPUT_DIR/tables"
"$PYTHON_BIN" "$BLOCK/09B_PLOT_FUNCTIONAL_REGIONAL_DIFFERENTIATION.py" \
  --tables "$OUTPUT_DIR/tables" --outdir "$OUTPUT_DIR/plots" --world-map "$MAP_FILE"
printf 'Completed manuscript Figure 6 (historical filename): %s\n' "$OUTPUT_DIR/plots/Figure_5_functional_regional_identity_combined.png"
