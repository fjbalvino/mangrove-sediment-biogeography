#!/usr/bin/env bash
# Reproduce the functional atlas, regional tests, coupling and Figure 4.
set -euo pipefail
PACKAGE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
EXECUTE=0
CONFIG=""
while (($#)); do
  case "$1" in
    --execute) EXECUTE=1; shift ;;
    --config)
      if (($# < 2)); then echo "--config requires a file" >&2; exit 2; fi
      CONFIG="$2"; shift 2 ;;
    --help|-h)
      echo "Usage: $0 [--config FILE] [--execute]"
      echo "Preview the commands, or run all functional-atlas analyses with --execute."
      exit 0 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done
source "$PACKAGE_ROOT/workflow/lib/config.sh"
magenta_load_config "$CONFIG"
cd "$PACKAGE_ROOT"
PYTHON_BIN="${PYTHON_BIN:-python3}"
FUNCTIONAL_INPUT_DIR="${FUNCTIONAL_INPUT_DIR:-data/functional}"
FUNCTIONAL_OUTPUT_DIR="${FUNCTIONAL_OUTPUT_DIR:-results/functional}"
FUNCTIONAL_FIGURE_DIR="${FUNCTIONAL_FIGURE_DIR:-figures/functional}"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
run() {
  printf '%q ' "$@"
  printf '\n'
  if ((EXECUTE)); then "$@"; fi
}
if ((EXECUTE)); then
  for name in functional_KO_counts_original.tsv.gz metadata_functional.csv taxonomic_matched_counts.npz; do
    if [[ ! -f "$FUNCTIONAL_INPUT_DIR/$name" ]]; then
      echo "Missing functional input: $FUNCTIONAL_INPUT_DIR/$name" >&2
      exit 2
    fi
  done
fi
run "$PYTHON_BIN" scripts/08_functional_atlas/08A_ANALYSE_FUNCTIONAL_ATLAS.py \
  --inputs "$FUNCTIONAL_INPUT_DIR" --output "$FUNCTIONAL_OUTPUT_DIR" \
  --sample-manifest "$FUNCTIONAL_SAMPLE_MANIFEST"
run "$PYTHON_BIN" scripts/08_functional_atlas/08B_PLOT_FUNCTIONAL_ATLAS.py \
  --results "$FUNCTIONAL_OUTPUT_DIR" --output "$FUNCTIONAL_FIGURE_DIR"
