#!/usr/bin/env bash
# Reproduce core composition and distance-decay for the 343 profiles with assigned marine realms.
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export BIOGEOGRAPHY_BASE="${BIOGEOGRAPHY_BASE:-$(cd -- "$script_dir/../.." && pwd)}"
export BIOGEOGRAPHY_RESULTS="${BIOGEOGRAPHY_RESULTS:-$BIOGEOGRAPHY_BASE/results/taxonomic_core}"
export MAGENTA_CORE_PHYLOSEQ="${MAGENTA_CORE_PHYLOSEQ:-$BIOGEOGRAPHY_BASE/data/taxonomic_core_counts.rds}"
export MAGENTA_METADATA_CSV="${MAGENTA_METADATA_CSV:-$BIOGEOGRAPHY_BASE/data/environmental_metadata.csv}"
export CORE_SAMPLE_MANIFEST="${CORE_SAMPLE_MANIFEST:-$BIOGEOGRAPHY_BASE/config/samples/taxonomic_core_assigned_realm_identifiers.txt}"
RSCRIPT_BIN="${RSCRIPT_BIN:-Rscript}"
# Stage 05C must use the classification produced by this run.
unset RUN_09A
for input in "$MAGENTA_CORE_PHYLOSEQ" "$MAGENTA_METADATA_CSV" "$CORE_SAMPLE_MANIFEST"; do
  if [[ ! -f "$input" ]]; then
    printf 'Required input not found: %s\n' "$input" >&2
    exit 1
  fi
done
if ! command -v "$RSCRIPT_BIN" >/dev/null 2>&1; then
  printf 'Rscript executable not found: %s\n' "$RSCRIPT_BIN" >&2
  exit 1
fi
for stage in \
  05A_CLASSIFY_TAXONOMIC_CORE_ASSOCIATED.R \
  05B_CLASSIFY_REALM_RESTRICTED_TAXA.R \
  05C_SUMMARIZE_CORE_ASSOCIATED_PHYLA.R \
  05D_CORE_NONCORE_DISTANCE_DECAY.R; do
  "$RSCRIPT_BIN" "$script_dir/$stage"
done
