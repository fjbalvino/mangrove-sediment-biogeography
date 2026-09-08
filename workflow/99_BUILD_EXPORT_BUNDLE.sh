#!/usr/bin/env bash

# Assemble the final shareable code + canonical-data bundle on landalab.
# Default is dry-run. Raw FASTQ and per-read classifications are not copied.

set +e

PACKAGE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG="$PACKAGE_ROOT/config/landalab.env"
EXECUTE=0
OUTPUT_PARENT="/home/fjbalvino/Biogeografia_2026/export"

while [ "$#" -gt 0 ]; do
  case "$1" in
    --config) CONFIG="$2"; shift 2 ;;
    --output-parent) OUTPUT_PARENT="$2"; shift 2 ;;
    --execute) EXECUTE=1; shift ;;
    *) printf "Argumento desconocido: %s\n" "$1"; exit 2 ;;
  esac
done

if [ ! -f "$CONFIG" ]; then
  printf "ERROR: no existe config: %s\n" "$CONFIG"
  exit 2
fi

# shellcheck disable=SC1090
source "$CONFIG"

STAMP="$(date -u +%Y%m%d_%H%M%S)"
BUNDLE_NAME="MAGENTA_GLOBAL_BIOGEOGRAPHY_REPRODUCIBILITY_v1.0.0_${STAMP}"
BUNDLE_DIR="$OUTPUT_PARENT/$BUNDLE_NAME"
ARCHIVE="$OUTPUT_PARENT/${BUNDLE_NAME}.tar.gz"

INPUTS=(
  "$PHYLOSEQ_405"
  "$PHYLOSEQ_401"
  "$MAGENTA_METADATA_CURATED"
  "$MAGENTA_METADATA_COMPLETE"
  "$ENVIRONMENTAL_METADATA_FINAL"
  "$ICAMP_GENUS_TREE"
  "$MAGENTA_KO_MATRIX"
)

STATUS=0
printf "\n=== EXPORT INPUTS ===\n"
for FILE in "${INPUTS[@]}"; do
  if [ -f "$FILE" ]; then
    printf "PASS  %s\n" "$FILE"
  else
    printf "FAIL  %s\n" "$FILE"
    STATUS=1
  fi
done

if [ -f "$ICAMP_GENUS_TREE" ]; then
  ACTUAL_TREE_SHA256="$(sha256sum "$ICAMP_GENUS_TREE" | awk '{print $1}')"
  if [ "$ACTUAL_TREE_SHA256" != "$ICAMP_GENUS_TREE_SHA256" ]; then
    printf "FAIL  tree checksum expected=%s observed=%s\n" "$ICAMP_GENUS_TREE_SHA256" "$ACTUAL_TREE_SHA256"
    STATUS=1
  else
    printf "PASS  tree checksum\n"
  fi
fi

printf "\nBundle directory: %s\n" "$BUNDLE_DIR"
printf "Archive: %s\n" "$ARCHIVE"

if [ "$STATUS" -ne 0 ]; then
  printf "EXPORT BLOCKED: faltan inputs o falló un checksum.\n"
  exit "$STATUS"
fi

if [ "$EXECUTE" -ne 1 ]; then
  printf "DRY RUN: use --execute después de revisar los inputs.\n"
  exit 0
fi

if [ -e "$BUNDLE_DIR" ] || [ -e "$ARCHIVE" ]; then
  printf "ERROR: el destino ya existe; no se sobrescribirá.\n"
  exit 3
fi

mkdir -p "$BUNDLE_DIR/code" "$BUNDLE_DIR/data/phyloseq" "$BUNDLE_DIR/data/metadata" "$BUNDLE_DIR/data/trees" "$BUNDLE_DIR/data/functional" "$BUNDLE_DIR/provenance"
cp -a "$PACKAGE_ROOT/." "$BUNDLE_DIR/code/"
cp --preserve=mode,timestamps "$PHYLOSEQ_405" "$BUNDLE_DIR/data/phyloseq/"
cp --preserve=mode,timestamps "$PHYLOSEQ_401" "$BUNDLE_DIR/data/phyloseq/"
cp --preserve=mode,timestamps "$MAGENTA_METADATA_CURATED" "$BUNDLE_DIR/data/metadata/"
cp --preserve=mode,timestamps "$MAGENTA_METADATA_COMPLETE" "$BUNDLE_DIR/data/metadata/"
cp --preserve=mode,timestamps "$ENVIRONMENTAL_METADATA_FINAL" "$BUNDLE_DIR/data/metadata/"
cp --preserve=mode,timestamps "$ICAMP_GENUS_TREE" "$BUNDLE_DIR/data/trees/"
cp --preserve=mode,timestamps "$MAGENTA_KO_MATRIX" "$BUNDLE_DIR/data/functional/"

(
  cd "$BUNDLE_DIR" || exit 1
  find code data -type f -print0 | sort -z | xargs -0 sha256sum > provenance/data_and_code_files.sha256
)
STATUS=$?
if [ "$STATUS" -ne 0 ]; then
  printf "FAILED: no se pudo generar el inventario.\n"
  exit "$STATUS"
fi

tar -C "$OUTPUT_PARENT" -czf "$ARCHIVE" "$BUNDLE_NAME"
STATUS=$?
if [ "$STATUS" -ne 0 ]; then
  printf "FAILED: tar exit=%s\n" "$STATUS"
  exit "$STATUS"
fi

sha256sum "$ARCHIVE" > "${ARCHIVE}.sha256"
printf "PASS: %s\n" "$ARCHIVE"
printf "La terminal permanece abierta.\n"
