#!/usr/bin/env bash

# Safe orchestrator for the taxonomic acquisition/profiling block.
# Default is dry-run. Pass --execute to start data acquisition and long jobs.

set +e

PACKAGE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG="$PACKAGE_ROOT/config/landalab.env"
EXECUTE=0
LOCAL_FASTQ_MAP=""

while [ "$#" -gt 0 ]; do
  case "$1" in
    --config)
      CONFIG="$2"
      shift 2
      ;;
    --local-fastq-map)
      LOCAL_FASTQ_MAP="$2"
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

STAMP="$(date -u +%Y%m%d_%H%M%S)"
RUN_DIR="$MAGENTA_PIPELINE_RESULTS/01_TAXONOMIC_PREPROCESSING_${STAMP}"
META_DIR="$RUN_DIR/metadata"
MANIFEST_DIR="$RUN_DIR/manifests"
QC_DIR="$RUN_DIR/fastqc_raw"
FASTP_REPORT_DIR="$RUN_DIR/fastp_reports"
LOG_DIR="$RUN_DIR/logs"

mkdir -p "$META_DIR" "$MANIFEST_DIR" "$QC_DIR" "$FASTP_REPORT_DIR" "$LOG_DIR"

METADATA_FROZEN="$META_DIR/01A_metadata_frozen_enriched.csv"
PUBLIC_METADATA="$META_DIR/01A_public_ENA_SRA_metadata.csv"
METADATA_CONFLICTS="$META_DIR/01A_curated_public_conflicts.csv"
RAW_MANIFEST="$MANIFEST_DIR/01B_raw_fastq_manifest.csv"
FASTQC_STATUS="$MANIFEST_DIR/01C_fastqc_raw_status.csv"
TRIM_MANIFEST="$MANIFEST_DIR/01D_trimmed_fastq_manifest.csv"
KRAKEN_MANIFEST="$MANIFEST_DIR/01E_kraken2_manifest.csv"
KRAKEN_COUNT_TABLE="$RUN_DIR/01F_kraken2_direct_otu_counts.tsv.gz"
KRAKEN_COUNT_AUDIT="$MANIFEST_DIR/01F_kraken2_direct_count_audit.tsv"

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
  run_step "01A_ACQUIRE_METADATA" \
    "$PYTHON_BIN" "$PACKAGE_ROOT/scripts/01_taxonomic_preprocessing/01A_ACQUIRE_COHORT_METADATA_ENA_SRA.py" \
    --curated-metadata "$MAGENTA_METADATA_CURATED" \
    --out "$METADATA_FROZEN" \
    --public-out "$PUBLIC_METADATA" \
    --conflicts-out "$METADATA_CONFLICTS" \
    --manifest-out "$META_DIR/01A_metadata_manifest.json" \
    --allow-missing-public
  RC=$?
fi

DOWNLOAD_COMMAND=(
  "$PYTHON_BIN" "$PACKAGE_ROOT/scripts/01_taxonomic_preprocessing/01B_DOWNLOAD_RAW_READS_ENA_SRA.py"
  --metadata "$METADATA_FROZEN"
  --out-dir "$MAGENTA_RAW_DIR"
  --manifest-out "$RAW_MANIFEST"
  --max-workers "$DOWNLOAD_WORKERS"
)
if [ -n "$LOCAL_FASTQ_MAP" ]; then
  DOWNLOAD_COMMAND+=(--local-fastq-map "$LOCAL_FASTQ_MAP")
fi

if [ "$RC" -eq 0 ]; then
  run_step "01B_DOWNLOAD_RAW_READS" "${DOWNLOAD_COMMAND[@]}"
  RC=$?
fi

if [ "$RC" -eq 0 ]; then
  run_step "01C_FASTQC_RAW" \
    "$PYTHON_BIN" "$PACKAGE_ROOT/scripts/01_taxonomic_preprocessing/01C_FASTQC_RAW_READS.py" \
    --manifest "$RAW_MANIFEST" \
    --out-dir "$QC_DIR" \
    --status-out "$FASTQC_STATUS" \
    --threads-per-job "$FASTQC_THREADS_PER_JOB" \
    --max-workers "$FASTQC_WORKERS"
  RC=$?
fi

if [ "$RC" -eq 0 ]; then
  run_step "01D_FASTP" \
    "$PYTHON_BIN" "$PACKAGE_ROOT/scripts/01_taxonomic_preprocessing/01D_TRIM_READS_FASTP.py" \
    --raw-manifest "$RAW_MANIFEST" \
    --out-dir "$MAGENTA_TRIMMED_DIR" \
    --manifest-out "$TRIM_MANIFEST" \
    --report-dir "$FASTP_REPORT_DIR" \
    --threads-per-sample "$FASTP_THREADS_PER_SAMPLE" \
    --max-workers "$FASTP_WORKERS"
  RC=$?
fi

if [ "$RC" -eq 0 ]; then
  run_step "01E_KRAKEN2" \
    "$PYTHON_BIN" "$PACKAGE_ROOT/scripts/01_taxonomic_preprocessing/01E_CLASSIFY_READS_KRAKEN2.py" \
    --trimmed-manifest "$TRIM_MANIFEST" \
    --db "$KRAKEN2_DB" \
    --out-dir "$MAGENTA_KRAKEN_DIR" \
    --manifest-out "$KRAKEN_MANIFEST" \
    --threads-per-sample "$KRAKEN_THREADS_PER_SAMPLE" \
    --max-workers "$KRAKEN_WORKERS"
  RC=$?
fi

if [ "$RC" -eq 0 ]; then
  run_step "01F_SUMMARIZE_KRAKEN2_DIRECT_COUNTS" \
    "$PYTHON_BIN" "$PACKAGE_ROOT/scripts/01_taxonomic_preprocessing/01F_SUMMARIZE_KRAKEN2_DIRECT_OTU_COUNTS.py" \
    --kraken-manifest "$KRAKEN_MANIFEST" \
    --out "$KRAKEN_COUNT_TABLE" \
    --audit-out "$KRAKEN_COUNT_AUDIT"
  RC=$?
fi

if [ "$RC" -eq 0 ]; then
  run_step "01G_BUILD_PHYLOSEQ_405" \
    "$R_BIN" "$PACKAGE_ROOT/scripts/01_taxonomic_preprocessing/01G_BUILD_TAXONOMIC_PHYLOSEQ.R" \
    --kraken-count-table "$KRAKEN_COUNT_TABLE" \
    --metadata "$METADATA_FROZEN" \
    --taxdump-dir "$NCBI_TAXDUMP_DIR" \
    --out-rds "$PHYLOSEQ_405" \
    --audit-dir "$RUN_DIR/01G_phyloseq_audit" \
    --rarefaction-depth "$RAREFACTION_DEPTH" \
    --seed "$RANDOM_SEED" \
    --expected-samples "$EXPECTED_PHYLOSEQ_405_SAMPLES" \
    --expected-taxa "$EXPECTED_PHYLOSEQ_405_TAXA"
  RC=$?
fi

if [ "$RC" -eq 0 ]; then
  run_step "01H_CREATE_SET401" \
    "$R_BIN" "$PACKAGE_ROOT/scripts/01_taxonomic_preprocessing/01H_EXCLUDE_GIRALIA4_CREATE_SET401.R" \
    --base "$BIOGEOGRAPHY_BASE" \
    --input-rds "$PHYLOSEQ_405" \
    --output-rds "$PHYLOSEQ_401"
  RC=$?
fi

if [ "$EXECUTE" -eq 1 ] && [ "$RC" -eq 0 ]; then
  printf "%s\n" "$RUN_DIR" > "$MAGENTA_PIPELINE_RESULTS/LATEST_01_TAXONOMIC_PREPROCESSING.txt"
  printf "PASS\n" > "$RUN_DIR/STATUS.txt"
elif [ "$EXECUTE" -eq 1 ]; then
  printf "FAILED\nexit_code=%s\n" "$RC" > "$RUN_DIR/STATUS.txt"
else
  printf "\nDRY RUN: no se ejecutó ninguna etapa. Use --execute después de revisar rutas y comandos.\n"
fi

printf "\nRUN_DIR=%s\n" "$RUN_DIR"
printf "STATUS_CODE=%s\n" "$RC"
printf "La terminal permanece abierta.\n"
exit "$RC"
