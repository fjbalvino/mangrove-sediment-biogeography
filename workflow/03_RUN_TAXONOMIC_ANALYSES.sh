#!/usr/bin/env bash
# One reference cohort and one current output directory per analysis.
set -euo pipefail
PACKAGE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
EXECUTE=0
WITH_CORE=0
CONFIG=""
while (($#)); do
  case "$1" in
    --execute) EXECUTE=1; shift ;;
    --core) WITH_CORE=1; shift ;;
    --config)
      if (($# < 2)); then echo "--config requires a file" >&2; exit 2; fi
      CONFIG="$2"; shift 2 ;;
    --assembly)
      echo "iCAMP retains its separate completed analysis; this workflow does not rerun it." >&2
      exit 2 ;;
    --help|-h)
      echo "Usage: $0 [--config FILE] [--execute] [--core]"
      echo "Repeat diversity, GAM, composition and optionally core analyses."
      exit 0 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done
source "$PACKAGE_ROOT/workflow/lib/config.sh"
magenta_load_config "$CONFIG"
export BIOGEOGRAPHY_BASE="${BIOGEOGRAPHY_BASE:-$PACKAGE_ROOT}"
export BIOGEOGRAPHY_RESULTS="${BIOGEOGRAPHY_RESULTS:-$BIOGEOGRAPHY_BASE/results}"
export MAGENTA_PHYLOSEQ="${MAGENTA_PHYLOSEQ:-$BIOGEOGRAPHY_BASE/data/phyloseq_taxonomic.rds}"
export MAGENTA_CORE_PHYLOSEQ="$MAGENTA_PHYLOSEQ"
export MAGENTA_METADATA_CSV="${MAGENTA_METADATA_CSV:-$BIOGEOGRAPHY_BASE/data/metadata_taxonomic.csv}"
export ALPHA_SAMPLE_MANIFEST="${ALPHA_SAMPLE_MANIFEST:-$PACKAGE_ROOT/config/samples/taxonomic_accessions.txt}"
export CORE_SAMPLE_MANIFEST="${CORE_SAMPLE_MANIFEST:-$ALPHA_SAMPLE_MANIFEST}"
R_BIN="${R_BIN:-Rscript}"
run_r() {
  printf '%q ' "$R_BIN" "$PACKAGE_ROOT/$1"
  printf '\n'
  if ((EXECUTE)); then "$R_BIN" "$PACKAGE_ROOT/$1"; fi
}
run_r scripts/03_alpha_diversity/03A_CALCULATE_HILL_DIVERSITY_BY_REALM.R
run_r scripts/03_alpha_diversity/03E_HILL_DIVERSITY_SITE_INFERENCE.R
run_r scripts/03_alpha_diversity/03F_TAXONOMIC_DOMINANCE.R
run_r scripts/03_alpha_diversity/03B_FIT_SHANNON_ENVIRONMENTAL_GAM.R
run_r scripts/03_alpha_diversity/03C_PLOT_TOP4_SHANNON_GAM_EFFECTS.R
run_r scripts/03_alpha_diversity/03D_PLOT_REMAINING_SIGNIFICANT_GAM_EFFECTS.R
run_r scripts/04_beta_diversity/04A_NMDS_AITCHISON_LONGITUDE.R
run_r scripts/04_beta_diversity/04B_DISTANCE_DECAY_AND_MANTEL.R
run_r scripts/04_beta_diversity/04C_PERMANOVA_REGIONAL_STRUCTURE.R
run_r scripts/04_beta_diversity/04D_DBRDA_VARIANCE_PARTITIONING.R
if ((WITH_CORE)); then
  (
    export BIOGEOGRAPHY_RESULTS="$BIOGEOGRAPHY_RESULTS/core"
    for script in 05A_CLASSIFY_TAXONOMIC_CORE_ASSOCIATED 05B_CLASSIFY_REALM_RESTRICTED_TAXA 05C_SUMMARIZE_CORE_ASSOCIATED_PHYLA 05D_CORE_NONCORE_DISTANCE_DECAY 05E_COMPOSE_FIGURE_3 05F_SUMMARIZE_CORE_RESULTS; do
      run_r "scripts/05_taxonomic_core_biogeography/${script}.R"
    done
  )
fi

export ICAMP_SUMMARY_TSV="${ICAMP_SUMMARY_TSV:-$MAGENTA_DATA_DIR/reference_results/icamp/icamp_process_contributions.tsv}"
run_r scripts/03_alpha_diversity/03G_COMPOSE_TAXONOMIC_FIGURE.R
