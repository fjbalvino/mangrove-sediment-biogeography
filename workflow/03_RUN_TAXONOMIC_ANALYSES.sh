#!/usr/bin/env bash
# Run the taxonomic analyses in dependency order. Commands are printed by default.
set -euo pipefail

PACKAGE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
EXECUTE=0
WITH_CORE=0
WITH_ASSEMBLY=0
CONFIG=""
while (($#)); do
  case "$1" in
    --execute) EXECUTE=1; shift ;;
    --core) WITH_CORE=1; shift ;;
    --assembly) WITH_ASSEMBLY=1; shift ;;
    --config)
      if (($# < 2)); then echo "--config requires a file" >&2; exit 2; fi
      CONFIG="$2"; shift 2 ;;
    --help|-h)
      echo "Usage: $0 [--config FILE] [--execute] [--core] [--assembly]"
      echo "Runs diversity, GAM and beta-diversity analyses. --core adds core analyses;"
      echo "--assembly adds iCAMP and its sampling-robustness analysis."
      exit 0 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done
if [[ -n "$CONFIG" ]]; then
  if [[ ! -f "$CONFIG" ]]; then echo "Configuration file not found: $CONFIG" >&2; exit 2; fi
  source "$CONFIG"
fi
export BIOGEOGRAPHY_BASE="${BIOGEOGRAPHY_BASE:-$PACKAGE_ROOT}"
export BIOGEOGRAPHY_RESULTS="${BIOGEOGRAPHY_RESULTS:-$BIOGEOGRAPHY_BASE/results}"
export MAGENTA_PHYLOSEQ_405="${MAGENTA_PHYLOSEQ_405:-$BIOGEOGRAPHY_BASE/data/taxonomic_counts.rds}"
export MAGENTA_PHYLOSEQ_401="${MAGENTA_PHYLOSEQ_401:-$BIOGEOGRAPHY_BASE/data/taxonomic_core_counts.rds}"
export MAGENTA_METADATA_CSV="${MAGENTA_METADATA_CSV:-$BIOGEOGRAPHY_BASE/data/environmental_metadata.csv}"
export ALPHA_SAMPLE_MANIFEST="${ALPHA_SAMPLE_MANIFEST:-$PACKAGE_ROOT/config/samples/taxonomic_diversity_accessions.txt}"
export CORE_SAMPLE_MANIFEST="${CORE_SAMPLE_MANIFEST:-$PACKAGE_ROOT/config/samples/taxonomic_core_accessions.txt}"
R_BIN="${R_BIN:-Rscript}"

run_r() {
  printf '%q ' "$R_BIN" "$PACKAGE_ROOT/$1"
  printf '\n'
  if ((EXECUTE)); then "$R_BIN" "$PACKAGE_ROOT/$1"; fi
}

# The 401-sample object is a sample subset; count values are not transformed here.
if [[ ! -f "$MAGENTA_PHYLOSEQ_401" ]]; then
  run_r scripts/01_taxonomic_preprocessing/01H_SELECT_CORE_ANALYSIS_SAMPLES.R
fi
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
  run_r scripts/05_taxonomic_core_biogeography/05A_CLASSIFY_TAXONOMIC_CORE_ASSOCIATED.R
  run_r scripts/05_taxonomic_core_biogeography/05B_CLASSIFY_REALM_RESTRICTED_TAXA.R
  run_r scripts/05_taxonomic_core_biogeography/05C_SUMMARIZE_CORE_ASSOCIATED_PHYLA.R
  run_r scripts/05_taxonomic_core_biogeography/05D_CORE_NONCORE_DISTANCE_DECAY.R
fi
if ((WITH_ASSEMBLY)); then
  run_r scripts/06_community_assembly/06A_ICAMP_GENUS_COMMUNITY_ASSEMBLY.R
  run_r scripts/06_community_assembly/06B_ICAMP_ROBUSTNESS_REALM_STUDY_LOCALITY.R
fi
