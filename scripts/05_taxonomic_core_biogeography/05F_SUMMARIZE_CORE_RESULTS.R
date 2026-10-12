#!/usr/bin/env Rscript
# Reconcile exported feature/read denominators and summarize Figure 3 statistics.
suppressPackageStartupMessages({library(readr); library(jsonlite)})
results <- Sys.getenv('BIOGEOGRAPHY_RESULTS', unset='results/taxonomic_core')
read_run <- function(tag) readLines(file.path(results,paste0('LATEST_',tag,'.txt')),warn=FALSE)[[1]]
a_run <- read_run('05A_TAXONOMIC_CORE_ASSOCIATED')
b_run <- read_run('05B_REALM_RESTRICTED_TAXA')
c_run <- read_run('05C_CORE_ASSOCIATED_PHYLA')
d_run <- read_run('05D_CORE_NONCORE_DISTANCE_DECAY')
read_table <- function(run,filename) as.data.frame(read_tsv(file.path(run,'tables',filename),show_col_types=FALSE))
features <- read_table(a_run,'09A_all_detected_feature_status.tsv')
components <- read_table(a_run,'09A_Figure_3a_component_summary.tsv')
realms <- read_table(a_run,'09A_realm_locality_counts.tsv')
inclusion <- read_table(a_run,'09A_locality_analysis_sample_inclusion.tsv')
phyla <- read_table(c_run,'11_Figure_3b_full_phylum_composition.tsv')
restriction <- read_table(b_run,'09B_global_restriction_summary.tsv')
restriction_by_realm <- read_table(b_run,'09B_realm_restricted_summary.tsv')
restriction_features <- read_table(b_run,'09B_realm_restricted_taxon_audit_unfiltered.tsv')
statistics <- read_table(d_run,'12_distance_decay_statistics.tsv')
contrasts <- read_table(d_run,'12_core_noncore_slope_contrasts.tsv')
obj <- readRDS(file.path(a_run,'rds','09A_Figure_3a_analysis_objects.rds'))
decay <- readRDS(file.path(d_run,'rds','12_Figure_3cd_analysis_objects.rds'))
core <- features[features$core_associated,,drop=FALSE]
eligible <- features[features$passes_global_eligibility,,drop=FALSE]
noncore <- eligible[!eligible$core_associated,,drop=FALSE]
ineligible <- features[!features$passes_global_eligibility,,drop=FALSE]
stopifnot(all(features$passes_global_eligibility == (features$total_reads >=50 & features$n_localities>=5)),
          all(!features$core_associated | features$passes_global_eligibility),
          sum(features$total_reads)==sum(obj$locality_counts),
          sum(features$total_reads)==18145*sum(inclusion$included_in_locality_analysis),
          sum(phyla$core_reads)==sum(core$total_reads),
          sum(phyla$n_core_features)==nrow(core),
          abs(sum(phyla$locality_mean_percent)-100)<1e-9)
for (i in seq_len(nrow(components))) {
  r <- components[i,]
  subset <- if(r$component=='Core-associated') core else noncore
  numerator <- if(r$metric_code=='richness') nrow(subset) else sum(subset$total_reads)
  denominator <- if(r$metric_code=='richness') nrow(eligible) else sum(eligible$total_reads)
  stopifnot(r$value==numerator,r$denominator==denominator,
            abs(r$percentage-100*numerator/denominator)<1e-9)
}
for (i in seq_len(nrow(restriction))) {
  r <- restriction[i,]
  pool <- if(r$scope=='All detected features') restriction_features else
    restriction_features[restriction_features$passes_primary_taxon_filter,,drop=FALSE]
  restricted <- pool[pool$realm_restricted_unfiltered,,drop=FALSE]
  stopifnot(r$n_features==nrow(pool),r$n_realm_restricted==nrow(restricted),
            r$total_reads==sum(pool$total_reads),r$realm_restricted_reads==sum(restricted$total_reads))
}
report <- list(
  status='PASS',
  cohort=list(n_reference_profiles=nrow(inclusion),
              n_locality_analysis_profiles=sum(inclusion$included_in_locality_analysis),
              excluded_from_locality_analysis=inclusion[!inclusion$included_in_locality_analysis,,drop=FALSE],
              n_realm_locality_units=sum(realms$n_localities),realm_units=realms),
  definitions=list(global_eligibility='At least 50 retained reads and detection in at least five realm-locality units.',
                   core='Detection in every sampled locality of at least one realm represented by at least three localities, with mean within-locality relative abundance >=0.0001 (0.01%) in that realm.',
                   eligible_read_denominator='Retained reads assigned to features passing global eligibility.',
                   locality_balanced_phyla='Pool counts by realm-locality unit, calculate each phylum as a proportion of its core-associated reads, then average with equal weight across units.',
                   realm_restriction='Detected in one sampled realm; same locality-resolved cohort as the core analysis, without global eligibility filtering.'),
  counts=list(n_detected_features=nrow(features),n_eligible_features=nrow(eligible),
              n_core_features=nrow(core),n_eligible_noncore_features=nrow(noncore),
              n_ineligible_features=nrow(ineligible),reads_all_retained=sum(features$total_reads),
              reads_eligible=sum(eligible$total_reads),reads_core=sum(core$total_reads),
              reads_eligible_noncore=sum(noncore$total_reads),reads_ineligible=sum(ineligible$total_reads),
              core_percent_of_eligible_richness=100*nrow(core)/nrow(eligible),
              core_percent_of_eligible_reads=100*sum(core$total_reads)/sum(eligible$total_reads),
              core_percent_of_all_retained_reads=100*sum(core$total_reads)/sum(features$total_reads)),
  n_core_phyla=nrow(phyla),core_phyla=phyla,
  realm_restriction=restriction,realm_restriction_by_realm=restriction_by_realm,
  distance_decay=list(statistics=statistics,slope_contrasts=contrasts,
                      geographic_mean_core_similarity=mean(decay$geographic$result$pairs$core_similarity),
                      geographic_mean_noncore_similarity=mean(decay$geographic$result$pairs$noncore_similarity),
                      climatic_mean_core_similarity=mean(decay$environmental$result$pairs$core_similarity),
                      climatic_mean_noncore_similarity=mean(decay$environmental$result$pairs$noncore_similarity),
                      permutations=9999),
  verification=c('Input cohort and locality exclusions explicitly reconciled',
                 'Read conservation at 18,145 retained reads per analyzed profile',
                 'Feature/read partition and global eligibility independently recalculated',
                 'All figure component numerators, denominators and percentages reconciled',
                 'Phylum counts and equal-locality percentages reconciled',
                 'Realm-restricted feature/read totals independently recalculated'))
write_json(report,file.path(results,'core_results_summary.json'),pretty=TRUE,auto_unbox=TRUE,digits=NA,na='null')
cat('CORE_SUMMARY_VERIFICATION_PASS\n')
