#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# Select the samples used for core and community-assembly analyses.
# Membership is defined by the accession manifest supplied with the repository.
# Only samples are removed; the original taxonomic feature space is preserved.

options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(phyloseq)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default) {
  exact <- which(args == flag)
  prefixed <- grep(paste0("^", flag, "="), args)
  if (length(exact) > 0L) {
    index <- exact[[length(exact)]]
    if (index == length(args)) stop(flag, " requires a value.", call. = FALSE)
    return(args[[index + 1L]])
  }
  if (length(prefixed) > 0L) {
    return(sub(paste0("^", flag, "="), "", args[[prefixed[[length(prefixed)]]]]))
  }
  default
}

BASE_DIR <- get_arg("--base", Sys.getenv("BIOGEOGRAPHY_BASE", unset = normalizePath(".")))
RESULTS_DIR <- Sys.getenv("BIOGEOGRAPHY_RESULTS", unset = file.path(BASE_DIR, "results"))
INPUT_RDS <- get_arg("--input-rds", Sys.getenv("MAGENTA_PHYLOSEQ_405", unset = file.path(BASE_DIR, "data", "taxonomic_counts.rds")))
OUTPUT_RDS <- get_arg("--output-rds", Sys.getenv("MAGENTA_PHYLOSEQ_401", unset = file.path(BASE_DIR, "data", "taxonomic_core_counts.rds")))
SAMPLE_MANIFEST <- get_arg("--manifest", Sys.getenv("CORE_SAMPLE_MANIFEST", unset = file.path(BASE_DIR, "config", "samples", "taxonomic_core_accessions.txt")))
if (!file.exists(SAMPLE_MANIFEST)) stop("Sample manifest not found: ", SAMPLE_MANIFEST, call. = FALSE)
KEEP_ACCESSIONS <- trimws(readLines(SAMPLE_MANIFEST, warn = FALSE))
KEEP_ACCESSIONS <- KEEP_ACCESSIONS[nzchar(KEEP_ACCESSIONS) & !grepl("^#", KEEP_ACCESSIONS)]
if (length(KEEP_ACCESSIONS) != 401L || anyDuplicated(KEEP_ACCESSIONS)) {
  stop("The core analysis manifest must contain 401 unique accessions.", call. = FALSE)
}

RUN_STAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")

OUT_DIR <- file.path(
  RESULTS_DIR,
  paste0(
    "01H_CORE_ANALYSIS_SAMPLES_",
    RUN_STAMP
  )
)

msg <- function(...) {
  cat(
    sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0(..., collapse = ""),
    "\n",
    sep = ""
  )
}

stop_with <- function(...) {
  stop(paste0(..., collapse = ""), call. = FALSE)
}

normalize_text <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "null", "NULL")] <- NA_character_
  x
}

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

msg("Selecting the core analysis samples")
msg("Input:  ", INPUT_RDS)
msg("Output: ", OUTPUT_RDS)
msg("Audit:  ", OUT_DIR)

if (!file.exists(INPUT_RDS)) {
  stop_with("Input phyloseq object not found: ", INPUT_RDS)
}

if (file.exists(OUTPUT_RDS)) {
  stop_with(
    "The derived object already exists and will not be overwritten: ",
    OUTPUT_RDS
  )
}

ps <- readRDS(INPUT_RDS)

if (!inherits(ps, "phyloseq")) {
  stop_with("The input RDS is not a phyloseq object.")
}

metadata <- as.data.frame(
  as(sample_data(ps), "data.frame")
)

metadata$.ps_sample_id <- rownames(metadata)

if (!"run_accession" %in% names(metadata)) {
  stop_with("sample_data(ps) does not contain run_accession.")
}

metadata$run_accession <- normalize_text(
  metadata$run_accession
)

invalid_ids <- metadata[
  is.na(metadata$run_accession) |
    duplicated(metadata$run_accession) |
    duplicated(metadata$run_accession, fromLast = TRUE),
  c(".ps_sample_id", "run_accession"),
  drop = FALSE
]

write.csv(
  invalid_ids,
  file.path(
    OUT_DIR,
    "01H_invalid_phyloseq_run_accessions.csv"
  ),
  row.names = FALSE
)

if (nrow(invalid_ids) > 0L) {
  stop_with(
    "Missing or duplicated run_accession values were found."
  )
}

if (nsamples(ps) != 405L) stop_with("Expected 405 samples in the source phyloseq object.")
missing_ids <- setdiff(KEEP_ACCESSIONS, metadata$run_accession)
if (length(missing_ids)) stop_with("Manifest accessions missing from phyloseq: ", paste(missing_ids, collapse = ", "))
EXCLUDE_ACCESSIONS <- setdiff(metadata$run_accession, KEEP_ACCESSIONS)
write.csv(data.frame(run_accession = metadata$run_accession,
                     included = metadata$run_accession %in% KEEP_ACCESSIONS),
          file.path(OUT_DIR, "01H_sample_membership.csv"), row.names = FALSE)

removed_metadata <- metadata[
  metadata$run_accession %in% EXCLUDE_ACCESSIONS,
  ,
  drop = FALSE
]

removed_metadata <- removed_metadata[
  match(
    EXCLUDE_ACCESSIONS,
    removed_metadata$run_accession
  ),
  ,
  drop = FALSE
]

write.csv(
  removed_metadata,
  file.path(
    OUT_DIR,
    "01H_removed_samples_full_metadata.csv"
  ),
  row.names = FALSE
)

samples_before <- nsamples(ps)
taxa_before <- ntaxa(ps)
reads_before <- sum(sample_sums(ps))

keep_sample_names <- metadata$.ps_sample_id[
  !metadata$run_accession %in% EXCLUDE_ACCESSIONS
]

ps_derived <- prune_samples(
  keep_sample_names,
  ps
)

# Do not prune taxa here. Downstream scripts will remove zero-sum taxa
# after their final cohort filters.
samples_after <- nsamples(ps_derived)
taxa_after <- ntaxa(ps_derived)
reads_after <- sum(sample_sums(ps_derived))

derived_metadata <- as.data.frame(
  as(sample_data(ps_derived), "data.frame")
)

derived_accessions <- normalize_text(
  derived_metadata$run_accession
)

if (samples_after != samples_before - length(EXCLUDE_ACCESSIONS)) {
  stop_with("Unexpected sample count after exclusion.")
}

if (taxa_after != taxa_before) {
  stop_with(
    "The feature space changed although only samples should be removed."
  )
}

if (any(EXCLUDE_ACCESSIONS %in% derived_accessions)) {
  stop_with(
    "At least one excluded accession remains in the derived object."
  )
}

dir.create(dirname(OUTPUT_RDS), recursive = TRUE, showWarnings = FALSE)
saveRDS(
  ps_derived,
  OUTPUT_RDS,
  compress = "xz"
)

# Independent post-save validation
ps_check <- readRDS(OUTPUT_RDS)

check_metadata <- as.data.frame(
  as(sample_data(ps_check), "data.frame")
)

check_accessions <- normalize_text(
  check_metadata$run_accession
)

if (
  nsamples(ps_check) != samples_after ||
  ntaxa(ps_check) != taxa_after ||
  any(EXCLUDE_ACCESSIONS %in% check_accessions)
) {
  stop_with("Post-save validation of the derived RDS failed.")
}

output_md5 <- unname(
  tools::md5sum(OUTPUT_RDS)
)

summary_table <- data.frame(
  metric = c(
    "input_rds",
    "output_rds",
    "samples_before",
    "samples_removed",
    "samples_after",
    "taxa_before",
    "taxa_after",
    "retained_reads_before",
    "retained_reads_removed",
    "retained_reads_after",
    "output_md5"
  ),
  value = c(
    INPUT_RDS,
    OUTPUT_RDS,
    samples_before,
    length(EXCLUDE_ACCESSIONS),
    samples_after,
    taxa_before,
    taxa_after,
    format(reads_before, scientific = FALSE),
    format(reads_before - reads_after, scientific = FALSE),
    format(reads_after, scientific = FALSE),
    output_md5
  ),
  stringsAsFactors = FALSE
)

write.csv(
  summary_table,
  file.path(
    OUT_DIR,
    "01H_run_summary.csv"
  ),
  row.names = FALSE
)

writeLines(
  c(
    "Canonical MAGENTA core analysis phyloseq",
    paste0(
      "Timestamp: ",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    ),
    paste0("Input RDS: ", INPUT_RDS),
    paste0("Output RDS: ", OUTPUT_RDS),
    paste0(
      "Samples: ",
      samples_before,
      " -> ",
      samples_after
    ),
    paste0(
      "Taxa retained unchanged: ",
      taxa_after
    ),
    paste0(
      "Excluded accessions: ",
      paste(EXCLUDE_ACCESSIONS, collapse = ", ")
    ),
    paste0("Output MD5: ", output_md5),
    paste0("Audit directory: ", OUT_DIR)
  ),
  file.path(
    OUT_DIR,
    "01H_run_summary.txt"
  )
)

writeLines(
  OUT_DIR,
  file.path(
    RESULTS_DIR,
    "LATEST_CORE_ANALYSIS_SAMPLES.txt"
  )
)

msg("Samples: ", samples_before, " -> ", samples_after)
msg("Taxa preserved: ", taxa_before, " -> ", taxa_after)
msg("Output MD5: ", output_md5)
msg("Completed successfully: ", OUT_DIR)
