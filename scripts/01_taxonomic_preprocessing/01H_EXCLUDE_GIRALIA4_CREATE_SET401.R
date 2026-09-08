#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# Create the manuscript Set401 phyloseq by excluding the four Giralia samples
# absent from the frozen analysis metadata.
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
    if (index == length(args)) stop(flag, " requiere un valor.", call. = FALSE)
    return(args[[index + 1L]])
  }
  if (length(prefixed) > 0L) {
    return(sub(paste0("^", flag, "="), "", args[[prefixed[[length(prefixed)]]]]))
  }
  default
}

BASE_DIR <- get_arg("--base", "/home/fjbalvino/Biogeografia_2026")

INPUT_RDS <- get_arg("--input-rds", paste0(
  "/data/Magenta/",
  "phyloseq_MAGENTA_EcologicalClass_",
  "MangroveSedimentOnly_TX1_3_OK.rds"
))

OUTPUT_RDS <- get_arg("--output-rds", paste0(
  "/data/Magenta/",
  "phyloseq_MAGENTA.rds"
))

EXCLUDE_ACCESSIONS <- c(
  "SRR29659508",
  "SRR29659599",
  "SRR29659629",
  "SRR29659854"
)

RUN_STAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")

OUT_DIR <- file.path(
  BASE_DIR,
  "resultados",
  paste0(
    "01H_EXCLUDE_GIRALIA4_CREATE_SET401_",
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

msg("Starting canonical Set401 derivation")
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

exclusion_check <- data.frame(
  run_accession = EXCLUDE_ACCESSIONS,
  n_phyloseq_rows = vapply(
    EXCLUDE_ACCESSIONS,
    function(id) sum(metadata$run_accession == id),
    integer(1)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  exclusion_check,
  file.path(
    OUT_DIR,
    "01H_exclusion_accession_check.csv"
  ),
  row.names = FALSE
)

if (any(exclusion_check$n_phyloseq_rows != 1L)) {
  stop_with(
    "Each requested accession must occur exactly once. ",
    "See 01H_exclusion_accession_check.csv."
  )
}

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
    "classified_reads_before",
    "classified_reads_removed",
    "classified_reads_after",
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
    "Canonical MAGENTA Set401 phyloseq",
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
    BASE_DIR,
    "resultados",
    "LATEST_01H_EXCLUDE_GIRALIA4_CREATE_SET401.txt"
  )
)

msg("Samples: ", samples_before, " -> ", samples_after)
msg("Taxa preserved: ", taxa_before, " -> ", taxa_after)
msg("Output MD5: ", output_md5)
msg("Completed successfully: ", OUT_DIR)
