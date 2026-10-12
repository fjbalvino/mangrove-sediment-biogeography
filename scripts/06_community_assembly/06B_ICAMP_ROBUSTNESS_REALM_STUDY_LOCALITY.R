#!/usr/bin/env Rscript

# =============================================================================
# 06B_ICAMP_ROBUSTNESS_REALM_STUDY_LOCALITY.R
#
# Robustness analysis for the completed genus-level iCAMP run.
#
# Primary claim tested
# --------------------
# Test whether the within-realm fraction not assigned to heterogeneous or
# homogeneous selection remains >0.5 after equal weighting of marine realms,
# studies and localities. The sum uses these output columns:
#
#   dispersal_limitation + homogenizing_dispersal + drift_weak_selection
#
# The last category includes drift, weak selection and other undominated
# processes; this sum must not be interpreted as a purely stochastic fraction.
#
# The script never treats the 24,692 pairwise turnovers as independent
# replicates. Inference is based on equal-weight group summaries, group-level
# bootstrap, sign-flip tests, leave-one-group-out analyses, and balanced
# post-hoc resampling of samples among adequately represented realms.
# =============================================================================

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)

SCRIPT_ID <- "08I_ICAMP_ROBUSTNESS_REALM_STUDY_LOCALITY"
SCRIPT_VERSION <- "1.1.0"

PROCESS_COLUMNS <- c(
  "heterogeneous_selection",
  "homogeneous_selection",
  "dispersal_limitation",
  "homogenizing_dispersal",
  "drift_weak_selection"
)

SELECTIVE_COLUMNS <- c(
  "heterogeneous_selection",
  "homogeneous_selection"
)

NONSELECTIVE_COLUMNS <- c(
  "dispersal_limitation",
  "homogenizing_dispersal",
  "drift_weak_selection"
)

get_arg <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  direct <- which(args == flag)
  if (length(direct) > 0L) {
    i <- direct[length(direct)]
    if (i < length(args)) return(args[i + 1L])
    return(default)
  }
  prefix <- paste0(flag, "=")
  inline <- args[startsWith(args, prefix)]
  if (length(inline) > 0L) {
    return(sub(prefix, "", inline[length(inline)], fixed = TRUE))
  }
  default
}

as_flag <- function(x) {
  tolower(trimws(as.character(x))) %in%
    c("1", "true", "t", "yes", "y", "si", "sí")
}

clean_chr <- function(x) {
  y <- trimws(as.character(x))
  bad <- is.na(y) | y == "" |
    tolower(y) %in% c("na", "n/a", "nan", "null", "none", "unknown")
  y[bad] <- NA_character_
  y
}

normalise_name <- function(x) {
  tolower(gsub("[^A-Za-z0-9]+", "", as.character(x)))
}

find_column <- function(x, candidates, label) {
  direct <- candidates[candidates %in% names(x)]
  if (length(direct) > 0L) return(direct[1L])
  idx <- match(normalise_name(candidates), normalise_name(names(x)))
  idx <- idx[!is.na(idx)]
  if (length(idx) > 0L) return(names(x)[idx[1L]])
  stop(
    sprintf(
      "Could not identify %s. Tried: %s",
      label,
      paste(candidates, collapse = ", ")
    ),
    call. = FALSE
  )
}

read_tsv <- function(path) {
  utils::read.delim(
    path,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

write_tsv <- function(x, path) {
  utils::write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = ""
  )
}

write_tsv_gz <- function(x, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  utils::write.table(
    x,
    file = con,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = ""
  )
}

write_status <- function(outdir, status, message = NULL) {
  lines <- status
  if (!is.null(message) && nzchar(message)) {
    lines <- c(lines, paste0("message=", message))
  }
  lines <- c(
    lines,
    paste0(
      "timestamp_utc=",
      format(Sys.time(), tz = "UTC", usetz = TRUE)
    )
  )
  writeLines(lines, file.path(outdir, "STATUS.txt"))
}

as_bool <- function(x) {
  toupper(trimws(as.character(x))) %in% c("TRUE", "T", "1", "YES")
}

long_value <- function(x, key, key_column = names(x)[1L]) {
  value <- x[[setdiff(names(x), key_column)[1L]]][x[[key_column]] == key]
  if (length(value) == 0L) return(NA_character_)
  as.character(value[1L])
}

sha256_file <- function(path) {
  executable <- Sys.which("sha256sum")
  if (!nzchar(executable) || !file.exists(path)) return(NA_character_)
  output <- suppressWarnings(system2(executable, shQuote(path), stdout = TRUE))
  if (length(output) == 0L) return(NA_character_)
  value <- strsplit(output[1L], "[[:space:]]+")[[1L]][1L]
  if (!grepl("^[0-9a-fA-F]{64}$", value)) NA_character_ else tolower(value)
}

add_fractions <- function(x) {
  x$selective_fraction <- rowSums(x[, SELECTIVE_COLUMNS, drop = FALSE])
  x$nonselective_fraction <- rowSums(
    x[, NONSELECTIVE_COLUMNS, drop = FALSE]
  )
  x
}

bootstrap_group_matrix <- function(group_matrix, n_boot, seed) {
  set.seed(seed)
  n_groups <- nrow(group_matrix)
  if (n_groups < 2L) {
    stop("At least two eligible groups are required for bootstrap.",
         call. = FALSE)
  }
  boot <- replicate(n_boot, {
    index <- sample.int(n_groups, n_groups, replace = TRUE)
    colMeans(group_matrix[index, , drop = FALSE])
  })
  if (is.null(dim(boot))) boot <- matrix(boot, nrow = ncol(group_matrix))
  data.frame(
    metric = colnames(group_matrix),
    estimate = colMeans(group_matrix),
    bootstrap_lwr = apply(boot, 1L, stats::quantile, probs = 0.025,
                          names = FALSE),
    bootstrap_upr = apply(boot, 1L, stats::quantile, probs = 0.975,
                          names = FALSE),
    n_groups = n_groups,
    n_boot = n_boot,
    stringsAsFactors = FALSE
  )
}

sign_flip_test <- function(values, null = 0.5, n_perm = 9999L, seed = 42L) {
  values <- values[is.finite(values)]
  differences <- values - null
  n <- length(differences)
  observed <- mean(differences)
  if (n < 2L) {
    return(data.frame(
      estimate = mean(values), null = null, difference = observed,
      p_one_sided = NA_real_, n_groups = n, method = "insufficient_groups",
      stringsAsFactors = FALSE
    ))
  }

  if (n <= 15L) {
    signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), n)))
    permuted <- as.vector(signs %*% differences / n)
    p_value <- mean(permuted >= observed - 1e-15)
    method <- paste0("exact_sign_flip_", 2^n, "_configurations")
  } else {
    set.seed(seed)
    permuted <- replicate(
      n_perm,
      mean(sample(c(-1, 1), n, replace = TRUE) * differences)
    )
    p_value <- (1 + sum(permuted >= observed - 1e-15)) / (n_perm + 1)
    method <- paste0("monte_carlo_sign_flip_", n_perm)
  }

  data.frame(
    estimate = mean(values),
    null = null,
    difference = observed,
    p_one_sided = p_value,
    n_groups = n,
    method = method,
    stringsAsFactors = FALSE
  )
}

build_group_sizes <- function(metadata, unit, group_column,
                              minimum_samples) {
  counts <- as.data.frame(
    table(metadata[[group_column]], useNA = "no"),
    stringsAsFactors = FALSE
  )
  names(counts) <- c("group", "n_samples")
  counts$n_samples <- as.integer(counts$n_samples)
  counts$unit <- unit
  counts$minimum_samples <- minimum_samples
  counts$sample_size_eligible <- counts$n_samples >= minimum_samples
  counts[, c(
    "unit", "group", "n_samples", "minimum_samples",
    "sample_size_eligible"
  )]
}

build_group_means <- function(pairwise, group_sizes, unit,
                              endpoint1, endpoint2,
                              minimum_pairs = 1L) {
  same <- !is.na(pairwise[[endpoint1]]) & !is.na(pairwise[[endpoint2]]) &
    pairwise[[endpoint1]] == pairwise[[endpoint2]]
  x <- pairwise[same, , drop = FALSE]
  x$.analysis_group <- x[[endpoint1]]
  split_rows <- split(seq_len(nrow(x)), x$.analysis_group, drop = TRUE)

  rows <- lapply(names(split_rows), function(group) {
    index <- split_rows[[group]]
    values <- x[index, , drop = FALSE]
    means <- colMeans(
      values[, c(PROCESS_COLUMNS, "selective_fraction",
                 "nonselective_fraction"), drop = FALSE]
    )
    data.frame(
      unit = unit,
      group = group,
      n_pairs = nrow(values),
      t(means),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })

  if (length(rows) == 0L) {
    stop(sprintf("No within-group pairs were available for unit %s.", unit),
         call. = FALSE)
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out <- merge(
    out,
    group_sizes[, c("group", "n_samples", "minimum_samples",
                    "sample_size_eligible")],
    by = "group",
    all.x = TRUE,
    sort = FALSE
  )
  out$minimum_pairs <- minimum_pairs
  out$eligible <- out$sample_size_eligible & out$n_pairs >= minimum_pairs
  out[, c(
    "unit", "group", "n_samples", "n_pairs", "minimum_samples",
    "minimum_pairs", "eligible", PROCESS_COLUMNS,
    "selective_fraction", "nonselective_fraction"
  )]
}

equal_group_summary <- function(group_means, n_boot, n_perm, seed) {
  eligible <- group_means[group_means$eligible, , drop = FALSE]
  metrics <- c(PROCESS_COLUMNS, "selective_fraction", "nonselective_fraction")
  matrix_values <- as.matrix(eligible[, metrics, drop = FALSE])
  storage.mode(matrix_values) <- "double"
  bootstrap <- bootstrap_group_matrix(matrix_values, n_boot, seed)
  bootstrap$unit <- unique(group_means$unit)
  bootstrap$weighting <- "equal_weight_per_eligible_group"

  sign_test <- sign_flip_test(
    eligible$nonselective_fraction,
    null = 0.5,
    n_perm = n_perm,
    seed = seed + 1000L
  )
  sign_test$unit <- unique(group_means$unit)
  sign_test$weighting <- "equal_weight_per_eligible_group"
  list(process = bootstrap, nonselective = sign_test)
}

equal_group_loo <- function(group_means) {
  eligible <- group_means[group_means$eligible, , drop = FALSE]
  metrics <- c(PROCESS_COLUMNS, "selective_fraction", "nonselective_fraction")
  if (nrow(eligible) < 3L) {
    stop(sprintf("Unit %s has fewer than three eligible groups for LOO.",
                 unique(group_means$unit)), call. = FALSE)
  }
  rows <- lapply(seq_len(nrow(eligible)), function(i) {
    remaining <- eligible[-i, , drop = FALSE]
    means <- colMeans(remaining[, metrics, drop = FALSE])
    data.frame(
      unit = unique(group_means$unit),
      loo_type = "equal_group_mean",
      omitted_group = eligible$group[i],
      omitted_group_n_samples = eligible$n_samples[i],
      omitted_group_n_pairs = eligible$n_pairs[i],
      n_groups_remaining = nrow(remaining),
      n_pairs_remaining = sum(remaining$n_pairs),
      t(means),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  do.call(rbind, rows)
}

endpoint_removal_loo <- function(pairwise, group_means,
                                 endpoint1, endpoint2) {
  eligible_groups <- group_means$group[group_means$eligible]
  metrics <- c(PROCESS_COLUMNS, "selective_fraction", "nonselective_fraction")
  rows <- lapply(eligible_groups, function(group) {
    keep <- pairwise[[endpoint1]] != group & pairwise[[endpoint2]] != group
    remaining <- pairwise[keep, , drop = FALSE]
    means <- colMeans(remaining[, metrics, drop = FALSE])
    group_row <- group_means[group_means$group == group, , drop = FALSE][1L, ]
    data.frame(
      unit = unique(group_means$unit),
      loo_type = "remove_all_endpoints_in_group",
      omitted_group = group,
      omitted_group_n_samples = group_row$n_samples,
      omitted_group_n_pairs = group_row$n_pairs,
      n_groups_remaining = NA_integer_,
      n_pairs_remaining = nrow(remaining),
      t(means),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  do.call(rbind, rows)
}

balanced_realm_resampling <- function(pairwise, metadata,
                                      realm_column,
                                      minimum_realm_samples,
                                      n_replicates, seed) {
  counts <- table(metadata[[realm_column]])
  eligible_realms <- names(counts)[counts >= minimum_realm_samples]
  if (length(eligible_realms) < 3L) {
    stop("Fewer than three realms meet the balanced-resampling threshold.",
         call. = FALSE)
  }
  target_n <- min(as.integer(counts[eligible_realms]))
  ids_by_realm <- split(metadata$run_accession, metadata[[realm_column]])
  pairs_by_realm <- split(pairwise, pairwise$realm1)

  set.seed(seed)
  result <- lapply(seq_len(n_replicates), function(replicate_id) {
    realm_values <- vapply(eligible_realms, function(realm) {
      selected <- sample(ids_by_realm[[realm]], target_n, replace = FALSE)
      x <- pairs_by_realm[[realm]]
      keep <- x$sample1 %in% selected & x$sample2 %in% selected
      if (sum(keep) != choose(target_n, 2)) {
        stop(sprintf(
          "Balanced subsample for %s produced %d pairs; expected %d.",
          realm, sum(keep), choose(target_n, 2)
        ), call. = FALSE)
      }
      mean(x$nonselective_fraction[keep])
    }, numeric(1))
    data.frame(
      replicate = replicate_id,
      n_realms = length(eligible_realms),
      samples_per_realm = target_n,
      pairs_per_realm = choose(target_n, 2),
      equal_realm_nonselective_fraction = mean(realm_values),
      minimum_realm_nonselective_fraction = min(realm_values),
      maximum_realm_nonselective_fraction = max(realm_values),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, result)
  attr(out, "eligible_realms") <- eligible_realms
  out
}

print_help <- function() {
  cat(paste(
    "Usage:",
    "  Rscript scripts/06_community_assembly/06B_ICAMP_ROBUSTNESS_REALM_STUDY_LOCALITY.R [options]",
    "",
    "Options:",
    "  --run PATH                 Completed 06A iCAMP output directory",
    "  --base PATH                Default BIOGEOGRAPHY_BASE or current directory",
    "  --results_dir PATH         Default BIOGEOGRAPHY_RESULTS or BASE/results",
    "  --n_boot INT               Default 10000",
    "  --n_perm INT               Default 9999",
    "  --n_balance INT            Default 2000",
    "  --min_realm_samples INT    Default 10",
    "  --min_study_samples INT    Default 3",
    "  --min_locality_samples INT Default 3",
    "  --seed INT                 Default 42",
    sep = "\n"
  ))
}

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args || "-h" %in% args) {
  print_help()
  quit(save = "no", status = 0L)
}

PROJECT_DIR <- get_arg("--base", Sys.getenv("BIOGEOGRAPHY_BASE", "."))
BASE_DIR <- PROJECT_DIR
RESULTS_DIR <- get_arg(
  "--results_dir", Sys.getenv("BIOGEOGRAPHY_RESULTS", file.path(PROJECT_DIR, "results"))
)
SOURCE_LABEL <- "08B_ICAMP_GENUS_COMMUNITY_ASSEMBLY"
LATEST_SOURCE <- file.path(RESULTS_DIR, paste0("LATEST_", SOURCE_LABEL, ".txt"))

SOURCE_RUN <- get_arg("--run", Sys.getenv("ICAMP_RUN", ""))
if (!nzchar(SOURCE_RUN)) {
  if (!file.exists(LATEST_SOURCE)) {
    stop("Source run was not supplied and latest pointer does not exist: ",
         LATEST_SOURCE, call. = FALSE)
  }
  SOURCE_RUN <- trimws(readLines(LATEST_SOURCE, warn = FALSE)[1L])
}

N_BOOT <- as.integer(get_arg("--n_boot", "10000"))
N_PERM <- as.integer(get_arg("--n_perm", "9999"))
N_BALANCE <- as.integer(get_arg("--n_balance", "2000"))
MIN_REALM_SAMPLES <- as.integer(get_arg("--min_realm_samples", "10"))
MIN_STUDY_SAMPLES <- as.integer(get_arg("--min_study_samples", "3"))
MIN_LOCALITY_SAMPLES <- as.integer(get_arg("--min_locality_samples", "3"))
SEED <- as.integer(get_arg("--seed", "42"))

numeric_parameters <- c(
  N_BOOT, N_PERM, N_BALANCE, MIN_REALM_SAMPLES,
  MIN_STUDY_SAMPLES, MIN_LOCALITY_SAMPLES, SEED
)
if (any(!is.finite(numeric_parameters)) || N_BOOT < 1000L ||
    N_PERM < 999L || N_BALANCE < 100L ||
    min(MIN_REALM_SAMPLES, MIN_STUDY_SAMPLES,
        MIN_LOCALITY_SAMPLES) < 2L) {
  stop("Invalid numerical parameters.", call. = FALSE)
}

STAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUT_DIR <- file.path(RESULTS_DIR, paste0(SCRIPT_ID, "_", STAMP))
TABLES_DIR <- file.path(OUT_DIR, "tables")
LOGS_DIR <- file.path(OUT_DIR, "logs")
dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOGS_DIR, recursive = TRUE, showWarnings = FALSE)

LATEST_OUT <- file.path(RESULTS_DIR, paste0("LATEST_", SCRIPT_ID, ".txt"))
writeLines(OUT_DIR, LATEST_OUT)
write_status(OUT_DIR, "RUNNING")

log_connection <- file(file.path(LOGS_DIR, "08I_log.txt"), open = "wt")
sink(log_connection, split = TRUE)
sink(log_connection, type = "message")

status_ok <- FALSE

main <- function() {
  cat("============================================================\n")
  cat(SCRIPT_ID, "\n")
  cat("Script version: ", SCRIPT_VERSION, "\n", sep = "")
  cat("Source run: ", SOURCE_RUN, "\n", sep = "")
  cat("Output: ", OUT_DIR, "\n", sep = "")
  cat("Started: ", format(Sys.time(), tz = "UTC", usetz = TRUE), "\n", sep = "")
  cat("============================================================\n\n")

  source_files <- c(
    status = file.path(SOURCE_RUN, "STATUS.txt"),
    parameters = file.path(SOURCE_RUN, "tables", "08B_parameters.tsv"),
    preflight = file.path(
      SOURCE_RUN, "tables", "08B_preflight_validation.tsv"
    ),
    pair_validation = file.path(
      SOURCE_RUN, "tables", "08B_icamp_pair_structure_validation.tsv"
    ),
    pairwise = file.path(
      SOURCE_RUN, "tables", "08B_icamp_pairwise_processes.tsv.gz"
    ),
    prepared = file.path(
      SOURCE_RUN, "rds", "08B_prepared_icamp_inputs.rds"
    )
  )
  missing <- names(source_files)[!file.exists(source_files)]
  if (length(missing) > 0L) {
    stop("Missing source files: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }

  status <- trimws(readLines(source_files["status"], warn = FALSE))
  parameters <- read_tsv(source_files["parameters"])
  preflight <- read_tsv(source_files["preflight"])
  pair_validation <- read_tsv(source_files["pair_validation"])
  pairwise <- read_tsv(source_files["pairwise"])
  prepared <- readRDS(source_files["prepared"])

  required_pair_columns <- c(
    "sample1", "sample2", "realm1", "realm2", "comparison_scope",
    PROCESS_COLUMNS
  )
  missing_pair_columns <- setdiff(required_pair_columns, names(pairwise))
  if (length(missing_pair_columns) > 0L) {
    stop("Pairwise table lacks: ",
         paste(missing_pair_columns, collapse = ", "), call. = FALSE)
  }
  if (is.null(prepared$metadata) || is.null(prepared$meta_group)) {
    stop("Prepared RDS lacks metadata or meta_group.", call. = FALSE)
  }

  metadata <- as.data.frame(prepared$metadata, stringsAsFactors = FALSE)
  run_column <- find_column(
    metadata,
    c("run_accession", "Run", ".sample_name", "sample_name"),
    "canonical run_accession"
  )
  study_column <- find_column(
    metadata,
    c(
      "study_accession", "Study_accession", "BioProject", "bioproject",
      "bioproject_accession", "project_accession", "Study", "study"
    ),
    "Study/BioProject"
  )
  locality_column <- find_column(
    metadata,
    c(
      "Locality", "locality", "LOCALITY", "sampling_locality",
      "sampling_location", "site", "Site", "location"
    ),
    "locality"
  )

  metadata$run_accession <- clean_chr(metadata[[run_column]])
  if (anyNA(metadata$run_accession) || anyDuplicated(metadata$run_accession)) {
    stop("Prepared metadata has missing or duplicated run_accession.",
         call. = FALSE)
  }

  realm_map <- setNames(
    clean_chr(prepared$meta_group[[1L]]),
    clean_chr(rownames(prepared$meta_group))
  )
  metadata$realm <- unname(realm_map[metadata$run_accession])
  if (anyNA(metadata$realm)) {
    stop("Realm mapping failed for prepared metadata.", call. = FALSE)
  }

  metadata$study_original <- clean_chr(metadata[[study_column]])
  metadata$locality_original <- clean_chr(metadata[[locality_column]])
  metadata$study_missing_original <- is.na(metadata$study_original)
  metadata$locality_missing_original <- is.na(metadata$locality_original)

  metadata$study_block <- metadata$study_original
  metadata$study_block[metadata$study_missing_original] <- paste0(
    "MISSING_STUDY::",
    metadata$run_accession[metadata$study_missing_original]
  )
  metadata$locality_block <- ifelse(
    metadata$locality_missing_original,
    paste0("MISSING_LOCALITY::", metadata$run_accession),
    paste(metadata$realm, metadata$locality_original, sep = " || ")
  )
  metadata$study_locality_block <- paste(
    metadata$study_block, metadata$locality_block, sep = " || "
  )

  sample_ids <- sort(unique(c(pairwise$sample1, pairwise$sample2)))
  missing_metadata <- setdiff(sample_ids, metadata$run_accession)
  extra_metadata <- setdiff(metadata$run_accession, sample_ids)
  if (length(missing_metadata) > 0L || length(extra_metadata) > 0L) {
    stop(
      sprintf(
        "Pairwise/metadata mismatch. Missing metadata=%d; extra metadata=%d.",
        length(missing_metadata), length(extra_metadata)
      ),
      call. = FALSE
    )
  }

  study_map <- setNames(metadata$study_block, metadata$run_accession)
  locality_map <- setNames(metadata$locality_block, metadata$run_accession)
  endpoint_realm_map <- setNames(metadata$realm, metadata$run_accession)
  study_locality_map <- setNames(
    metadata$study_locality_block,
    metadata$run_accession
  )

  pairwise$study1 <- unname(study_map[pairwise$sample1])
  pairwise$study2 <- unname(study_map[pairwise$sample2])
  pairwise$realm_metadata1 <- unname(endpoint_realm_map[pairwise$sample1])
  pairwise$realm_metadata2 <- unname(endpoint_realm_map[pairwise$sample2])
  pairwise$locality1 <- unname(locality_map[pairwise$sample1])
  pairwise$locality2 <- unname(locality_map[pairwise$sample2])
  pairwise$study_locality1 <- unname(study_locality_map[pairwise$sample1])
  pairwise$study_locality2 <- unname(study_locality_map[pairwise$sample2])
  pairwise <- add_fractions(pairwise)

  process_matrix <- as.matrix(pairwise[, PROCESS_COLUMNS, drop = FALSE])
  storage.mode(process_matrix) <- "double"
  pair_keys <- ifelse(
    pairwise$sample1 <= pairwise$sample2,
    paste(pairwise$sample1, pairwise$sample2, sep = "||"),
    paste(pairwise$sample2, pairwise$sample1, sep = "||")
  )

  input_audit <- data.frame(
    check = c(
      "source_STATUS_PASS", "source_rand_999",
      "source_preflight_all_TRUE", "source_pair_validation_all_TRUE",
      "metadata_run_accession_unique", "pairwise_metadata_exact_match",
      "pairwise_realms_match_prepared_metadata",
      "pairwise_rows_79800", "pairwise_samples_400",
      "no_duplicate_unordered_pairs", "within_realm_rows_24692",
      "process_values_finite", "process_values_in_unit_interval",
      "processes_sum_to_one", "selection_nonselection_sum_to_one"
    ),
    pass = c(
      any(status == "PASS"),
      identical(
        suppressWarnings(as.integer(long_value(parameters, "rand", "parameter"))),
        999L
      ),
      all(as_bool(preflight$pass)),
      all(as_bool(pair_validation$pass)),
      !anyNA(metadata$run_accession) && !anyDuplicated(metadata$run_accession),
      length(missing_metadata) == 0L && length(extra_metadata) == 0L,
      identical(as.character(pairwise$realm1), pairwise$realm_metadata1) &&
        identical(as.character(pairwise$realm2), pairwise$realm_metadata2),
      nrow(pairwise) == 79800L,
      length(sample_ids) == 400L,
      !anyDuplicated(pair_keys),
      sum(pairwise$comparison_scope == "within_realm") == 24692L,
      all(is.finite(process_matrix)),
      all(process_matrix >= 0 & process_matrix <= 1),
      max(abs(rowSums(process_matrix) - 1)) <= 1e-6,
      max(abs(
        pairwise$selective_fraction + pairwise$nonselective_fraction - 1
      )) <= 1e-6
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(input_audit, file.path(TABLES_DIR, "08I_input_audit.tsv"))
  if (any(!input_audit$pass)) {
    stop("Input audit failed: ",
         paste(input_audit$check[!input_audit$pass], collapse = ", "),
         call. = FALSE)
  }

  within <- pairwise[
    pairwise$comparison_scope == "within_realm",
    , drop = FALSE
  ]

  realm_sizes <- build_group_sizes(
    metadata, "realm", "realm", MIN_REALM_SAMPLES
  )
  study_sizes <- build_group_sizes(
    metadata, "study", "study_block", MIN_STUDY_SAMPLES
  )
  locality_sizes <- build_group_sizes(
    metadata, "locality", "locality_block", MIN_LOCALITY_SAMPLES
  )
  study_locality_sizes <- build_group_sizes(
    metadata,
    "study_locality",
    "study_locality_block",
    MIN_LOCALITY_SAMPLES
  )
  group_sizes <- rbind(
    realm_sizes, study_sizes, locality_sizes, study_locality_sizes
  )

  realm_means <- build_group_means(
    within, realm_sizes, "realm", "realm1", "realm2", minimum_pairs = 3L
  )
  study_means <- build_group_means(
    within, study_sizes, "study", "study1", "study2", minimum_pairs = 3L
  )
  locality_means <- build_group_means(
    within, locality_sizes, "locality", "locality1", "locality2",
    minimum_pairs = 3L
  )
  study_locality_means <- build_group_means(
    within, study_locality_sizes, "study_locality",
    "study_locality1", "study_locality2", minimum_pairs = 3L
  )
  all_group_means <- rbind(
    realm_means, study_means, locality_means, study_locality_means
  )

  write_tsv(
    metadata[, c(
      "run_accession", "realm", "study_original", "study_block",
      "locality_original", "locality_block", "study_locality_block",
      "study_missing_original", "locality_missing_original"
    )],
    file.path(TABLES_DIR, "08I_sample_metadata_used.tsv")
  )
  write_tsv(group_sizes, file.path(TABLES_DIR, "08I_group_sizes.tsv"))
  write_tsv(
    all_group_means,
    file.path(TABLES_DIR, "08I_group_process_means.tsv")
  )

  group_mean_list <- list(
    realm = realm_means,
    study = study_means,
    locality = locality_means,
    study_locality = study_locality_means
  )

  inference_list <- lapply(seq_along(group_mean_list), function(i) {
    equal_group_summary(
      group_mean_list[[i]],
      n_boot = N_BOOT,
      n_perm = N_PERM,
      seed = SEED + i * 100L
    )
  })

  equal_process <- do.call(rbind, lapply(inference_list, `[[`, "process"))
  nonselective_inference <- do.call(
    rbind,
    lapply(inference_list, `[[`, "nonselective")
  )

  pair_weighted <- data.frame(
    metric = c(PROCESS_COLUMNS, "selective_fraction", "nonselective_fraction"),
    estimate = colMeans(
      within[, c(PROCESS_COLUMNS, "selective_fraction",
                 "nonselective_fraction"), drop = FALSE]
    ),
    bootstrap_lwr = NA_real_,
    bootstrap_upr = NA_real_,
    n_groups = NA_integer_,
    n_boot = NA_integer_,
    unit = "all_within_realm_pairs",
    weighting = "pair_weighted_descriptive_only",
    stringsAsFactors = FALSE
  )
  weighting_process_summary <- rbind(pair_weighted, equal_process)
  write_tsv(
    weighting_process_summary,
    file.path(TABLES_DIR, "08I_weighting_process_summary.tsv")
  )
  write_tsv(
    nonselective_inference,
    file.path(TABLES_DIR, "08I_nonselective_group_inference.tsv")
  )

  equal_loo <- do.call(
    rbind,
    lapply(group_mean_list, equal_group_loo)
  )
  endpoint_loo <- rbind(
    endpoint_removal_loo(within, realm_means, "realm1", "realm2"),
    endpoint_removal_loo(within, study_means, "study1", "study2"),
    endpoint_removal_loo(within, locality_means, "locality1", "locality2"),
    endpoint_removal_loo(
      within, study_locality_means,
      "study_locality1", "study_locality2"
    )
  )
  all_loo <- rbind(equal_loo, endpoint_loo)
  write_tsv(all_loo, file.path(TABLES_DIR, "08I_leave_one_group_out.tsv"))

  loo_summary <- do.call(
    rbind,
    lapply(
      split(all_loo, interaction(all_loo$unit, all_loo$loo_type, drop = TRUE)),
      function(x) data.frame(
        unit = x$unit[1L],
        loo_type = x$loo_type[1L],
        n_iterations = nrow(x),
        minimum_nonselective_fraction = min(x$nonselective_fraction),
        median_nonselective_fraction = stats::median(x$nonselective_fraction),
        maximum_nonselective_fraction = max(x$nonselective_fraction),
        all_iterations_above_0.5 = all(x$nonselective_fraction > 0.5),
        stringsAsFactors = FALSE
      )
    )
  )
  rownames(loo_summary) <- NULL
  write_tsv(
    loo_summary,
    file.path(TABLES_DIR, "08I_leave_one_group_out_summary.tsv")
  )

  balanced <- balanced_realm_resampling(
    within,
    metadata,
    realm_column = "realm",
    minimum_realm_samples = MIN_REALM_SAMPLES,
    n_replicates = N_BALANCE,
    seed = SEED + 5000L
  )
  eligible_balanced_realms <- attr(balanced, "eligible_realms")
  write_tsv_gz(
    balanced,
    file.path(TABLES_DIR, "08I_balanced_realm_subsampling.tsv.gz")
  )
  balanced_summary <- data.frame(
    n_replicates = nrow(balanced),
    eligible_realms = paste(eligible_balanced_realms, collapse = ";"),
    n_realms = balanced$n_realms[1L],
    samples_per_realm = balanced$samples_per_realm[1L],
    pairs_per_realm = balanced$pairs_per_realm[1L],
    mean_nonselective_fraction = mean(
      balanced$equal_realm_nonselective_fraction
    ),
    q025_nonselective_fraction = as.numeric(stats::quantile(
      balanced$equal_realm_nonselective_fraction, 0.025, names = FALSE
    )),
    median_nonselective_fraction = stats::median(
      balanced$equal_realm_nonselective_fraction
    ),
    q975_nonselective_fraction = as.numeric(stats::quantile(
      balanced$equal_realm_nonselective_fraction, 0.975, names = FALSE
    )),
    minimum_nonselective_fraction = min(
      balanced$equal_realm_nonselective_fraction
    ),
    proportion_replicates_above_0.5 = mean(
      balanced$equal_realm_nonselective_fraction > 0.5
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(
    balanced_summary,
    file.path(TABLES_DIR, "08I_balanced_realm_subsampling_summary.tsv")
  )

  primary_units <- c("realm", "study", "locality")
  primary_boot <- equal_process[
    equal_process$unit %in% primary_units &
      equal_process$metric == "nonselective_fraction",
    , drop = FALSE
  ]
  primary_loo <- loo_summary[loo_summary$unit %in% primary_units, , drop = FALSE]

  decision <- data.frame(
    criterion = c(
      "pair_weighted_nonselective_above_0.5",
      "equal_realm_bootstrap_lower_above_0.5",
      "equal_study_bootstrap_lower_above_0.5",
      "equal_locality_bootstrap_lower_above_0.5",
      "all_primary_LOO_iterations_above_0.5",
      "balanced_realm_q025_above_0.5",
      "balanced_realm_all_replicates_above_0.5"
    ),
    pass = c(
      pair_weighted$estimate[
        pair_weighted$metric == "nonselective_fraction"
      ] > 0.5,
      primary_boot$bootstrap_lwr[primary_boot$unit == "realm"] > 0.5,
      primary_boot$bootstrap_lwr[primary_boot$unit == "study"] > 0.5,
      primary_boot$bootstrap_lwr[primary_boot$unit == "locality"] > 0.5,
      all(primary_loo$all_iterations_above_0.5),
      balanced_summary$q025_nonselective_fraction > 0.5,
      identical(balanced_summary$proportion_replicates_above_0.5, 1)
    ),
    stringsAsFactors = FALSE
  )
  decision$observed <- c(
    pair_weighted$estimate[
      pair_weighted$metric == "nonselective_fraction"
    ],
    primary_boot$bootstrap_lwr[primary_boot$unit == "realm"],
    primary_boot$bootstrap_lwr[primary_boot$unit == "study"],
    primary_boot$bootstrap_lwr[primary_boot$unit == "locality"],
    min(primary_loo$minimum_nonselective_fraction),
    balanced_summary$q025_nonselective_fraction,
    balanced_summary$proportion_replicates_above_0.5
  )
  decision$expected <- c(
    ">0.5", ">0.5", ">0.5", ">0.5", ">0.5", ">0.5", "1"
  )
  write_tsv(decision, file.path(TABLES_DIR, "08I_robustness_decision.tsv"))

  overall_pass <- all(decision$pass)
  summary_table <- data.frame(
    item = c(
      "status", "script", "script_version", "source_run",
      "source_run_sha256_pairwise", "n_samples", "n_within_realm_pairs",
      "n_realms_total", "n_realms_primary", "n_studies_total",
      "n_studies_eligible", "n_localities_total",
      "n_localities_eligible", "n_missing_study_original",
      "n_missing_locality_original", "primary_claim",
      "overall_robustness_pass", "output_directory"
    ),
    value = c(
      ifelse(overall_pass, "PASS", "REVIEW"),
      SCRIPT_ID, SCRIPT_VERSION, SOURCE_RUN,
      sha256_file(source_files["pairwise"]),
      nrow(metadata), nrow(within), nrow(realm_sizes),
      sum(realm_means$eligible), nrow(study_sizes),
      sum(study_means$eligible), nrow(locality_sizes),
      sum(locality_means$eligible),
      sum(metadata$study_missing_original),
      sum(metadata$locality_missing_original),
      "nonselective_fraction_above_0.5",
      overall_pass, OUT_DIR
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(summary_table, file.path(TABLES_DIR, "08I_summary.tsv"))

  parameters_out <- data.frame(
    parameter = c(
      "script", "script_version", "source_run", "run_accession_column",
      "study_column", "locality_column", "n_boot", "n_perm",
      "n_balance", "min_realm_samples", "min_study_samples",
      "min_locality_samples", "seed", "selective_definition",
      "nonselective_definition"
    ),
    value = c(
      SCRIPT_ID, SCRIPT_VERSION, SOURCE_RUN, run_column,
      study_column, locality_column, N_BOOT, N_PERM, N_BALANCE,
      MIN_REALM_SAMPLES, MIN_STUDY_SAMPLES, MIN_LOCALITY_SAMPLES, SEED,
      paste(SELECTIVE_COLUMNS, collapse = "+"),
      paste(NONSELECTIVE_COLUMNS, collapse = "+")
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(parameters_out, file.path(TABLES_DIR, "08I_parameters.tsv"))
  capture.output(sessionInfo(), file = file.path(LOGS_DIR, "08I_sessionInfo.txt"))

  cat("\n=== GROUP INFERENCE FOR THE FRACTION NOT ASSIGNED TO SELECTION ===\n")
  print(nonselective_inference, row.names = FALSE)
  cat("\n=== LEAVE-ONE-GROUP-OUT SUMMARY ===\n")
  print(loo_summary, row.names = FALSE)
  cat("\n=== BALANCED REALM SUBSAMPLING ===\n")
  print(balanced_summary, row.names = FALSE)
  cat("\n=== ROBUSTNESS DECISION ===\n")
  print(decision, row.names = FALSE)

  write_status(OUT_DIR, ifelse(overall_pass, "PASS", "REVIEW"))
  cat("\nFINAL STATUS: ", ifelse(overall_pass, "PASS", "REVIEW"), "\n", sep = "")
  invisible(overall_pass)
}

run_error <- NULL
tryCatch(
  {
    main()
    status_ok <- TRUE
  },
  error = function(e) {
    run_error <<- conditionMessage(e)
    message("ERROR: ", run_error)
  }
)

sink(type = "message")
sink()
close(log_connection)

if (!status_ok) {
  write_status(OUT_DIR, "FAILED", run_error)
  writeLines(run_error, file.path(OUT_DIR, "ERROR.txt"))
  stop(run_error, call. = FALSE)
}
