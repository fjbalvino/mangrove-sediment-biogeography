#!/usr/bin/env Rscript

# =============================================================================
# 156_PREPARE_KO_CLR_MATRICES_FOR_ASSOCIATION_AND_NETWORK.R
#
# Purpose
# -------
# Build the validated KO count, relative-abundance and CLR matrices used by the
# functional association models and the final KO consensus network.
#
# Reproducibility corrections relative to the historical script
# ---------------------------------------------------------------
# 1. Prevalent/Shell/Cloud classes are computed directly from the aligned raw
#    KO count matrix. The workflow no longer depends on the obsolete historical
#    script "154_build_KO_functional_prevalence_classes" (the 154 identifier is
#    now used by the reconciled functional PCA analysis).
# 2. run_accession is the canonical sample key. Relaxed ID matching is disabled
#    by default and, when explicitly enabled, every substitution is audited.
# 3. Missing-like matrix cells are counted and reported. Their conversion to
#    zero must be explicitly allowed (default TRUE for sparse exported tables).
# 4. Canonical and legacy output names are written together so that scripts
#    171, 178 and 179 remain compatible.
#
# Canonical thresholds
# --------------------
# Presence:     count > 0
# Prevalent:    prevalence >= 0.95
# Shell:        0.15 <= prevalence < 0.95
# Cloud:        prevalence < 0.15
# Association: prevalence_n >= 10 and total_count >= 10
# Network:     prevalence_n >= 20 and total_count >= 50
# Network cap: 5,000 KOs, ranked deterministically by CLR variance
# CLR:         pseudocount = 1
#
# Required input
# --------------
# --ko_matrix PATH    Feature-by-sample KO count table (KO in the first or a
#                     named KO column; sample columns identified by run_accession)
#
# Metadata is supplied with --metadata PATH or discovered in the latest script
# 134 run. Use --help for all options.
# =============================================================================

options(warn = 1)

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

stopf <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)
warnf <- function(fmt, ...) warning(sprintf(fmt, ...), call. = FALSE)

clean_chr <- function(x) {
  y <- trimws(as.character(x))
  y[y %in% c("", "NA", "NaN", "NULL", "None", ".")] <- NA_character_
  y
}

is_missing_like <- function(x) {
  is.na(x) | trimws(as.character(x)) %in%
    c("", "NA", "NaN", "NULL", "None", ".")
}

read_latest <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  x <- clean_chr(readLines(path, warn = FALSE))
  x <- x[!is.na(x)]
  if (length(x) == 0L) NA_character_ else x[1L]
}

read_tsv_auto <- function(path) {
  if (!file.exists(path)) stopf("Input file does not exist: %s", path)
  if (requireNamespace("data.table", quietly = TRUE)) {
    return(as.data.frame(
      data.table::fread(path, sep = "\t", header = TRUE, data.table = FALSE,
                        check.names = FALSE, na.strings = c("", "NA", "NaN")),
      check.names = FALSE
    ))
  }
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else path
  if (inherits(con, "connection")) on.exit(close(con), add = TRUE)
  utils::read.delim(con, sep = "\t", header = TRUE, check.names = FALSE,
                    stringsAsFactors = FALSE, quote = "", comment.char = "",
                    na.strings = c("", "NA", "NaN"))
}

write_tsv <- function(x, path) {
  utils::write.table(x, file = path, sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = TRUE, na = "")
}

write_matrix_tsv_gz <- function(mat, path, row_id = "KO") {
  df <- data.frame(feature_id = rownames(mat), mat, check.names = FALSE)
  names(df)[1L] <- row_id
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  utils::write.table(df, file = con, sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = TRUE, na = "")
}

require_columns <- function(x, cols, label) {
  absent <- setdiff(cols, names(x))
  if (length(absent) > 0L) {
    stopf("%s is missing required column(s): %s",
          label, paste(absent, collapse = ", "))
  }
}

find_column <- function(x, candidates, label, required = TRUE) {
  exact <- candidates[candidates %in% names(x)]
  if (length(exact) > 0L) return(exact[1L])
  lower_names <- tolower(names(x))
  idx <- match(tolower(candidates), lower_names)
  idx <- idx[!is.na(idx)]
  if (length(idx) > 0L) return(names(x)[idx[1L]])
  if (required) {
    stopf("Could not identify %s. Tried: %s",
          label, paste(candidates, collapse = ", "))
  }
  NA_character_
}

normalise_id <- function(x) {
  tolower(gsub("[^a-zA-Z0-9]", "", trimws(as.character(x))))
}

match_sample_ids <- function(target_ids, source_ids, allow_relaxed = FALSE) {
  target_ids <- clean_chr(target_ids)
  source_ids <- clean_chr(source_ids)
  if (anyNA(target_ids) || anyDuplicated(target_ids)) {
    stop("Target run_accession values must be non-missing and unique.", call. = FALSE)
  }
  if (anyNA(source_ids) || anyDuplicated(source_ids)) {
    stop("KO-matrix sample-column names must be non-missing and unique.", call. = FALSE)
  }

  index <- match(target_ids, source_ids)
  method <- ifelse(is.na(index), "unmatched", "exact")

  if (allow_relaxed && anyNA(index)) {
    norm_source <- normalise_id(source_ids)
    for (i in which(is.na(index))) {
      candidates <- which(norm_source == normalise_id(target_ids[i]))
      if (length(candidates) > 1L) {
        stopf("Relaxed ID match for '%s' is ambiguous: %s",
              target_ids[i], paste(source_ids[candidates], collapse = ", "))
      }
      if (length(candidates) == 1L) {
        index[i] <- candidates
        method[i] <- "relaxed_normalised"
      }
    }
  }

  audit <- data.frame(
    run_accession = target_ids,
    matched_matrix_column = ifelse(is.na(index), NA_character_, source_ids[index]),
    match_method = method,
    stringsAsFactors = FALSE
  )
  list(index = index, audit = audit)
}

numeric_matrix_from_df <- function(x, sample_names, missing_as_zero) {
  out <- matrix(NA_real_, nrow = nrow(x), ncol = ncol(x),
                dimnames = list(NULL, sample_names))
  audit <- vector("list", ncol(x))
  bad_rows <- list()

  for (j in seq_len(ncol(x))) {
    raw <- x[[j]]
    missing <- is_missing_like(raw)
    numeric <- suppressWarnings(as.numeric(as.character(raw)))
    bad <- !missing & is.na(numeric)
    if (any(bad)) {
      bad_rows[[length(bad_rows) + 1L]] <- data.frame(
        matrix_column = sample_names[j],
        row_number = which(bad),
        raw_value = as.character(raw[bad]),
        stringsAsFactors = FALSE
      )
    }
    if (any(missing) && !missing_as_zero) {
      stopf("Matrix column '%s' contains %d missing-like cell(s). Re-run with --missing_as_zero 1 only if zeros are intended.",
            sample_names[j], sum(missing))
    }
    numeric[missing] <- 0
    out[, j] <- numeric
    audit[[j]] <- data.frame(
      matrix_column = sample_names[j],
      n_cells = length(raw),
      n_missing_converted_to_zero = sum(missing),
      n_non_numeric = sum(bad),
      stringsAsFactors = FALSE
    )
  }

  bad_df <- if (length(bad_rows) == 0L) {
    data.frame(matrix_column = character(), row_number = integer(),
               raw_value = character(), stringsAsFactors = FALSE)
  } else {
    do.call(rbind, bad_rows)
  }
  list(matrix = out, audit = do.call(rbind, audit), bad = bad_df)
}

collapse_duplicate_KOs <- function(mat, ko) {
  ko <- toupper(clean_chr(ko))
  if (anyNA(ko)) stop("KO identifiers contain missing values.", call. = FALSE)
  invalid <- !grepl("^K[0-9]{5}$", ko)
  if (any(invalid)) {
    stopf("Invalid KO identifiers detected (first values: %s)",
          paste(utils::head(unique(ko[invalid]), 10L), collapse = ", "))
  }
  if (!anyDuplicated(ko)) {
    rownames(mat) <- ko
    return(list(matrix = mat, duplicated_rows_collapsed = 0L))
  }
  collapsed <- rowsum(mat, group = ko, reorder = TRUE, na.rm = FALSE)
  list(matrix = collapsed,
       duplicated_rows_collapsed = length(ko) - length(unique(ko)))
}

read_KO_matrix_features_rows <- function(path, missing_as_zero) {
  raw <- read_tsv_auto(path)
  if (nrow(raw) == 0L || ncol(raw) < 2L) {
    stop("The KO matrix must contain at least one KO row and one sample column.",
         call. = FALSE)
  }
  ko_col <- find_column(
    raw,
    c("KO", "ko", "KEGG_Ortholog", "KEGG_ortholog", "feature_id", "#OTU ID"),
    "KO identifier column",
    required = FALSE
  )
  if (is.na(ko_col)) ko_col <- names(raw)[1L]
  sample_cols <- setdiff(names(raw), ko_col)
  converted <- numeric_matrix_from_df(raw[sample_cols], sample_cols, missing_as_zero)
  collapsed <- collapse_duplicate_KOs(converted$matrix, raw[[ko_col]])
  list(
    matrix = collapsed$matrix,
    conversion_audit = converted$audit,
    non_numeric = converted$bad,
    ko_column = ko_col,
    duplicated_rows_collapsed = collapsed$duplicated_rows_collapsed
  )
}

calc_relative_abundance <- function(counts) {
  totals <- colSums(counts)
  if (any(!is.finite(totals) | totals <= 0)) {
    stopf("%d sample(s) have non-positive or non-finite KO library sizes.",
          sum(!is.finite(totals) | totals <= 0))
  }
  sweep(counts, 2L, totals, "/")
}

calc_clr <- function(counts, pseudocount) {
  if (!is.finite(pseudocount) || pseudocount <= 0) {
    stop("pseudocount must be finite and > 0.", call. = FALSE)
  }
  logged <- log(counts + pseudocount)
  sweep(logged, 2L, colMeans(logged), "-")
}

assign_prevalence_class <- function(prevalence_fraction,
                                    prevalent_threshold,
                                    shell_threshold) {
  ifelse(
    prevalence_fraction >= prevalent_threshold,
    "Prevalent",
    ifelse(prevalence_fraction >= shell_threshold, "Shell", "Cloud")
  )
}

make_filter_stats <- function(counts, presence_threshold,
                              prevalent_threshold, shell_threshold) {
  prevalence_n <- rowSums(counts > presence_threshold)
  prevalence_fraction <- prevalence_n / ncol(counts)
  total_count <- rowSums(counts)
  data.frame(
    KO = rownames(counts),
    prevalence_n = as.integer(prevalence_n),
    prevalence_fraction = as.numeric(prevalence_fraction),
    total_count = as.numeric(total_count),
    functional_prevalence_class = assign_prevalence_class(
      prevalence_fraction, prevalent_threshold, shell_threshold
    ),
    stringsAsFactors = FALSE
  )
}

count_by_class <- function(x, filter_set) {
  levels <- c("Prevalent", "Shell", "Cloud")
  counts <- table(factor(x$functional_prevalence_class, levels = levels))
  data.frame(filter_set = filter_set,
             functional_prevalence_class = levels,
             n_KOs = as.integer(counts), stringsAsFactors = FALSE)
}

standardise_reference_annotation <- function(path) {
  x <- read_tsv_auto(path)
  ko_col <- find_column(x, c("KO", "ko", "KEGG_Ortholog"), "reference KO column")
  class_col <- find_column(
    x,
    c("functional_prevalence_class", "prevalence_class", "KO_class",
      "core_class", "class"),
    "reference prevalence-class column"
  )
  out <- data.frame(
    KO = toupper(clean_chr(x[[ko_col]])),
    reference_class = clean_chr(x[[class_col]]),
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$KO), , drop = FALSE]
  if (anyDuplicated(out$KO)) {
    duplicated_KOs <- unique(out$KO[duplicated(out$KO) | duplicated(out$KO, fromLast = TRUE)])
    conflict <- vapply(duplicated_KOs, function(k) {
      length(unique(stats::na.omit(out$reference_class[out$KO == k]))) > 1L
    }, logical(1))
    if (any(conflict)) {
      stopf("Reference annotation assigns conflicting classes to: %s",
            paste(utils::head(duplicated_KOs[conflict], 10L), collapse = ", "))
    }
    out <- out[!duplicated(out$KO), , drop = FALSE]
  }
  out
}

validate_parameters <- function(cfg) {
  if (!is.finite(cfg$presence_threshold) || cfg$presence_threshold < 0) {
    stop("presence_threshold must be finite and >= 0.", call. = FALSE)
  }
  if (!is.finite(cfg$shell_threshold) || !is.finite(cfg$prevalent_threshold) ||
      cfg$shell_threshold < 0 || cfg$prevalent_threshold > 1 ||
      cfg$shell_threshold >= cfg$prevalent_threshold) {
    stop("Require 0 <= shell_threshold < prevalent_threshold <= 1.", call. = FALSE)
  }
  integer_fields <- c("min_prev_assoc", "min_prev_network", "max_network_kos")
  for (field in integer_fields) {
    if (is.na(cfg[[field]]) || cfg[[field]] < 1L) {
      stopf("%s must be an integer >= 1.", field)
    }
  }
  if (!is.finite(cfg$min_total_assoc) || cfg$min_total_assoc < 0 ||
      !is.finite(cfg$min_total_network) || cfg$min_total_network < 0) {
    stop("Total-count thresholds must be finite and >= 0.", call. = FALSE)
  }
}

print_help <- function() {
  cat(paste(
    "Usage:",
    "  Rscript 156_PREPARE_KO_CLR_MATRICES_FOR_ASSOCIATION_AND_NETWORK.R \\",
    "    --ko_matrix PATH [--metadata PATH] [options]",
    "",
    "Key options:",
    "  --base PATH                    Project base directory",
    "  --results_dir PATH             Results directory (default: BASE/resultados)",
    "  --run134 PATH                  Script-134 output directory",
    "  --metadata PATH                Curated functional metadata TSV",
    "  --metadata_id_col NAME         Metadata sample-key column (auto-detected)",
    "  --prevalence_annotation PATH   Optional independent class reference",
    "  --allow_relaxed_id_match 0|1   Default 0",
    "  --missing_as_zero 0|1          Default 1; conversion is audited",
    "  --strict_reference 0|1         Fail on reference-class discrepancies",
    "  --strict_expected 0|1          Enforce manuscript dimensions (default 1)",
    sep = "\n"
  ))
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if ("--help" %in% args || "-h" %in% args) {
    print_help()
    return(invisible(NULL))
  }

  base <- get_arg(
    "--base",
    Sys.getenv("BIOGEOGRAPHY_BASE", "/home/fjbalvino/Biogeografia_2026")
  )
  results_dir <- get_arg("--results_dir", file.path(base, "resultados"))
  run134 <- get_arg(
    "--run134",
    read_latest(file.path(results_dir,
                          "LATEST_134_curate_metadata_blocks_and_axes.txt"))
  )
  metadata_path <- get_arg("--metadata", NA_character_)
  if (is.na(metadata_path) || metadata_path == "") {
    if (!is.na(run134) && dir.exists(run134)) {
      metadata_path <- file.path(
        run134, "tables", "134_metadata_functional_curated_axes.tsv"
      )
    }
  }

  cfg <- list(
    presence_threshold = as.numeric(get_arg("--presence_threshold", "0")),
    prevalent_threshold = as.numeric(get_arg("--prevalent_threshold", "0.95")),
    shell_threshold = as.numeric(get_arg("--shell_threshold", "0.15")),
    min_prev_assoc = as.integer(get_arg("--min_prev_assoc", "10")),
    min_total_assoc = as.numeric(get_arg("--min_total_assoc", "10")),
    min_prev_network = as.integer(get_arg("--min_prev_network", "20")),
    min_total_network = as.numeric(get_arg("--min_total_network", "50")),
    max_network_kos = as.integer(get_arg("--max_network_kos", "5000")),
    pseudocount = as.numeric(get_arg("--pseudocount", "1")),
    allow_relaxed_id_match = as_flag(get_arg("--allow_relaxed_id_match", "0")),
    missing_as_zero = as_flag(get_arg("--missing_as_zero", "1")),
    strict_reference = as_flag(get_arg("--strict_reference", "1")),
    strict_expected = as_flag(get_arg("--strict_expected", "1")),
    expected_n_samples = as.integer(get_arg("--expected_n_samples", "194")),
    expected_n_kos = as.integer(get_arg("--expected_n_kos", "17508")),
    expected_n_assoc = as.integer(get_arg("--expected_n_assoc", "13351")),
    expected_n_network = as.integer(get_arg("--expected_n_network", "5000")),
    expected_all_prevalent = as.integer(get_arg("--expected_all_prevalent", "1106")),
    expected_all_shell = as.integer(get_arg("--expected_all_shell", "9853")),
    expected_all_cloud = as.integer(get_arg("--expected_all_cloud", "6549")),
    expected_assoc_prevalent = as.integer(get_arg("--expected_assoc_prevalent", "1106")),
    expected_assoc_shell = as.integer(get_arg("--expected_assoc_shell", "9853")),
    expected_assoc_cloud = as.integer(get_arg("--expected_assoc_cloud", "2392"))
  )
  validate_parameters(cfg)

  ko_matrix_path <- get_arg("--ko_matrix", NA_character_)
  reference_path <- get_arg("--prevalence_annotation", NA_character_)
  metadata_id_col_arg <- get_arg("--metadata_id_col", NA_character_)

  if (is.na(ko_matrix_path) || ko_matrix_path == "") {
    stop("--ko_matrix is required; the raw KO matrix is no longer inferred from an obsolete script-154 run.",
         call. = FALSE)
  }
  if (is.na(metadata_path) || metadata_path == "" || !file.exists(metadata_path)) {
    stop("Curated metadata was not found. Supply --metadata PATH or a valid --run134 directory.",
         call. = FALSE)
  }
  ko_matrix_path <- normalizePath(ko_matrix_path, winslash = "/", mustWork = TRUE)
  metadata_path <- normalizePath(metadata_path, winslash = "/", mustWork = TRUE)
  if (!is.na(reference_path) && reference_path != "") {
    reference_path <- normalizePath(reference_path, winslash = "/", mustWork = TRUE)
  }

  script_name <- "156_PREPARE_KO_CLR_MATRICES_FOR_ASSOCIATION_AND_NETWORK"
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  outdir <- file.path(results_dir, paste0(script_name, "_", stamp))
  tables_dir <- file.path(outdir, "tables")
  logs_dir <- file.path(outdir, "logs")
  rds_dir <- file.path(outdir, "rds")
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
  latest_files <- c(
    file.path(results_dir, paste0("LATEST_", script_name, ".txt")),
    file.path(results_dir,
              "LATEST_156_build_KO_CLR_matrix_prevalence_classes.txt")
  )
  invisible(lapply(latest_files, function(path) writeLines(outdir, path)))
  writeLines(c("status\tRUNNING", paste0("started_utc\t", format(
    Sys.time(), tz = "UTC", usetz = TRUE
  ))), file.path(logs_dir, "156_status.tsv"))
  writeLines("RUNNING", file.path(outdir, "STATUS.txt"))

  execute <- function() {
    cat("[1/8] Reading and validating curated metadata...\n")
    metadata <- read_tsv_auto(metadata_path)
    id_col <- if (!is.na(metadata_id_col_arg) && metadata_id_col_arg != "") {
      if (!metadata_id_col_arg %in% names(metadata)) {
        stopf("Requested metadata ID column is absent: %s", metadata_id_col_arg)
      }
      metadata_id_col_arg
    } else {
      find_column(metadata,
                  c("run_accession", ".sample_name", "sample_name", "Run"),
                  "metadata run_accession column")
    }
    metadata$run_accession <- clean_chr(metadata[[id_col]])
    missing_metadata_ids <- which(is.na(metadata$run_accession))
    if (length(missing_metadata_ids) > 0L) {
      write_tsv(data.frame(row_number = missing_metadata_ids),
                file.path(tables_dir, "156_metadata_missing_run_accession.tsv"))
      stopf("Metadata contains %d missing run_accession value(s).",
            length(missing_metadata_ids))
    }
    if (anyDuplicated(metadata$run_accession)) {
      dup <- metadata$run_accession[
        duplicated(metadata$run_accession) |
          duplicated(metadata$run_accession, fromLast = TRUE)
      ]
      write_tsv(data.frame(run_accession = unique(dup)),
                file.path(tables_dir, "156_metadata_duplicated_run_accession.tsv"))
      stop("Metadata run_accession values are not unique.", call. = FALSE)
    }

    cat("[2/8] Reading the KO count matrix and auditing numeric cells...\n")
    ko_res <- read_KO_matrix_features_rows(ko_matrix_path, cfg$missing_as_zero)
    write_tsv(ko_res$conversion_audit,
              file.path(tables_dir, "156_KO_matrix_cell_conversion_audit.tsv"))
    write_tsv(ko_res$non_numeric,
              file.path(tables_dir, "156_non_numeric_KO_matrix_values.tsv"))
    if (nrow(ko_res$non_numeric) > 0L) {
      stopf("The KO matrix contains %d non-numeric, non-missing cell(s).",
            nrow(ko_res$non_numeric))
    }
    counts_raw <- ko_res$matrix
    if (any(!is.finite(counts_raw)) || any(counts_raw < 0)) {
      stop("KO counts must be finite and non-negative.", call. = FALSE)
    }

    cat("[3/8] Aligning samples by canonical run_accession...\n")
    matching <- match_sample_ids(
      metadata$run_accession,
      colnames(counts_raw),
      allow_relaxed = cfg$allow_relaxed_id_match
    )
    write_tsv(matching$audit,
              file.path(tables_dir, "156_KO_matrix_sample_matching_all.tsv"))
    missing_runs <- matching$audit$run_accession[is.na(matching$index)]
    if (length(missing_runs) > 0L) {
      write_tsv(data.frame(run_accession = missing_runs),
                file.path(tables_dir,
                          "156_metadata_runs_missing_from_KO_matrix.tsv"))
      stopf("%d metadata sample(s) are absent from the KO matrix.",
            length(missing_runs))
    }
    used_source <- matching$audit$matched_matrix_column
    unmatched_matrix <- setdiff(colnames(counts_raw), used_source)
    write_tsv(data.frame(matrix_column = unmatched_matrix),
              file.path(tables_dir, "156_KO_matrix_unmatched_columns.tsv"))

    counts_aligned <- counts_raw[, matching$index, drop = FALSE]
    colnames(counts_aligned) <- metadata$run_accession
    metadata <- metadata[match(colnames(counts_aligned), metadata$run_accession),,
                         drop = FALSE]
    if (!identical(metadata$run_accession, colnames(counts_aligned))) {
      stop("Internal error: metadata and KO matrix are not identically ordered.",
           call. = FALSE)
    }
    write_tsv(metadata, file.path(tables_dir, "156_metadata_aligned.tsv"))
    library_stats <- data.frame(
      run_accession = colnames(counts_aligned),
      KO_library_size = colSums(counts_aligned),
      n_detected_KOs = colSums(counts_aligned > cfg$presence_threshold),
      stringsAsFactors = FALSE
    )
    write_tsv(library_stats,
              file.path(tables_dir, "156_functional_library_size_by_sample.tsv"))

    cat("[4/8] Computing KO prevalence classes directly from the aligned matrix...\n")
    nonzero <- rowSums(counts_aligned > cfg$presence_threshold) > 0L
    counts_all <- counts_aligned[nonzero, , drop = FALSE]
    if (nrow(counts_all) == 0L) stop("No detected KOs remain.", call. = FALSE)
    filter_stats <- make_filter_stats(
      counts_all,
      cfg$presence_threshold,
      cfg$prevalent_threshold,
      cfg$shell_threshold
    )
    write_tsv(filter_stats,
              file.path(tables_dir, "156_KO_filter_stats_all_nonzero.tsv"))
    prevalence_annotation <- filter_stats[c(
      "KO", "prevalence_n", "prevalence_fraction", "total_count",
      "functional_prevalence_class"
    )]
    write_tsv(prevalence_annotation,
              file.path(tables_dir, "156_KO_prevalence_class_annotation.tsv"))

    reference_validation <- data.frame(
      KO = character(), computed_class = character(), reference_class = character(),
      match = logical(), stringsAsFactors = FALSE
    )
    if (!is.na(reference_path) && reference_path != "") {
      reference <- standardise_reference_annotation(reference_path)
      reference_validation <- merge(
        prevalence_annotation[c("KO", "functional_prevalence_class")],
        reference,
        by = "KO", all = TRUE, sort = TRUE
      )
      names(reference_validation)[2L] <- "computed_class"
      reference_validation$match <-
        !is.na(reference_validation$computed_class) &
        !is.na(reference_validation$reference_class) &
        reference_validation$computed_class == reference_validation$reference_class
      write_tsv(reference_validation,
                file.path(tables_dir,
                          "156_prevalence_annotation_reference_validation.tsv"))
      discrepancy <- is.na(reference_validation$computed_class) |
        is.na(reference_validation$reference_class) |
        !reference_validation$match
      if (any(discrepancy) && cfg$strict_reference) {
        stopf("Reference prevalence annotation has %d discrepancy/discrepancies.",
              sum(discrepancy))
      }
    } else {
      write_tsv(reference_validation,
                file.path(tables_dir,
                          "156_prevalence_annotation_reference_validation.tsv"))
    }

    cat("[5/8] Building association count, relative-abundance and CLR matrices...\n")
    keep_assoc <- filter_stats$prevalence_n >= cfg$min_prev_assoc &
      filter_stats$total_count >= cfg$min_total_assoc
    counts_assoc <- counts_all[filter_stats$KO[keep_assoc], , drop = FALSE]
    assoc_stats <- filter_stats[keep_assoc, , drop = FALSE]
    if (nrow(counts_assoc) == 0L) stop("No KOs pass the association filter.", call. = FALSE)
    relative_assoc <- calc_relative_abundance(counts_assoc)
    clr_assoc <- calc_clr(counts_assoc, cfg$pseudocount)

    write_matrix_tsv_gz(counts_assoc,
                        file.path(tables_dir, "156_KO_counts_association_filtered.tsv.gz"))
    write_matrix_tsv_gz(relative_assoc,
                        file.path(tables_dir, "156_KO_relative_association_filtered.tsv.gz"))
    write_matrix_tsv_gz(
      relative_assoc,
      file.path(
        tables_dir,
        "156_KO_relative_abundance_association_filtered.tsv.gz"
      )
    )
    write_matrix_tsv_gz(clr_assoc,
                        file.path(tables_dir, "156_KO_CLR_association_filtered.tsv.gz"))
    write_tsv(assoc_stats,
              file.path(tables_dir, "156_KO_association_filter_stats.tsv"))

    cat("[6/8] Selecting the deterministic network KO universe...\n")
    keep_network <- assoc_stats$prevalence_n >= cfg$min_prev_network &
      assoc_stats$total_count >= cfg$min_total_network
    candidate_KOs <- assoc_stats$KO[keep_network]
    candidate_clr <- clr_assoc[candidate_KOs, , drop = FALSE]
    candidate_variance <- apply(candidate_clr, 1L, stats::var)
    rank_order <- order(-candidate_variance, rownames(candidate_clr))
    if (length(rank_order) > cfg$max_network_kos) {
      rank_order <- rank_order[seq_len(cfg$max_network_kos)]
    }
    selected_KOs <- rownames(candidate_clr)[rank_order]
    clr_network <- candidate_clr[selected_KOs, , drop = FALSE]
    counts_network <- counts_assoc[selected_KOs, , drop = FALSE]
    network_stats <- assoc_stats[
      match(selected_KOs, assoc_stats$KO), , drop = FALSE
    ]
    network_stats$CLR_variance <- as.numeric(candidate_variance[selected_KOs])
    network_stats$network_rank_by_CLR_variance <- seq_len(nrow(network_stats))
    network_stats$variance_rank <- network_stats$network_rank_by_CLR_variance
    network_stats$selected_for_network <- TRUE

    # Legacy filenames below are consumed by scripts 171, 178 and 179.
    write_matrix_tsv_gz(clr_network,
                        file.path(tables_dir, "156_KO_CLR_network_filtered.tsv.gz"))
    write_matrix_tsv_gz(counts_network,
                        file.path(tables_dir, "156_KO_counts_network_filtered.tsv.gz"))
    write_tsv(network_stats,
              file.path(tables_dir, "156_KO_network_filter_stats.tsv"))

    cat("[7/8] Validating dimensions, transformations and class counts...\n")
    class_counts <- rbind(
      count_by_class(filter_stats, "all_nonzero"),
      count_by_class(assoc_stats, "association_filtered"),
      count_by_class(network_stats, "network_filtered")
    )
    write_tsv(class_counts,
              file.path(tables_dir,
                        "156_KO_prevalence_class_counts_after_filters.tsv"))

    class_n <- function(set, class) {
      z <- class_counts$n_KOs[
        class_counts$filter_set == set &
          class_counts$functional_prevalence_class == class
      ]
      if (length(z) == 0L) 0L else as.integer(z[1L])
    }
    validation <- data.frame(
      metric = c(
        "n_samples", "n_KOs_all_nonzero", "all_Prevalent", "all_Shell",
        "all_Cloud", "n_KOs_association_filtered", "association_Prevalent",
        "association_Shell", "association_Cloud", "n_KOs_network_filtered"
      ),
      observed = c(
        ncol(counts_all), nrow(counts_all),
        class_n("all_nonzero", "Prevalent"),
        class_n("all_nonzero", "Shell"),
        class_n("all_nonzero", "Cloud"),
        nrow(counts_assoc),
        class_n("association_filtered", "Prevalent"),
        class_n("association_filtered", "Shell"),
        class_n("association_filtered", "Cloud"),
        nrow(clr_network)
      ),
      expected = c(
        cfg$expected_n_samples, cfg$expected_n_kos,
        cfg$expected_all_prevalent, cfg$expected_all_shell,
        cfg$expected_all_cloud, cfg$expected_n_assoc,
        cfg$expected_assoc_prevalent, cfg$expected_assoc_shell,
        cfg$expected_assoc_cloud, cfg$expected_n_network
      ),
      stringsAsFactors = FALSE
    )
    validation$pass <- validation$observed == validation$expected
    write_tsv(validation,
              file.path(tables_dir, "156_expected_dimensions_validation.tsv"))

    transform_validation <- data.frame(
      metric = c(
        "max_relative_column_sum_error",
        "max_abs_CLR_column_mean_association",
        "n_nonfinite_CLR_association",
        "n_nonfinite_CLR_network"
      ),
      observed = c(
        max(abs(colSums(relative_assoc) - 1)),
        max(abs(colMeans(clr_assoc))),
        sum(!is.finite(clr_assoc)),
        sum(!is.finite(clr_network))
      ),
      tolerance = c(1e-10, 1e-10, 0, 0),
      stringsAsFactors = FALSE
    )
    transform_validation$pass <- c(
      transform_validation$observed[1L] <= transform_validation$tolerance[1L],
      transform_validation$observed[2L] <= transform_validation$tolerance[2L],
      transform_validation$observed[3L] == 0,
      transform_validation$observed[4L] == 0
    )
    write_tsv(transform_validation,
              file.path(tables_dir, "156_transformation_validation.tsv"))
    write_tsv(transform_validation,
              file.path(tables_dir, "156_CLR_transform_validation.tsv"))

    if (cfg$strict_expected && any(!validation$pass)) {
      failed <- validation$metric[!validation$pass]
      stopf("Manuscript dimension validation failed: %s",
            paste(failed, collapse = ", "))
    }
    if (any(!transform_validation$pass)) {
      stop("Relative-abundance or CLR validation failed.", call. = FALSE)
    }

    cat("[8/8] Saving provenance, objects and completion status...\n")
    parameters <- data.frame(
      parameter = c(
        "script", "base", "results_dir", "run134", "metadata_input",
        "metadata_id_column_source", "ko_matrix_input",
        "prevalence_annotation_reference", "presence_threshold",
        "prevalent_threshold", "shell_threshold", "min_prev_assoc",
        "min_total_assoc", "min_prev_network", "min_total_network",
        "max_network_kos", "pseudocount", "allow_relaxed_id_match",
        "missing_as_zero", "strict_reference", "strict_expected"
      ),
      value = c(
        script_name, base, results_dir, run134, metadata_path, id_col,
        ko_matrix_path, reference_path, cfg$presence_threshold,
        cfg$prevalent_threshold, cfg$shell_threshold, cfg$min_prev_assoc,
        cfg$min_total_assoc, cfg$min_prev_network, cfg$min_total_network,
        cfg$max_network_kos, cfg$pseudocount,
        cfg$allow_relaxed_id_match, cfg$missing_as_zero,
        cfg$strict_reference, cfg$strict_expected
      ),
      stringsAsFactors = FALSE
    )
    write_tsv(parameters, file.path(tables_dir, "156_parameters.tsv"))
    summary <- data.frame(
      item = c(
        "status", "run134", "ko_matrix", "prevalence_class_annotation",
        "canonical_sample_id", "metadata_run_accession_column",
        "n_metadata_samples", "n_input_matrix_columns",
        "n_unmatched_matrix_columns", "n_duplicate_KO_rows_collapsed",
        "n_missing_cells_converted_to_zero", "n_KOs_nonzero",
        "presence_threshold", "prevalent_threshold", "shell_threshold",
        "min_prev_assoc", "min_total_assoc", "n_KOs_association_filtered",
        "min_prev_network", "min_total_network", "max_network_kos",
        "n_KOs_network_candidates_before_cap", "n_KOs_network_filtered",
        "pseudocount", "network_CLR_reference", "association_CLR_matrix",
        "network_CLR_matrix", "metadata_aligned", "notes"
      ),
      value = c(
        "PASS", run134, ko_matrix_path,
        file.path(tables_dir, "156_KO_prevalence_class_annotation.tsv"),
        "run_accession", id_col, nrow(metadata), ncol(counts_raw),
        length(unmatched_matrix), ko_res$duplicated_rows_collapsed,
        sum(ko_res$conversion_audit$n_missing_converted_to_zero),
        nrow(counts_all), cfg$presence_threshold, cfg$prevalent_threshold,
        cfg$shell_threshold, cfg$min_prev_assoc, cfg$min_total_assoc,
        nrow(counts_assoc), cfg$min_prev_network, cfg$min_total_network,
        cfg$max_network_kos, length(candidate_KOs), nrow(clr_network),
        cfg$pseudocount, "association_filtered_feature_space",
        file.path(tables_dir, "156_KO_CLR_association_filtered.tsv.gz"),
        file.path(tables_dir, "156_KO_CLR_network_filtered.tsv.gz"),
        file.path(tables_dir, "156_metadata_aligned.tsv"),
        paste(
          "Prevalent-Shell-Cloud classes computed directly from the aligned",
          "raw matrix; classes are annotations, not filtering rules."
        )
      ),
      stringsAsFactors = FALSE
    )
    write_tsv(summary, file.path(tables_dir, "156_summary.tsv"))
    saveRDS(
      list(
        metadata = metadata,
        counts_all_nonzero = counts_all,
        prevalence_annotation = prevalence_annotation,
        counts_association = counts_assoc,
        relative_association = relative_assoc,
        clr_association = clr_assoc,
        counts_network = counts_network,
        clr_network = clr_network,
        parameters = cfg
      ),
      file.path(rds_dir, "156_KO_CLR_matrix_metadata_objects.rds")
    )
    saveRDS(clr_assoc,
            file.path(rds_dir, "156_KO_CLR_association_filtered.rds"))
    saveRDS(clr_network,
            file.path(rds_dir, "156_KO_CLR_network_filtered.rds"))
    writeLines("PASS", file.path(outdir, "STATUS.txt"))
    writeLines(c(
      "status\tPASS",
      paste0("completed_utc\t", format(Sys.time(), tz = "UTC", usetz = TRUE)),
      paste0("output_directory\t", outdir)
    ), file.path(logs_dir, "156_status.tsv"))
    cat("PASS: ", outdir, "\n", sep = "")
  }

  tryCatch(
    execute(),
    error = function(e) {
      writeLines(c(
        "status\tFAIL",
        paste0("failed_utc\t", format(Sys.time(), tz = "UTC", usetz = TRUE)),
        paste0("message\t", conditionMessage(e)),
        paste0("output_directory\t", outdir)
      ), file.path(logs_dir, "156_status.tsv"))
      writeLines("FAIL", file.path(outdir, "STATUS.txt"))
      stop(e)
    }
  )
  invisible(outdir)
}

main()
