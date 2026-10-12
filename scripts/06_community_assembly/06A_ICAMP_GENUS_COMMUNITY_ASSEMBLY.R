#!/usr/bin/env Rscript

# =============================================================================
# 06A_ICAMP_GENUS_COMMUNITY_ASSEMBLY.R
# Estimate bacterial genus-level assembly contributions using iCAMP.
#
# Inputs
# ------
# data/icamp_taxonomic_counts.rds: the separate iCAMP input with run_accession, REALM,
# EcologicalClass and genus annotations.
# Its embedded metadata are the metadata used in the reported analysis.
# The retained result used 400 metagenomes, 1,892 genera and 999 randomizations;
# it was not rerun with the taxonomic diversity/core cohort.
# data/genus_phylogeny_GTDB_R226_bac120.newick: external bacterial genus tree
# derived from the GTDB R226 bac120 reference phylogeny, with branch lengths.
# The reference file used for the manuscript has 1,905 tips and SHA-256:
# 15761713cbfbe4edcf0600bb65afbea0f39cfe6a01a7b111a92aeeba59ffd3fd
# The tree is an independent input, not a tree inferred from Kraken assignments.
# Supply a different external Newick with --tree PATH when needed.
#
# Name matching removes genus prefixes, punctuation and case differences.
# Ambiguous tree-tip keys are excluded. This is name harmonization, not a
# validated taxonomic synonym crosswalk or a test of phylogenetic niche signal.
# Inspect tables/08B_genus_tree_matching.tsv and 08B_tree_tip_key_audit.tsv.
#
# Primary settings
# ----------------
# EcologicalClass == "Mangrove sediment"; metacommunities defined by REALM;
# genus prevalence >=10 retained samples; no abundance rank cap; bMNTD with RC;
# 999 randomizations; minimum bin size 24; seed 42.
# Tree coverage is the fraction of all genus-assigned counts matched to
# unambiguous tree tips, calculated BEFORE genus-prevalence filtering. Samples
# with coverage <0.50 are excluded. Coverage is not the classified-read fraction
# of a complete metagenome. Outputs separate within-realm and between-realm pairs;
# the manuscript reports the within-realm summary.
#
# --smoke_test 1 selects up to three high-coverage samples per realm and up to
# 60 genera (prevalence >=3), with five randomizations and minimum bin size 12.
# --preflight_only 1 prepares and validates inputs without running the null model.
#
# Main outputs
# ------------
# tables/08B_icamp_pairwise_processes.tsv.gz
# tables/08B_icamp_process_summary_within_realm.tsv
# tables/08B_icamp_process_summary_between_realm_pair.tsv
# tables/08B_tree_coverage_by_sample.tsv
# tables/08B_genus_tree_matching.tsv
# tables/08B_parameters.tsv
# rds/08B_icamp_cm_result_genus.rds
# rds/08B_prepared_icamp_inputs.rds
# logs/08B_run.log; STATUS.txt
# =============================================================================

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)

SCRIPT_VERSION <- "2.3.0"
ICAMP_PROCESS_NAMES <- c(
  "heterogeneous_selection", "homogeneous_selection",
  "dispersal_limitation", "homogenizing_dispersal",
  "drift_weak_selection"
)

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

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
  missing_like <- is.na(y) | y == "" |
    tolower(y) %in% c("na", "n/a", "nan", "null", "none", "unknown")
  y[missing_like] <- NA_character_
  y
}

stopf <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)
warnf <- function(fmt, ...) warning(sprintf(fmt, ...), call. = FALSE)

sha256_file <- function(path) {
  executable <- Sys.which("sha256sum")
  if (!nzchar(executable) || is.na(path) || !file.exists(path)) {
    return(NA_character_)
  }
  output <- suppressWarnings(system2(executable, shQuote(path), stdout = TRUE,
                                     stderr = TRUE))
  if (length(output) == 0L) return(NA_character_)
  value <- strsplit(output[1L], "[[:space:]]+")[[1L]][1L]
  if (!grepl("^[0-9a-fA-F]{64}$", value)) NA_character_ else tolower(value)
}

normalise_name <- function(x) {
  tolower(gsub("[^a-zA-Z0-9]+", "", as.character(x)))
}

find_column <- function(x, candidates, label, required = TRUE) {
  direct <- candidates[candidates %in% names(x)]
  if (length(direct) > 0L) return(direct[1L])
  idx <- match(normalise_name(candidates), normalise_name(names(x)))
  idx <- idx[!is.na(idx)]
  if (length(idx) > 0L) return(names(x)[idx[1L]])
  if (required) {
    stopf("Could not identify %s. Tried: %s",
          label, paste(candidates, collapse = ", "))
  }
  NA_character_
}

write_tsv <- function(x, path) {
  utils::write.table(
    x, file = path, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE, na = ""
  )
}

write_tsv_gz <- function(x, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  utils::write.table(
    x, file = con, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE, na = ""
  )
}

write_matrix_tsv_gz <- function(x, path, row_id = "run_accession") {
  out <- data.frame(row_key = rownames(x), x, check.names = FALSE)
  names(out)[1L] <- row_id
  write_tsv_gz(out, path)
}

write_status <- function(outdir, status, message = NULL) {
  lines <- status
  if (!is.null(message) && nzchar(message)) {
    lines <- c(lines, paste0("message=", message))
  }
  lines <- c(lines, paste0("timestamp_utc=", format(
    Sys.time(), tz = "UTC", usetz = TRUE
  )))
  writeLines(lines, file.path(outdir, "STATUS.txt"))
}

package_versions <- function(pkgs) {
  data.frame(
    package = pkgs,
    version = vapply(pkgs, function(pkg) {
      if (requireNamespace(pkg, quietly = TRUE)) {
        as.character(utils::packageVersion(pkg))
      } else {
        NA_character_
      }
    }, character(1)),
    stringsAsFactors = FALSE
  )
}

strip_genus_prefix <- function(x) {
  y <- clean_chr(x)
  y <- gsub("^[\"']|[\"']$", "", y)

  has_gtdb_genus <- !is.na(y) & grepl("g__", y, fixed = TRUE)
  if (any(has_gtdb_genus)) {
    y[has_gtdb_genus] <- sub(".*g__", "", y[has_gtdb_genus])
    y[has_gtdb_genus] <- sub("[;|].*$", "", y[has_gtdb_genus])
  }

  y <- sub("^[gG]__", "", y)
  y <- trimws(y)
  bad <- is.na(y) | y == "" |
    grepl("^(unclassified|uncultured|unidentified|unknown|ambiguous|metagenome)$",
          y, ignore.case = TRUE) |
    grepl("^(bacteria|archaea)$", y, ignore.case = TRUE)
  y[bad] <- NA_character_
  y
}

genus_key <- function(x) {
  y <- strip_genus_prefix(x)
  key <- tolower(gsub("[^a-zA-Z0-9]+", "", y))
  key[is.na(y) | key == ""] <- NA_character_
  key
}

collapse_unique <- function(x, sep = "; ") {
  x <- sort(unique(clean_chr(x)))
  x <- x[!is.na(x)]
  if (length(x) == 0L) NA_character_ else paste(x, collapse = sep)
}

resolve_phyloseq_path <- function(requested) {
  if (!file.exists(requested)) {
    stopf("Phyloseq input does not exist: %s. Supply --phyloseq PATH or ICAMP_PHYLOSEQ.",
          requested)
  }
  normalizePath(requested, winslash = "/", mustWork = TRUE)
}

resolve_tree <- function(tree_arg) {
  if (is.na(tree_arg) || !nzchar(tree_arg) || !file.exists(tree_arg)) {
    stop("An external genus-level Newick tree is required. Supply --tree PATH or TREE_FILE.",
         call. = FALSE)
  }
  list(
    tree = ape::read.tree(tree_arg),
    source = "external_newick_argument",
    path = normalizePath(tree_arg, winslash = "/", mustWork = TRUE)
  )
}

validate_tree <- function(tree) {
  if (!inherits(tree, "phylo")) stop("The selected tree is not class 'phylo'.", call. = FALSE)
  if (length(tree$tip.label) < 10L) stop("The selected tree has fewer than 10 tips.", call. = FALSE)
  if (anyNA(tree$tip.label) || anyDuplicated(tree$tip.label)) {
    stop("Tree tip labels must be non-missing and unique.", call. = FALSE)
  }
  if (is.null(tree$edge.length)) {
    stop("The genus tree lacks branch lengths; iCAMP requires phylogenetic distances.",
         call. = FALSE)
  }
  if (any(!is.finite(tree$edge.length)) || any(tree$edge.length < 0)) {
    stop("Tree branch lengths must be finite and non-negative.", call. = FALSE)
  }
  if (!any(tree$edge.length > 0)) {
    stop("All tree branch lengths are zero.", call. = FALSE)
  }
  invisible(TRUE)
}

aggregate_to_genus <- function(ps) {
  taxonomy <- as(phyloseq::tax_table(ps), "matrix")
  genus_col <- find_column(
    as.data.frame(taxonomy, check.names = FALSE),
    c("Genus", "genus", "Rank_G", "G", "g"),
    "genus column in tax_table"
  )
  raw_genus <- taxonomy[, genus_col]
  cleaned_genus <- strip_genus_prefix(raw_genus)
  key <- genus_key(raw_genus)

  otu <- as(phyloseq::otu_table(ps), "matrix")
  if (!phyloseq::taxa_are_rows(ps)) otu <- t(otu)
  if (!identical(rownames(otu), rownames(taxonomy))) {
    taxonomy <- taxonomy[match(rownames(otu), rownames(taxonomy)), , drop = FALSE]
    raw_genus <- taxonomy[, genus_col]
    cleaned_genus <- strip_genus_prefix(raw_genus)
    key <- genus_key(raw_genus)
  }

  valid <- !is.na(key)
  if (sum(valid) < 10L) {
    stopf("Only %d taxonomic features have a usable genus annotation.", sum(valid))
  }

  counts_by_key <- rowsum(otu[valid, , drop = FALSE], group = key[valid],
                          reorder = TRUE, na.rm = FALSE)
  labels_by_key <- split(cleaned_genus[valid], key[valid])
  source_features_by_key <- split(rownames(otu)[valid], key[valid])
  genus_map <- data.frame(
    genus_key = names(labels_by_key),
    genus_label = vapply(labels_by_key, collapse_unique, character(1)),
    n_label_variants = vapply(labels_by_key, function(z) {
      length(unique(stats::na.omit(clean_chr(z))))
    }, integer(1)),
    n_source_features = vapply(source_features_by_key, length, integer(1)),
    stringsAsFactors = FALSE
  )
  genus_map <- genus_map[match(rownames(counts_by_key), genus_map$genus_key),,
                         drop = FALSE]

  comm <- t(counts_by_key)
  storage.mode(comm) <- "numeric"
  list(
    comm = comm,
    genus_map = genus_map,
    genus_column = genus_col,
    n_input_features = nrow(otu),
    n_features_without_genus = sum(!valid)
  )
}

match_genera_to_tree <- function(comm, genus_map, tree) {
  tip_map <- data.frame(
    tree_tip = tree$tip.label,
    genus_key = genus_key(tree$tip.label),
    stringsAsFactors = FALSE
  )
  tip_map$key_ambiguous <- duplicated(tip_map$genus_key) |
    duplicated(tip_map$genus_key, fromLast = TRUE) |
    is.na(tip_map$genus_key)

  unique_tips <- tip_map[!tip_map$key_ambiguous, , drop = FALSE]
  match_index <- match(genus_map$genus_key, unique_tips$genus_key)
  genus_map$tree_tip <- unique_tips$tree_tip[match_index]
  genus_map$tree_match <- !is.na(genus_map$tree_tip)
  genus_map$total_count_before_filters <- colSums(comm)[genus_map$genus_key]
  genus_map$prevalence_n_before_filters <- colSums(comm > 0)[genus_map$genus_key]

  list(genus_map = genus_map, tip_map = tip_map)
}

select_smoke_samples <- function(metadata, comm, matched_keys,
                                 group_col, n_per_group) {
  matched_keys <- intersect(colnames(comm), matched_keys)
  if (length(matched_keys) == 0L) {
    stop("No tree-matched genera are available for smoke sample selection.",
         call. = FALSE)
  }
  assigned_count <- rowSums(comm)
  matched_count <- rowSums(comm[, matched_keys, drop = FALSE])
  matched_fraction <- matched_count / assigned_count
  matched_fraction[!is.finite(matched_fraction)] <- -Inf
  selected <- unlist(lapply(split(seq_len(nrow(metadata)), metadata[[group_col]]),
    function(index) {
      ordered <- index[order(
        -matched_fraction[index],
        -matched_count[index],
        metadata$run_accession[index]
      )]
      utils::head(ordered, n_per_group)
    }
  ), use.names = FALSE)
  selected <- sort(unique(selected))
  metadata[selected, , drop = FALSE]
}

detect_process_columns <- function(x) {
  aliases <- list(
    heterogeneous_selection = c("hes", "heterogeneousselection", "heterosel"),
    homogeneous_selection = c("hos", "homogeneousselection", "homosel"),
    dispersal_limitation = c("dl", "dispersallimitation", "dispersallimited"),
    homogenizing_dispersal = c("hd", "homogenizingdispersal", "homogenisingdispersal"),
    drift_weak_selection = c("dr", "drift", "driftandothers", "driftorweakselection")
  )
  normal <- normalise_name(names(x))
  found <- vapply(aliases, function(values) {
    hit <- which(normal %in% values)
    if (length(hit) == 0L) NA_character_ else names(x)[hit[1L]]
  }, character(1))
  if (anyNA(found)) return(NULL)
  found
}

extract_pairwise_processes <- function(icamp_result) {
  candidates <- list()
  if (!is.null(icamp_result$qp.comm)) {
    candidates[["qp.comm"]] <- icamp_result$qp.comm
  }
  top_names <- names(icamp_result)
  if (is.null(top_names)) top_names <- paste0("element_", seq_along(icamp_result))
  for (i in seq_along(icamp_result)) {
    object <- icamp_result[[i]]
    if (is.matrix(object) || is.data.frame(object)) {
      candidates[[top_names[i]]] <- object
    }
  }

  for (candidate_name in names(candidates)) {
    x <- as.data.frame(candidates[[candidate_name]], check.names = FALSE,
                       stringsAsFactors = FALSE)
    process_columns <- detect_process_columns(x)
    if (is.null(process_columns) || ncol(x) < 7L) next

    id_columns <- setdiff(names(x), unname(process_columns))
    if (length(id_columns) < 2L) next
    out <- data.frame(
      sample1 = clean_chr(x[[id_columns[1L]]]),
      sample2 = clean_chr(x[[id_columns[2L]]]),
      heterogeneous_selection = suppressWarnings(as.numeric(
        x[[process_columns["heterogeneous_selection"]]]
      )),
      homogeneous_selection = suppressWarnings(as.numeric(
        x[[process_columns["homogeneous_selection"]]]
      )),
      dispersal_limitation = suppressWarnings(as.numeric(
        x[[process_columns["dispersal_limitation"]]]
      )),
      homogenizing_dispersal = suppressWarnings(as.numeric(
        x[[process_columns["homogenizing_dispersal"]]]
      )),
      drift_weak_selection = suppressWarnings(as.numeric(
        x[[process_columns["drift_weak_selection"]]]
      )),
      stringsAsFactors = FALSE
    )
    process_names <- names(out)[3:7]
    sums <- rowSums(out[process_names])
    if (stats::median(sums, na.rm = TRUE) > 2) {
      out[process_names] <- out[process_names] / 100
      sums <- rowSums(out[process_names])
    }
    out$process_sum <- sums
    out$dominant_process <- process_names[max.col(
      as.matrix(out[process_names]), ties.method = "first"
    )]
    attr(out, "source_element") <- candidate_name
    return(out)
  }
  NULL
}

annotate_pairwise_geography <- function(pairwise, meta_group) {
  if (is.null(rownames(meta_group)) || ncol(meta_group) < 1L) {
    stop("meta_group must have sample IDs as row names and one group column.",
         call. = FALSE)
  }

  sample_ids <- clean_chr(rownames(meta_group))
  group_values <- clean_chr(meta_group[[1L]])
  if (anyNA(sample_ids) || anyDuplicated(sample_ids)) {
    stop("meta_group has missing or duplicated sample IDs.", call. = FALSE)
  }
  if (anyNA(group_values)) {
    stop("meta_group has missing metacommunity labels.", call. = FALSE)
  }

  group_map <- stats::setNames(group_values, sample_ids)
  turnover_ids <- unique(c(pairwise$sample1, pairwise$sample2))
  missing_ids <- setdiff(turnover_ids, names(group_map))
  if (length(missing_ids) > 0L) {
    stopf(
      "Pairwise iCAMP output contains sample IDs absent from meta_group: %s",
      paste(utils::head(missing_ids, 20L), collapse = ", ")
    )
  }

  realm1 <- unname(group_map[pairwise$sample1])
  realm2 <- unname(group_map[pairwise$sample2])
  first_realm <- ifelse(realm1 <= realm2, realm1, realm2)
  second_realm <- ifelse(realm1 <= realm2, realm2, realm1)

  geography <- data.frame(
    realm1 = realm1,
    realm2 = realm2,
    comparison_scope = ifelse(
      realm1 == realm2, "within_realm", "between_realms"
    ),
    realm_pair = paste(first_realm, second_realm, sep = " || "),
    stringsAsFactors = FALSE
  )

  remaining_columns <- setdiff(names(pairwise), c("sample1", "sample2"))
  cbind(
    pairwise[c("sample1", "sample2")],
    geography,
    pairwise[remaining_columns],
    stringsAsFactors = FALSE
  )
}

audit_pairwise_structure <- function(pairwise, meta_group) {
  sample_ids <- rownames(meta_group)
  group_values <- clean_chr(meta_group[[1L]])
  group_sizes <- table(group_values)
  group_levels <- sort(names(group_sizes))

  expected_list <- list()
  counter <- 0L
  for (i in seq_along(group_levels)) {
    for (j in i:length(group_levels)) {
      counter <- counter + 1L
      group1 <- group_levels[i]
      group2 <- group_levels[j]
      expected_n <- if (identical(group1, group2)) {
        choose(as.integer(group_sizes[group1]), 2)
      } else {
        as.integer(group_sizes[group1]) * as.integer(group_sizes[group2])
      }
      expected_list[[counter]] <- data.frame(
        realm_pair = paste(group1, group2, sep = " || "),
        comparison_scope = ifelse(
          identical(group1, group2), "within_realm", "between_realms"
        ),
        expected_pairs = as.integer(expected_n),
        stringsAsFactors = FALSE
      )
    }
  }
  expected_counts <- do.call(rbind, expected_list)

  observed_counts <- as.data.frame(
    table(pairwise$realm_pair), stringsAsFactors = FALSE
  )
  names(observed_counts) <- c("realm_pair", "observed_pairs")
  observed_counts$observed_pairs <- as.integer(observed_counts$observed_pairs)

  pair_counts <- merge(
    expected_counts, observed_counts, by = "realm_pair",
    all = TRUE, sort = FALSE
  )
  pair_counts$expected_pairs[is.na(pair_counts$expected_pairs)] <- 0L
  pair_counts$observed_pairs[is.na(pair_counts$observed_pairs)] <- 0L
  pair_counts <- pair_counts[
    order(pair_counts$comparison_scope, pair_counts$realm_pair),
    , drop = FALSE
  ]
  rownames(pair_counts) <- NULL

  expected_total <- choose(length(sample_ids), 2)
  expected_within <- sum(choose(as.integer(group_sizes), 2))
  expected_between <- expected_total - expected_within
  observed_within <- sum(pairwise$comparison_scope == "within_realm")
  observed_between <- sum(pairwise$comparison_scope == "between_realms")

  pair_keys <- vapply(seq_len(nrow(pairwise)), function(i) {
    paste(sort(c(pairwise$sample1[i], pairwise$sample2[i])), collapse = "||")
  }, character(1))

  process_matrix <- as.matrix(pairwise[ICAMP_PROCESS_NAMES])
  process_sums <- rowSums(process_matrix)
  all_process_values_finite <- all(is.finite(process_matrix))
  process_values_in_range <- all_process_values_finite &&
    all(process_matrix >= -1e-12 & process_matrix <= 1 + 1e-12)
  maximum_sum_error <- if (length(process_sums) > 0L) {
    suppressWarnings(max(abs(process_sums - 1), na.rm = TRUE))
  } else {
    Inf
  }

  validation <- data.frame(
    check = c(
      "all_samples_present", "all_sample_pairs_present",
      "within_realm_pair_count", "between_realm_pair_count",
      "all_realm_pair_counts_exact", "all_pair_scopes_classified",
      "no_self_comparisons", "no_duplicate_unordered_pairs",
      "all_process_values_finite", "process_values_between_zero_and_one",
      "five_processes_sum_to_one"
    ),
    pass = c(
      all(c(pairwise$sample1, pairwise$sample2) %in% sample_ids),
      nrow(pairwise) == expected_total,
      observed_within == expected_within,
      observed_between == expected_between,
      all(pair_counts$expected_pairs == pair_counts$observed_pairs),
      all(pairwise$comparison_scope %in% c("within_realm", "between_realms")),
      all(pairwise$sample1 != pairwise$sample2),
      !anyDuplicated(pair_keys),
      all_process_values_finite,
      process_values_in_range,
      is.finite(maximum_sum_error) && maximum_sum_error <= 1e-6
    ),
    observed = as.character(c(
      length(unique(c(pairwise$sample1, pairwise$sample2))),
      nrow(pairwise), observed_within, observed_between,
      sum(pair_counts$expected_pairs == pair_counts$observed_pairs),
      sum(pairwise$comparison_scope %in% c("within_realm", "between_realms")),
      sum(pairwise$sample1 == pairwise$sample2),
      anyDuplicated(pair_keys), all_process_values_finite,
      paste(range(process_matrix, na.rm = TRUE), collapse = " to "),
      maximum_sum_error
    )),
    expected = c(
      as.character(length(sample_ids)), as.character(expected_total),
      as.character(expected_within), as.character(expected_between),
      as.character(nrow(pair_counts)), as.character(nrow(pairwise)),
      "0", "0", "TRUE", "0 to 1", "<= 1e-6"
    ),
    stringsAsFactors = FALSE
  )

  list(
    validation = validation,
    pair_counts = pair_counts,
    expected_total = as.integer(expected_total),
    expected_within = as.integer(expected_within),
    expected_between = as.integer(expected_between),
    observed_within = as.integer(observed_within),
    observed_between = as.integer(observed_between)
  )
}

summarise_processes <- function(pairwise, group_cols = character()) {
  missing_processes <- setdiff(ICAMP_PROCESS_NAMES, names(pairwise))
  if (length(missing_processes) > 0L) {
    stopf("Missing iCAMP process columns: %s",
          paste(missing_processes, collapse = ", "))
  }
  missing_groups <- setdiff(group_cols, names(pairwise))
  if (length(missing_groups) > 0L) {
    stopf("Missing process-summary grouping columns: %s",
          paste(missing_groups, collapse = ", "))
  }
  if (nrow(pairwise) == 0L) {
    stop("Cannot summarize an empty pairwise process table.", call. = FALSE)
  }

  if (length(group_cols) == 0L) {
    split_rows <- list(all_pairs = seq_len(nrow(pairwise)))
  } else {
    group_key <- do.call(
      paste,
      c(lapply(pairwise[group_cols], as.character), sep = "\034")
    )
    split_rows <- split(seq_len(nrow(pairwise)), group_key, drop = TRUE)
  }

  pieces <- lapply(split_rows, function(index) {
    statistics <- do.call(rbind, lapply(ICAMP_PROCESS_NAMES, function(process) {
      values <- pairwise[[process]][index]
      values <- values[is.finite(values)]
      data.frame(
        process = process,
        n_turnovers = length(values),
        mean_fraction = mean(values),
        sd_fraction = stats::sd(values),
        median_fraction = stats::median(values),
        q025_fraction = as.numeric(
          stats::quantile(values, 0.025, names = FALSE)
        ),
        q975_fraction = as.numeric(
          stats::quantile(values, 0.975, names = FALSE)
        ),
        mean_percent = 100 * mean(values),
        sd_percent = 100 * stats::sd(values),
        stringsAsFactors = FALSE
      )
    }))

    if (length(group_cols) == 0L) return(statistics)
    group_values <- pairwise[index[1L], group_cols, drop = FALSE]
    cbind(
      group_values[rep(1L, nrow(statistics)), , drop = FALSE],
      statistics,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}

summarise_dominant_process <- function(pairwise, group_cols = character()) {
  missing_groups <- setdiff(group_cols, names(pairwise))
  if (length(missing_groups) > 0L) {
    stopf("Missing dominant-process grouping columns: %s",
          paste(missing_groups, collapse = ", "))
  }
  if (nrow(pairwise) == 0L) {
    stop("Cannot summarize an empty dominant-process table.", call. = FALSE)
  }

  if (length(group_cols) == 0L) {
    split_rows <- list(all_pairs = seq_len(nrow(pairwise)))
  } else {
    group_key <- do.call(
      paste,
      c(lapply(pairwise[group_cols], as.character), sep = "\034")
    )
    split_rows <- split(seq_len(nrow(pairwise)), group_key, drop = TRUE)
  }

  pieces <- lapply(split_rows, function(index) {
    counts <- table(pairwise$dominant_process[index])
    result <- data.frame(
      dominant_process = names(counts),
      n_turnovers = as.integer(counts),
      proportion = as.integer(counts) / sum(counts),
      stringsAsFactors = FALSE
    )
    if (length(group_cols) == 0L) return(result)
    group_values <- pairwise[index[1L], group_cols, drop = FALSE]
    cbind(
      group_values[rep(1L, nrow(result)), , drop = FALSE],
      result,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}

print_help <- function() {
  cat(paste(
    "Usage:",
    "  Rscript scripts/06_community_assembly/06A_ICAMP_GENUS_COMMUNITY_ASSEMBLY.R --smoke_test 1",
    "",
    "Main options:",
    "  --phyloseq PATH          Default ICAMP_PHYLOSEQ or BASE/data/icamp_taxonomic_counts.rds",
    "  --tree PATH              Default BASE/data/genus_phylogeny_GTDB_R226_bac120.newick",
    "  --base PATH              Default BIOGEOGRAPHY_BASE or current directory",
    "  --results_dir PATH       Default BIOGEOGRAPHY_RESULTS or BASE/results",
    "  --meta_group_var NAME    Default REALM",
    "  --analysis_label TEXT    Optional filesystem-safe run label",
    "  --smoke_test 0|1         Default 0",
    "  --preflight_only 0|1     Prepare/validate inputs without running iCAMP",
    "  --rand INT               Full default 999",
    "  --n_top_genus INT        Full default 0 (all eligible genera)",
    "  --min_prev_n INT         Full default 10; sensitivity: 5 and 20",
    "  --bin_size_limit INT     Full default 24",
    "  --nworker INT            Full default min(20, detected cores - 1)",
    "  --sig_index NAME         Default RC; sensitivity alternative SES.RC",
    "  --min_sample_tree_coverage NUM  Default 0.50",
    "  --filter_samples_by_tree_coverage 0|1  Default 1",
    "  --allow_low_tree_coverage 0|1   Default 0",
    sep = "\n"
  ))
}

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

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

PS_REQUESTED <- get_arg(
  "--phyloseq",
  Sys.getenv("ICAMP_PHYLOSEQ", file.path(PROJECT_DIR, "data", "icamp_taxonomic_counts.rds"))
)
TREE_DEFAULT <- file.path(
  PROJECT_DIR, "data", "genus_phylogeny_GTDB_R226_bac120.newick"
)
CANONICAL_TREE_SHA256 <- "15761713cbfbe4edcf0600bb65afbea0f39cfe6a01a7b111a92aeeba59ffd3fd"
TREE_ARG <- get_arg("--tree", Sys.getenv("TREE_FILE", TREE_DEFAULT))

SMOKE_TEST <- as_flag(get_arg("--smoke_test", Sys.getenv("SMOKE_TEST", "0")))
PREFLIGHT_ONLY <- as_flag(get_arg(
  "--preflight_only", Sys.getenv("PREFLIGHT_ONLY", "0")
))
FILTER_ECOLOGICAL_CLASS <- as_flag(get_arg(
  "--filter_ecological_class", "1"
))
ECOLOGICAL_CLASS <- get_arg("--ecological_class", "Mangrove sediment")
META_GROUP_VAR_REQUESTED <- get_arg("--meta_group_var", "REALM")
DROP_UNASSIGNED_REALM <- as_flag(get_arg("--drop_unassigned_realm", "1"))
ALLOW_LOW_TREE_COVERAGE <- as_flag(get_arg(
  "--allow_low_tree_coverage", "0"
))
FILTER_SAMPLES_BY_TREE_COVERAGE <- as_flag(get_arg(
  "--filter_samples_by_tree_coverage", "1"
))
MIN_SAMPLE_TREE_COVERAGE <- as.numeric(get_arg(
  "--min_sample_tree_coverage", "0.50"
))

ANALYSIS_LABEL <- clean_chr(get_arg("--analysis_label", ""))
if (length(ANALYSIS_LABEL) == 0L || is.na(ANALYSIS_LABEL[1L])) {
  ANALYSIS_LABEL <- ""
} else {
  ANALYSIS_LABEL <- gsub(
    "[^A-Za-z0-9._-]+", "_", ANALYSIS_LABEL[1L]
  )
  ANALYSIS_LABEL <- gsub("^_+|_+$", "", ANALYSIS_LABEL)
  if (!nzchar(ANALYSIS_LABEL)) {
    stop("--analysis_label does not contain a usable label.", call. = FALSE)
  }
}

EXCLUDE_RUNS <- clean_chr(strsplit(get_arg(
  "--exclude_runs",
  paste(c("SRR29659508", "SRR29659599", "SRR29659629", "SRR29659854"),
        collapse = ",")
), "[,;]")[[1L]])
EXCLUDE_RUNS <- EXCLUDE_RUNS[!is.na(EXCLUDE_RUNS)]

detected_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
if (!is.finite(detected_cores)) detected_cores <- 2L
default_workers <- max(1L, min(20L, as.integer(detected_cores) - 1L))

if (SMOKE_TEST) {
  RAND <- as.integer(get_arg("--rand", "5"))
  N_TOP_GENUS <- as.integer(get_arg("--n_top_genus", "60"))
  MIN_PREV_N <- as.integer(get_arg("--min_prev_n", "3"))
  BIN_SIZE_LIMIT <- as.integer(get_arg("--bin_size_limit", "12"))
  NWORKER <- as.integer(get_arg("--nworker", "2"))
  MEMORY_G <- as.numeric(get_arg("--memory_g", "10"))
  SMOKE_N_PER_REALM <- as.integer(get_arg("--smoke_n_per_realm", "3"))
  DETAIL_SAVE <- as_flag(get_arg("--detail_save", "0"))
} else {
  RAND <- as.integer(get_arg("--rand", "999"))
  N_TOP_GENUS <- as.integer(get_arg("--n_top_genus", "0"))
  MIN_PREV_N <- as.integer(get_arg("--min_prev_n", "10"))
  BIN_SIZE_LIMIT <- as.integer(get_arg("--bin_size_limit", "24"))
  NWORKER <- as.integer(get_arg("--nworker", as.character(default_workers)))
  MEMORY_G <- as.numeric(get_arg("--memory_g", "50"))
  SMOKE_N_PER_REALM <- NA_integer_
  DETAIL_SAVE <- as_flag(get_arg("--detail_save", "1"))
}

SEED <- as.integer(get_arg("--seed", "42"))
SIG_INDEX <- get_arg("--sig_index", "RC")
PHYLO_METRIC <- get_arg("--phylo_metric", "bMNTD")
STRICT_OUTPUT <- as_flag(get_arg("--strict_output", "1"))

if (!is.finite(RAND) || RAND < 2L) stop("--rand must be >= 2.", call. = FALSE)
if (!is.finite(N_TOP_GENUS) || N_TOP_GENUS < 0L ||
    (N_TOP_GENUS > 0L && N_TOP_GENUS < 10L)) {
  stop("--n_top_genus must be 0 or >= 10.", call. = FALSE)
}
if (!is.finite(MIN_PREV_N) || MIN_PREV_N < 1L) {
  stop("--min_prev_n must be >= 1.", call. = FALSE)
}
if (!is.finite(BIN_SIZE_LIMIT) || BIN_SIZE_LIMIT < 5L) {
  stop("--bin_size_limit must be >= 5.", call. = FALSE)
}
if (!is.finite(NWORKER) || NWORKER < 1L) NWORKER <- 1L
if (!is.finite(MIN_SAMPLE_TREE_COVERAGE) ||
    MIN_SAMPLE_TREE_COVERAGE < 0 || MIN_SAMPLE_TREE_COVERAGE > 1) {
  stop("--min_sample_tree_coverage must be between 0 and 1.",
       call. = FALSE)
}
if (!SIG_INDEX %in% c("Confidence", "SES.RC", "SES", "RC")) {
  stop("--sig_index must be Confidence, SES.RC, SES, or RC.", call. = FALSE)
}
if (!PHYLO_METRIC %in% c("bMPD", "bMNTD", "both", "bNRI", "bNTI")) {
  stop("Unsupported --phylo_metric.", call. = FALSE)
}

# Prevent BLAS x iCAMP worker oversubscription.
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  BLIS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)
set.seed(SEED)

SCRIPT_ID <- "08B_ICAMP_GENUS_COMMUNITY_ASSEMBLY"
STAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
RUN_LABEL_PARTS <- c(
  SCRIPT_ID,
  if (SMOKE_TEST) "SMOKE" else character(0),
  if (nzchar(ANALYSIS_LABEL)) ANALYSIS_LABEL else character(0)
)
RUN_LABEL <- paste(RUN_LABEL_PARTS, collapse = "_")
OUT_DIR <- file.path(RESULTS_DIR, paste0(RUN_LABEL, "_", STAMP))
TABLES_DIR <- file.path(OUT_DIR, "tables")
RDS_DIR <- file.path(OUT_DIR, "rds")
LOGS_DIR <- file.path(OUT_DIR, "logs")
PD_WD <- file.path(OUT_DIR, "pd_big")

dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOGS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PD_WD, recursive = TRUE, showWarnings = FALSE)

LATEST_FILE <- file.path(RESULTS_DIR, paste0("LATEST_", RUN_LABEL, ".txt"))
writeLines(OUT_DIR, LATEST_FILE)
write_status(OUT_DIR, "RUNNING")

# -----------------------------------------------------------------------------
# Main analysis
# -----------------------------------------------------------------------------

status_ok <- FALSE
log_connection <- file(file.path(LOGS_DIR, "08B_run.log"), open = "wt")
sink(log_connection, split = TRUE)
sink(log_connection, type = "message")

main <- function() {
  cat("============================================================\n")
  cat(SCRIPT_ID, "\n")
  cat("Mode: ", ifelse(SMOKE_TEST, "SMOKE", "FULL"), "\n", sep = "")
  cat("Analysis label: ",
      ifelse(nzchar(ANALYSIS_LABEL), ANALYSIS_LABEL, "none"),
      "\n", sep = "")
  cat("Started: ", format(Sys.time(), tz = "UTC", usetz = TRUE), "\n", sep = "")
  cat("Output: ", OUT_DIR, "\n", sep = "")
  cat("============================================================\n\n")

  required_packages <- c("phyloseq", "ape", "iCAMP")
  versions <- package_versions(c(required_packages, "parallel"))
  write_tsv(versions, file.path(TABLES_DIR, "08B_package_versions.tsv"))
  missing_packages <- versions$package[is.na(versions$version) &
    versions$package %in% required_packages]
  if (length(missing_packages) > 0L) {
    stopf("Missing required R package(s): %s",
          paste(missing_packages, collapse = ", "))
  }
  icamp_version <- utils::packageVersion("iCAMP")
  icamp_cm_formals <- names(formals(iCAMP::icamp.cm))
  required_icamp_cm_formals <- c(
    "comm", "tree", "meta.group", "pd.desc", "pd.spname", "pd.wd",
    "rand", "prefix", "phylo.rand.scale", "taxa.rand.scale",
    "phylo.metric", "sig.index", "bin.size.limit", "nworker",
    "memory.G", "rtree.save", "detail.save", "qp.save", "detail.null",
    "ignore.zero", "output.wd", "correct.special", "unit.sum",
    "special.method", "ses.cut", "rc.cut", "conf.cut", "omit.option",
    "taxo.metric", "transform.method", "dirichlet"
  )
  missing_icamp_cm_formals <- setdiff(
    required_icamp_cm_formals, icamp_cm_formals
  )
  if (length(missing_icamp_cm_formals) > 0L) {
    stopf(
      "Installed iCAMP::icamp.cm lacks required argument(s): %s",
      paste(missing_icamp_cm_formals, collapse = ", ")
    )
  }
  if (icamp_version < "1.8.2") {
    warnf(
      paste(
        "iCAMP %s is older than 1.8.2. The official changelog records",
        "fixes in 1.8.1 for Confidence calculations and in 1.8.2 for",
        "functions that may affect icamp.cm. This run uses sig.index=%s;",
        "the installed version is recorded."
      ),
      as.character(icamp_version), SIG_INDEX
    )
    if (!PREFLIGHT_ONLY) {
      stop(paste(
        "Executable iCAMP runs with multiple metacommunities require",
        "iCAMP >= 1.8.2. Version 1.5.12 contains an RC.bin.cm scoping bug",
        "that attempts to use an undefined object named 'comm' when",
        "meta.group is supplied. Preflight-only mode remains available."
      ), call. = FALSE)
    }
  }

  cat("[1/9] Resolving and reading taxonomic counts...\n")
  ps_path <- resolve_phyloseq_path(PS_REQUESTED)
  ps <- readRDS(ps_path)
  if (!inherits(ps, "phyloseq")) {
    stop("Input RDS is not a phyloseq object.", call. = FALSE)
  }
  cat("    Input: ", ps_path, "\n", sep = "")
  cat("    Initial samples: ", phyloseq::nsamples(ps), "\n", sep = "")
  cat("    Initial taxonomic features: ", phyloseq::ntaxa(ps), "\n", sep = "")

  metadata <- as(phyloseq::sample_data(ps), "data.frame")
  metadata$.phyloseq_sample_id <- rownames(metadata)
  run_col <- find_column(
    metadata,
    c("run_accession", ".sample_name", "sample_name", "Run"),
    "canonical run_accession column"
  )
  metadata$run_accession <- clean_chr(metadata[[run_col]])
  if (anyNA(metadata$run_accession) || anyDuplicated(metadata$run_accession)) {
    stop("run_accession is missing or duplicated in sample_data.", call. = FALSE)
  }

  cat("[2/9] Applying ecological-class and accession filters...\n")
  keep <- rep(TRUE, nrow(metadata))
  filter_reason <- rep("retained", nrow(metadata))
  if (FILTER_ECOLOGICAL_CLASS) {
    ecological_col <- find_column(
      metadata, c("EcologicalClass", "ecological_class"),
      "EcologicalClass column"
    )
    ecological_value <- clean_chr(metadata[[ecological_col]])
    wrong_class <- is.na(ecological_value) | ecological_value != ECOLOGICAL_CLASS
    keep[wrong_class] <- FALSE
    filter_reason[wrong_class] <- "not_mangrove_sediment"
  }
  excluded_accession <- metadata$run_accession %in% EXCLUDE_RUNS
  keep[excluded_accession] <- FALSE
  filter_reason[excluded_accession] <- "excluded_accession"

  sample_audit <- data.frame(
    phyloseq_sample_id = metadata$.phyloseq_sample_id,
    run_accession = metadata$run_accession,
    retained_after_canonical_filters = keep,
    filter_reason = filter_reason,
    stringsAsFactors = FALSE
  )
  write_tsv(sample_audit, file.path(TABLES_DIR, "08B_sample_filter_audit.tsv"))
  if (sum(keep) < 10L) stop("Fewer than 10 samples remain after filtering.", call. = FALSE)

  retained_ps_ids <- metadata$.phyloseq_sample_id[keep]
  ps <- phyloseq::prune_samples(retained_ps_ids, ps)
  ps <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps) > 0, ps)
  metadata <- metadata[keep, , drop = FALSE]
  metadata <- metadata[match(phyloseq::sample_names(ps), metadata$.phyloseq_sample_id),,
                       drop = FALSE]

  group_col <- find_column(
    metadata,
    c(META_GROUP_VAR_REQUESTED, "REALM", "marine_realm", "realm"),
    "marine-realm metacommunity column"
  )
  metadata$meta_group <- clean_chr(metadata[[group_col]])
  missing_group <- is.na(metadata$meta_group) |
    tolower(metadata$meta_group) %in% c("unassigned", "unknown")
  if (any(missing_group)) {
    if (DROP_UNASSIGNED_REALM) {
      ps <- phyloseq::prune_samples(
        metadata$.phyloseq_sample_id[!missing_group], ps
      )
      metadata <- metadata[!missing_group, , drop = FALSE]
    } else {
      metadata$meta_group[missing_group] <- "Unassigned"
    }
  }
  if (nrow(metadata) < 10L) stop("Fewer than 10 realm-assigned samples remain.", call. = FALSE)

  cat("    Samples after ecological-class and accession filters: ", nrow(metadata), "\n", sep = "")
  cat("    Metacommunities: ",
      paste(sort(unique(metadata$meta_group)), collapse = "; "), "\n", sep = "")

  cat("[3/9] Aggregating taxonomic counts to genus...\n")
  genus_result <- aggregate_to_genus(ps)
  comm <- genus_result$comm
  rownames(comm) <- metadata$run_accession[
    match(rownames(comm), metadata$.phyloseq_sample_id)
  ]
  if (anyNA(rownames(comm)) || anyDuplicated(rownames(comm))) {
    stop("Failed to assign unique run_accession row names to genus matrix.",
         call. = FALSE)
  }
  metadata <- metadata[match(rownames(comm), metadata$run_accession),,
                       drop = FALSE]
  if (!identical(rownames(comm), metadata$run_accession)) {
    stop("Metadata and genus matrix are not identically ordered.", call. = FALSE)
  }
  if (any(!is.finite(comm)) || any(comm < 0)) {
    stop("Genus count matrix contains negative or non-finite values.", call. = FALSE)
  }

  cat("[4/9] Resolving genus phylogeny and matching genus identifiers...\n")
  tree_result <- resolve_tree(TREE_ARG)
  tree <- tree_result$tree
  validate_tree(tree)
  tree_sha256 <- sha256_file(tree_result$path)
  canonical_tree_in_use <- !is.na(tree_result$path) &&
    file.exists(TREE_DEFAULT) &&
    identical(
      normalizePath(tree_result$path, winslash = "/", mustWork = TRUE),
      normalizePath(TREE_DEFAULT, winslash = "/", mustWork = TRUE)
    )
  if (canonical_tree_in_use &&
      (is.na(tree_sha256) || tree_sha256 != CANONICAL_TREE_SHA256)) {
    stopf(
      "GTDB reference genus tree SHA-256 mismatch. Expected %s; observed %s.",
      CANONICAL_TREE_SHA256,
      ifelse(is.na(tree_sha256), "unavailable", tree_sha256)
    )
  }
  matching <- match_genera_to_tree(comm, genus_result$genus_map, tree)
  genus_map <- matching$genus_map
  write_tsv(matching$tip_map,
            file.path(TABLES_DIR, "08B_tree_tip_key_audit.tsv"))

  matched_keys <- genus_map$genus_key[genus_map$tree_match]
  matched_keys <- intersect(colnames(comm), matched_keys)
  if (length(matched_keys) < BIN_SIZE_LIMIT) {
    stopf(
      "Only %d genera match unambiguous tree tips; fewer than bin.size.limit=%d.",
      length(matched_keys), BIN_SIZE_LIMIT
    )
  }

  if (SMOKE_TEST) {
    cat("[5/9] Selecting deterministic high-coverage smoke samples by realm...\n")
    smoke_metadata <- select_smoke_samples(
      metadata, comm, matched_keys, "meta_group", SMOKE_N_PER_REALM
    )
    comm <- comm[smoke_metadata$run_accession, , drop = FALSE]
    metadata <- smoke_metadata
    cat("    Smoke samples: ", nrow(metadata), "\n", sep = "")
  } else {
    cat("[5/9] Full mode: retaining all realm-assigned samples...\n")
  }

  n_samples_before_tree_coverage_filter <- nrow(comm)
  # Coverage uses all genus-assigned counts before the prevalence filter.
  genus_assigned_count <- rowSums(comm)
  tree_matched_count <- rowSums(comm[, matched_keys, drop = FALSE])
  tree_coverage_fraction <- tree_matched_count / genus_assigned_count
  tree_coverage_fraction[!is.finite(tree_coverage_fraction)] <- NA_real_

  meets_required_coverage <-
    is.finite(tree_coverage_fraction) &
    tree_coverage_fraction >= MIN_SAMPLE_TREE_COVERAGE
  retained_after_tree_coverage_filter <- if (
    FILTER_SAMPLES_BY_TREE_COVERAGE
  ) {
    meets_required_coverage
  } else {
    rep(TRUE, length(meets_required_coverage))
  }
  tree_coverage_filter_reason <- if (
    FILTER_SAMPLES_BY_TREE_COVERAGE
  ) {
    ifelse(
      !is.finite(tree_coverage_fraction),
      "non_finite_tree_coverage",
      ifelse(
        tree_coverage_fraction < MIN_SAMPLE_TREE_COVERAGE,
        "below_tree_coverage_threshold",
        "retained"
      )
    )
  } else {
    rep("filter_not_applied", length(meets_required_coverage))
  }

  coverage_audit <- data.frame(
    run_accession = rownames(comm),
    meta_group = metadata$meta_group,
    genus_assigned_count = genus_assigned_count,
    tree_matched_count = tree_matched_count,
    tree_coverage_fraction = tree_coverage_fraction,
    meets_required_coverage = meets_required_coverage,
    below_required_coverage = !meets_required_coverage,
    tree_coverage_filter_applied = FILTER_SAMPLES_BY_TREE_COVERAGE,
    retained_after_tree_coverage_filter =
      retained_after_tree_coverage_filter,
    tree_coverage_filter_reason = tree_coverage_filter_reason,
    stringsAsFactors = FALSE
  )
  write_tsv(
    coverage_audit,
    file.path(TABLES_DIR, "08B_tree_coverage_by_sample.tsv")
  )

  tree_coverage_exclusions <- coverage_audit[
    !coverage_audit$retained_after_tree_coverage_filter,
    ,
    drop = FALSE
  ]
  write_tsv(
    tree_coverage_exclusions,
    file.path(TABLES_DIR, "08B_tree_coverage_excluded_samples.tsv")
  )

  min_tree_coverage_before_filter <- min(
    tree_coverage_fraction,
    na.rm = TRUE
  )
  median_tree_coverage_before_filter <- stats::median(
    tree_coverage_fraction,
    na.rm = TRUE
  )
  low_coverage_ids <- coverage_audit$run_accession[
    is.na(coverage_audit$tree_coverage_fraction) |
      coverage_audit$below_required_coverage
  ]

  if (
    FILTER_SAMPLES_BY_TREE_COVERAGE &&
    any(!retained_after_tree_coverage_filter)
  ) {
    comm <- comm[
      retained_after_tree_coverage_filter,
      ,
      drop = FALSE
    ]
    metadata <- metadata[
      match(rownames(comm), metadata$run_accession),
      ,
      drop = FALSE
    ]
    if (nrow(comm) < 10L) {
      stop(
        "Fewer than 10 samples remain after tree-coverage filtering.",
        call. = FALSE
      )
    }
    cat(
      "    Samples excluded by tree-coverage filter: ",
      sum(!retained_after_tree_coverage_filter),
      "\n",
      sep = ""
    )
  }

  n_samples_excluded_by_tree_coverage <-
    n_samples_before_tree_coverage_filter - nrow(comm)
  analysis_coverage_index <- match(
    rownames(comm),
    coverage_audit$run_accession
  )
  tree_coverage_fraction <-
    coverage_audit$tree_coverage_fraction[analysis_coverage_index]
  min_tree_coverage <- min(tree_coverage_fraction, na.rm = TRUE)
  median_tree_coverage <- stats::median(
    tree_coverage_fraction,
    na.rm = TRUE
  )

  cat(
    "    Minimum sample tree coverage before filter: ",
    sprintf("%.2f%%", 100 * min_tree_coverage_before_filter),
    "\n",
    sep = ""
  )
  cat(
    "    Median sample tree coverage before filter: ",
    sprintf("%.2f%%", 100 * median_tree_coverage_before_filter),
    "\n",
    sep = ""
  )
  cat(
    "    Minimum sample tree coverage used: ",
    sprintf("%.2f%%", 100 * min_tree_coverage),
    "\n",
    sep = ""
  )
  cat(
    "    Median sample tree coverage used: ",
    sprintf("%.2f%%", 100 * median_tree_coverage),
    "\n",
    sep = ""
  )

  if (
    length(low_coverage_ids) > 0L &&
    !FILTER_SAMPLES_BY_TREE_COVERAGE &&
    !ALLOW_LOW_TREE_COVERAGE
  ) {
    stopf(
      paste(
        "%d sample(s) have tree coverage below %.1f%% or non-finite coverage.",
        "Examples: %s. Supply a more complete tree; use",
        "--allow_low_tree_coverage 1 only for an explicitly documented",
        "sensitivity analysis."
      ),
      length(low_coverage_ids), 100 * MIN_SAMPLE_TREE_COVERAGE,
      paste(utils::head(low_coverage_ids, 10L), collapse = ", ")
    )
  }

  group_counts <- as.data.frame(
    table(metadata$meta_group),
    stringsAsFactors = FALSE
  )
  names(group_counts) <- c("meta_group", "n_samples")
  write_tsv(
    group_counts,
    file.path(TABLES_DIR, "08B_meta_group_sample_counts.tsv")
  )
  if (any(group_counts$n_samples < 2L)) {
    stopf(
      "Metacommunity group(s) with fewer than two samples: %s",
      paste(
        group_counts$meta_group[group_counts$n_samples < 2L],
        collapse = ", "
      )
    )
  }

  current_prevalence <- colSums(comm > 0)
  current_total <- colSums(comm)
  genus_map$prevalence_n_analysis_samples <- current_prevalence[genus_map$genus_key]
  genus_map$total_count_analysis_samples <- current_total[genus_map$genus_key]
  genus_map$passes_prevalence <-
    genus_map$prevalence_n_analysis_samples >= MIN_PREV_N
  genus_map$eligible_for_ranking <-
    genus_map$tree_match & genus_map$passes_prevalence &
    genus_map$total_count_analysis_samples > 0
  genus_map$abundance_rank_among_eligible <- NA_integer_

  eligible_index <- which(genus_map$eligible_for_ranking)
  eligible_index <- eligible_index[order(
    -genus_map$total_count_analysis_samples[eligible_index],
    genus_map$genus_key[eligible_index]
  )]
  if (length(eligible_index) < BIN_SIZE_LIMIT) {
    stopf(
      paste(
        "Only %d genera pass tree matching and prevalence >= %d;",
        "fewer than bin.size.limit=%d."
      ),
      length(eligible_index), MIN_PREV_N, BIN_SIZE_LIMIT
    )
  }
  genus_map$abundance_rank_among_eligible[eligible_index] <-
    seq_along(eligible_index)
  selected_index <- if (N_TOP_GENUS == 0L) {
    eligible_index
  } else {
    utils::head(eligible_index, N_TOP_GENUS)
  }
  genus_map$selected_for_icamp <- seq_len(nrow(genus_map)) %in% selected_index
  write_tsv(genus_map, file.path(TABLES_DIR, "08B_genus_tree_matching.tsv"))

  selected_keys <- genus_map$genus_key[selected_index]
  selected_tips <- genus_map$tree_tip[selected_index]
  comm <- comm[, selected_keys, drop = FALSE]
  colnames(comm) <- selected_tips
  nonzero_taxa <- colSums(comm) > 0
  comm <- comm[, nonzero_taxa, drop = FALSE]
  tree <- ape::keep.tip(tree, colnames(comm))
  comm <- comm[, tree$tip.label, drop = FALSE]
  if (!identical(colnames(comm), tree$tip.label)) {
    stop("Internal error aligning community columns and tree tips.", call. = FALSE)
  }
  if (ncol(comm) < BIN_SIZE_LIMIT) {
    stop("Too few genera remain after final alignment.", call. = FALSE)
  }

  cat("    Tree source: ", tree_result$source, "\n", sep = "")
  cat("    Genera retained for iCAMP: ", ncol(comm), "\n", sep = "")
  cat("    Samples retained for iCAMP: ", nrow(comm), "\n", sep = "")

  meta_group <- data.frame(
    meta.com = metadata$meta_group,
    row.names = metadata$run_accession,
    stringsAsFactors = FALSE
  )
  if (!identical(rownames(comm), rownames(meta_group))) {
    stop("comm and meta.group row names are not identically ordered.",
         call. = FALSE)
  }

  cat("[6/9] Writing prepared inputs and provenance...\n")
  ape::write.tree(tree, file.path(TABLES_DIR, "08B_genus_tree_used.nwk"))
  write_matrix_tsv_gz(
    comm, file.path(TABLES_DIR, "08B_genus_count_matrix_used.tsv.gz")
  )
  write_tsv(
    data.frame(run_accession = rownames(meta_group), meta_group = meta_group$meta.com),
    file.path(TABLES_DIR, "08B_meta_group_used.tsv")
  )
  saveRDS(
    list(comm = comm, tree = tree, meta_group = meta_group,
         metadata = metadata, genus_map = genus_map,
         tree_coverage_audit = coverage_audit),
    file.path(RDS_DIR, "08B_prepared_icamp_inputs.rds")
  )

  parameters <- data.frame(
    parameter = c(
      "script", "script_version", "mode", "analysis_label",
      "phyloseq_input",
      "tree_source", "tree_input",
      "tree_sha256",
      "run_accession_source_column", "ecological_class",
      "excluded_accessions", "meta_group_source_column",
      "n_samples_before_tree_coverage_filter",
      "n_samples_excluded_by_tree_coverage", "n_samples",
      "n_metacommunities", "genus_taxonomy_column", "n_input_taxa",
      "n_features_without_genus", "min_prev_n", "n_top_genus_requested",
      "n_genera_used", "rand", "bin_size_limit", "nworker", "memory_G",
      "seed", "phylo_metric", "sig_index", "phylo_rand_scale",
      "taxa_rand_scale", "taxo_metric", "correct_special", "preflight_only",
      "filter_samples_by_tree_coverage",
      "minimum_sample_tree_coverage_required",
      "minimum_sample_tree_coverage_before_filter",
      "median_sample_tree_coverage_before_filter",
      "minimum_sample_tree_coverage_observed",
      "median_sample_tree_coverage_observed",
      "allow_low_tree_coverage",
      "tree_coverage_excluded_accessions",
      "iCAMP_version", "output_directory"
    ),
    value = c(
      SCRIPT_ID, SCRIPT_VERSION, ifelse(SMOKE_TEST, "SMOKE", "FULL"),
      ifelse(nzchar(ANALYSIS_LABEL), ANALYSIS_LABEL, "none"), ps_path,
      tree_result$source, tree_result$path, tree_sha256, run_col, ECOLOGICAL_CLASS,
      paste(EXCLUDE_RUNS, collapse = ";"), group_col,
      n_samples_before_tree_coverage_filter,
      n_samples_excluded_by_tree_coverage, nrow(comm),
      length(unique(meta_group$meta.com)), genus_result$genus_column,
      genus_result$n_input_features, genus_result$n_features_without_genus,
      MIN_PREV_N, N_TOP_GENUS, ncol(comm), RAND, BIN_SIZE_LIMIT, NWORKER,
      MEMORY_G, SEED, PHYLO_METRIC, SIG_INDEX, "within.bin", "across.all",
      "bray", TRUE, PREFLIGHT_ONLY, FILTER_SAMPLES_BY_TREE_COVERAGE,
      MIN_SAMPLE_TREE_COVERAGE,
      min_tree_coverage_before_filter, median_tree_coverage_before_filter,
      min_tree_coverage, median_tree_coverage, ALLOW_LOW_TREE_COVERAGE,
      paste(tree_coverage_exclusions$run_accession, collapse = ";"),
      as.character(icamp_version), OUT_DIR
    ),
    stringsAsFactors = FALSE
  )
  parameters <- rbind(parameters, data.frame(
    parameter = c("tree_coverage_denominator", "genus_name_matching", "undominated_definition"),
    value = c(
      "All genus-assigned counts before genus-prevalence filtering",
      "Normalized genus names; ambiguous tree-tip keys excluded; no synonym crosswalk",
      "Drift, weak selection or other undominated processes"
    ),
    stringsAsFactors = FALSE
  ))
  write_tsv(parameters, file.path(TABLES_DIR, "08B_parameters.tsv"))
  capture.output(sessionInfo(), file = file.path(LOGS_DIR, "08B_sessionInfo.txt"))

  validation <- data.frame(
    check = c(
      "run_accession_unique", "comm_metadata_identical_order",
      "comm_tree_identical_taxa", "comm_values_finite_nonnegative",
      "all_samples_positive_library", "n_genera_at_least_bin_size",
      "all_meta_groups_at_least_two_samples",
      "sample_tree_coverage_finite",
      "sample_tree_coverage_threshold_or_documented_override",
      "tree_coverage_filter_count_consistent",
      "tree_coverage_filter_retained_only_eligible"
    ),
    pass = c(
      !anyNA(rownames(comm)) && !anyDuplicated(rownames(comm)),
      identical(rownames(comm), rownames(meta_group)),
      identical(colnames(comm), tree$tip.label),
      all(is.finite(comm)) && all(comm >= 0),
      all(rowSums(comm) > 0),
      ncol(comm) >= BIN_SIZE_LIMIT,
      all(group_counts$n_samples >= 2L),
      all(is.finite(tree_coverage_fraction)),
      all(is.finite(tree_coverage_fraction)) &&
        (
          all(tree_coverage_fraction >= MIN_SAMPLE_TREE_COVERAGE) ||
          (!FILTER_SAMPLES_BY_TREE_COVERAGE && ALLOW_LOW_TREE_COVERAGE)
        ),
      nrow(comm) ==
        n_samples_before_tree_coverage_filter -
        n_samples_excluded_by_tree_coverage,
      !FILTER_SAMPLES_BY_TREE_COVERAGE ||
        (
          all(tree_coverage_fraction >= MIN_SAMPLE_TREE_COVERAGE) &&
          n_samples_excluded_by_tree_coverage ==
            nrow(tree_coverage_exclusions)
        )
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(validation, file.path(TABLES_DIR, "08B_preflight_validation.tsv"))
  if (any(!validation$pass)) {
    stopf("Preflight validation failed: %s",
          paste(validation$check[!validation$pass], collapse = ", "))
  }

  preflight_summary <- data.frame(
    item = c(
      "status", "mode", "analysis_label",
      "n_samples_before_tree_coverage_filter",
      "n_samples_excluded_by_tree_coverage", "n_samples",
      "n_genera", "n_metacommunities", "min_prev_n",
      "filter_samples_by_tree_coverage",
      "minimum_sample_tree_coverage_required",
      "minimum_sample_tree_coverage_before_filter",
      "minimum_sample_tree_coverage_observed",
      "median_sample_tree_coverage_observed",
      "tree_sha256", "output_directory"
    ),
    value = c(
      "PASS_PREFLIGHT_INPUTS", ifelse(SMOKE_TEST, "SMOKE", "FULL"),
      ifelse(nzchar(ANALYSIS_LABEL), ANALYSIS_LABEL, "none"),
      n_samples_before_tree_coverage_filter,
      n_samples_excluded_by_tree_coverage, nrow(comm), ncol(comm),
      length(unique(meta_group$meta.com)), MIN_PREV_N,
      FILTER_SAMPLES_BY_TREE_COVERAGE, MIN_SAMPLE_TREE_COVERAGE,
      min_tree_coverage_before_filter, min_tree_coverage,
      median_tree_coverage, tree_sha256, OUT_DIR
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(
    preflight_summary,
    file.path(TABLES_DIR, "08B_preflight_summary.tsv")
  )

  if (PREFLIGHT_ONLY) {
    cat("[7/9] PREFLIGHT_ONLY=1: iCAMP execution intentionally skipped.\n")
    write_status(OUT_DIR, "PASS_PREFLIGHT")
    return(invisible(list(preflight_only = TRUE)))
  }

  cat("[7/9] Running iCAMP::icamp.cm...\n")
  icamp_result <- iCAMP::icamp.cm(
    comm = comm,
    tree = tree,
    meta.group = meta_group,
    pd.desc = NULL,
    pd.spname = NULL,
    pd.wd = PD_WD,
    rand = RAND,
    prefix = ifelse(SMOKE_TEST, "iCAMP_genus_smoke", "iCAMP_genus"),
    phylo.rand.scale = "within.bin",
    taxa.rand.scale = "across.all",
    phylo.metric = PHYLO_METRIC,
    sig.index = SIG_INDEX,
    bin.size.limit = BIN_SIZE_LIMIT,
    nworker = NWORKER,
    memory.G = MEMORY_G,
    rtree.save = TRUE,
    detail.save = DETAIL_SAVE,
    qp.save = TRUE,
    detail.null = FALSE,
    ignore.zero = TRUE,
    output.wd = OUT_DIR,
    correct.special = TRUE,
    unit.sum = rowSums(comm),
    special.method = "depend",
    ses.cut = 1.96,
    rc.cut = 0.95,
    conf.cut = 0.975,
    omit.option = "no",
    taxo.metric = "bray",
    transform.method = NULL,
    dirichlet = FALSE
  )

  saveRDS(icamp_result, file.path(RDS_DIR, "08B_icamp_cm_result_genus.rds"))
  result_names <- names(icamp_result)
  if (is.null(result_names) || any(result_names == "")) {
    result_names <- paste0("element_", seq_along(icamp_result))
  }
  structure_table <- data.frame(
    element = result_names,
    class = vapply(icamp_result, function(x) paste(class(x), collapse = ";"),
                   character(1)),
    dimensions = vapply(icamp_result, function(x) {
      d <- dim(x)
      if (is.null(d)) as.character(length(x)) else paste(d, collapse = "x")
    }, character(1)),
    stringsAsFactors = FALSE
  )
  write_tsv(structure_table, file.path(TABLES_DIR, "08B_icamp_result_structure.tsv"))

  cat("[8/9] Extracting pairwise process contributions...\n")
  pairwise <- extract_pairwise_processes(icamp_result)
  n_pairwise_turnovers <- NA_integer_
  n_within_realm_turnovers <- NA_integer_
  n_between_realm_turnovers <- NA_integer_
  if (is.null(pairwise)) {
    message <- paste(
      "Could not identify a pairwise five-process table in the iCAMP output.",
      "Inspect 08B_icamp_result_structure.tsv and the result RDS."
    )
    if (STRICT_OUTPUT) stop(message, call. = FALSE) else warnf("%s", message)
  } else {
    source_element <- attr(pairwise, "source_element")
    pairwise <- annotate_pairwise_geography(pairwise, meta_group)
    pair_audit <- audit_pairwise_structure(pairwise, meta_group)

    write_tsv(
      pair_audit$validation,
      file.path(TABLES_DIR, "08B_icamp_pair_structure_validation.tsv")
    )
    write_tsv(
      pair_audit$pair_counts,
      file.path(TABLES_DIR, "08B_icamp_realm_pair_counts.tsv")
    )
    if (any(!pair_audit$validation$pass)) {
      stopf(
        "Pairwise iCAMP structure validation failed: %s",
        paste(
          pair_audit$validation$check[!pair_audit$validation$pass],
          collapse = ", "
        )
      )
    }

    write_tsv_gz(
      pairwise,
      file.path(TABLES_DIR, "08B_icamp_pairwise_processes.tsv.gz")
    )

    process_summary_all <- summarise_processes(pairwise)
    process_summary_all <- cbind(
      comparison_scope = "all_pairs",
      process_summary_all,
      source_element = source_element,
      stringsAsFactors = FALSE
    )
    write_tsv(
      process_summary_all,
      file.path(TABLES_DIR, "08B_icamp_process_summary_all_pairs.tsv")
    )

    process_summary_scope <- summarise_processes(
      pairwise, group_cols = "comparison_scope"
    )
    process_summary_scope$source_element <- source_element
    write_tsv(
      process_summary_scope,
      file.path(TABLES_DIR, "08B_icamp_process_summary.tsv")
    )

    within_realm <- pairwise[
      pairwise$comparison_scope == "within_realm", , drop = FALSE
    ]
    within_realm$realm <- within_realm$realm1
    process_summary_within <- summarise_processes(
      within_realm, group_cols = "realm"
    )
    process_summary_within$source_element <- source_element
    write_tsv(
      process_summary_within,
      file.path(TABLES_DIR, "08B_icamp_process_summary_within_realm.tsv")
    )

    between_realms <- pairwise[
      pairwise$comparison_scope == "between_realms", , drop = FALSE
    ]
    process_summary_between <- summarise_processes(
      between_realms, group_cols = "realm_pair"
    )
    process_summary_between$source_element <- source_element
    write_tsv(
      process_summary_between,
      file.path(
        TABLES_DIR,
        "08B_icamp_process_summary_between_realm_pair.tsv"
      )
    )

    dominant_scope <- summarise_dominant_process(
      pairwise, group_cols = "comparison_scope"
    )
    write_tsv(
      dominant_scope,
      file.path(TABLES_DIR, "08B_dominant_process_counts.tsv")
    )

    dominant_within <- summarise_dominant_process(
      within_realm, group_cols = "realm"
    )
    write_tsv(
      dominant_within,
      file.path(
        TABLES_DIR,
        "08B_dominant_process_counts_within_realm.tsv"
      )
    )

    n_pairwise_turnovers <- nrow(pairwise)
    n_within_realm_turnovers <- pair_audit$observed_within
    n_between_realm_turnovers <- pair_audit$observed_between
    cat("    Pairwise turnovers: ", n_pairwise_turnovers, "\n", sep = "")
    cat("    Within-realm turnovers: ", n_within_realm_turnovers,
        " (primary interpretation scope)\n", sep = "")
    cat("    Between-realm turnovers: ", n_between_realm_turnovers,
        "\n", sep = "")
  }

  cat("[9/9] Finalizing run...\n")
  run_summary <- data.frame(
    item = c(
      "status", "mode", "analysis_label",
      "n_samples_before_tree_coverage_filter",
      "n_samples_excluded_by_tree_coverage", "n_samples",
      "n_genera", "n_metacommunities",
      "rand", "min_prev_n", "bin_size_limit", "phylo_metric", "sig_index",
      "filter_samples_by_tree_coverage",
      "minimum_sample_tree_coverage_required",
      "iCAMP_version", "n_pairwise_turnovers",
      "n_within_realm_turnovers", "n_between_realm_turnovers",
      "primary_interpretation_scope", "output_directory"
    ),
    value = c(
      "PASS", ifelse(SMOKE_TEST, "SMOKE", "FULL"),
      ifelse(nzchar(ANALYSIS_LABEL), ANALYSIS_LABEL, "none"),
      n_samples_before_tree_coverage_filter,
      n_samples_excluded_by_tree_coverage, nrow(comm), ncol(comm),
      length(unique(meta_group$meta.com)), RAND, MIN_PREV_N, BIN_SIZE_LIMIT,
      PHYLO_METRIC, SIG_INDEX, FILTER_SAMPLES_BY_TREE_COVERAGE,
      MIN_SAMPLE_TREE_COVERAGE, as.character(icamp_version),
      n_pairwise_turnovers, n_within_realm_turnovers,
      n_between_realm_turnovers, "within_realm", OUT_DIR
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(run_summary, file.path(TABLES_DIR, "08B_summary.tsv"))
  write_status(OUT_DIR, "PASS")
  cat("\nPASS: ", OUT_DIR, "\n", sep = "")
  invisible(icamp_result)
}

run_error <- NULL
tryCatch(
  {
    main()
    status_ok <- TRUE
  },
  error = function(e) {
    message("\n[FATAL] ", conditionMessage(e))
    write_status(OUT_DIR, "FAILED", conditionMessage(e))
    run_error <<- e
  },
  finally = {
    if (!status_ok && is.null(run_error)) {
      try(write_status(OUT_DIR, "FAILED", "See logs/08B_run.log"), silent = TRUE)
    }
    try(sink(type = "message"), silent = TRUE)
    try(sink(type = "output"), silent = TRUE)
    try(close(log_connection), silent = TRUE)
  }
)

if (!is.null(run_error)) stop(run_error)
