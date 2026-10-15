#!/usr/bin/env Rscript

# =============================================================================
# 144_BUILD_KO_FUNCTIONAL_ANNOTATION_SNAPSHOT.R
#
# Purpose
# -------
# Create the frozen KO functional-annotation snapshot consumed by the final
# network-annotation script 173.
#
# Reproducibility corrections relative to the historical script
# ---------------------------------------------------------------
# 1. This script annotates a KO universe only. It no longer rebuilds or plots an
#    obsolete 140/143 network; final network construction, curation, annotation
#    and plotting belong to scripts 171-175.
# 2. KEGG MODULE mappings must come from a versioned local snapshot. A live
#    KEGGREST query is possible only with --allow_kegg_query 1, and the returned
#    mapping is immediately frozen in the output directory.
# 3. Functional categories use the historical rule set but now record category
#    scores and ties, making heuristic assignments auditable.
# 4. Canonical and legacy filenames are written together so script 173 remains
#    compatible without modification.
#
# Required inputs
# ---------------
# KO universe: supplied with --ko_table or discovered in the latest script-156
#              run (156_KO_network_filter_stats.tsv).
# KO functions: supplied with --ko_dictionary or discovered in script 137.
# KEGG mapping: supplied with --kegg_snapshot, unless a deliberate live query
#               is authorised with --allow_kegg_query 1.
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

clean_chr <- function(x) {
  y <- trimws(as.character(x))
  y[y %in% c("", "NA", "NaN", "NULL", "None", ".")] <- NA_character_
  y
}

read_latest <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  x <- clean_chr(readLines(path, warn = FALSE))
  x <- x[!is.na(x)]
  if (length(x) == 0L) NA_character_ else x[1L]
}

first_existing_latest <- function(paths) {
  for (path in paths) {
    value <- read_latest(path)
    if (!is.na(value) && dir.exists(value)) return(value)
  }
  NA_character_
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

find_column <- function(x, candidates, label, required = TRUE) {
  exact <- candidates[candidates %in% names(x)]
  if (length(exact) > 0L) return(exact[1L])
  idx <- match(tolower(candidates), tolower(names(x)))
  idx <- idx[!is.na(idx)]
  if (length(idx) > 0L) return(names(x)[idx[1L]])
  if (required) {
    stopf("Could not identify %s. Tried: %s",
          label, paste(candidates, collapse = ", "))
  }
  NA_character_
}

normalise_KO <- function(x) {
  y <- toupper(clean_chr(x))
  y <- sub("^KO:", "", y)
  y <- sub("^ORTHOLOGY:", "", y)
  y
}

normalise_module <- function(x) {
  y <- toupper(clean_chr(x))
  y <- sub("^MD:", "", y)
  y <- sub("^MODULE:", "", y)
  y
}

collapse_unique <- function(x) {
  x <- sort(unique(clean_chr(x)))
  x <- x[!is.na(x)]
  if (length(x) == 0L) NA_character_ else paste(x, collapse = "; ")
}

collapse_annotation_by_KO <- function(x, value_col, output_col) {
  x <- x[!is.na(x$KO), c("KO", value_col), drop = FALSE]
  if (nrow(x) == 0L) {
    out <- data.frame(KO = character(), value = character(),
                      n_distinct_values = integer(), stringsAsFactors = FALSE)
  } else {
    split_values <- split(x[[value_col]], x$KO)
    out <- data.frame(
      KO = names(split_values),
      value = vapply(split_values, collapse_unique, character(1)),
      n_distinct_values = vapply(split_values, function(z) {
        length(unique(stats::na.omit(clean_chr(z))))
      }, integer(1)),
      stringsAsFactors = FALSE
    )
  }
  names(out)[2L] <- output_col
  out
}

standardise_KO_universe <- function(path) {
  x <- read_tsv_auto(path)
  ko_col <- find_column(x, c(
    "KO", "ko", "ko_id", "kegg_ko", "KEGG_Ortholog", "orthology",
    "feature_id", "Id", "id"
  ),
                        "KO-universe identifier column")
  ko <- normalise_KO(x[[ko_col]])
  ko <- ko[!is.na(ko)]
  invalid <- !grepl("^K[0-9]{5}$", ko)
  if (any(invalid)) {
    stopf("KO universe contains invalid identifiers (first values: %s)",
          paste(utils::head(unique(ko[invalid]), 10L), collapse = ", "))
  }
  data.frame(KO = unique(ko), source_order = seq_along(unique(ko)),
             stringsAsFactors = FALSE)
}

standardise_KO_dictionary <- function(path) {
  x <- read_tsv_auto(path)
  ko_col <- find_column(x, c(
    "KO", "ko", "ko_id", "kegg_ko", "KEGG_Ortholog", "orthology",
    "feature_id", "Id", "id"
  ),
                        "KO-dictionary identifier column")
  function_col <- find_column(
    x,
    c("KO_function", "ko_function", "function", "Function", "description",
      "KO_description", "KO_definition", "gene_name", "definition", "name"),
    "KO function/description column"
  )
  z <- data.frame(KO = normalise_KO(x[[ko_col]]),
                  KO_function = clean_chr(x[[function_col]]),
                  stringsAsFactors = FALSE)
  collapse_annotation_by_KO(z, "KO_function", "KO_function")
}

standardise_core_classes <- function(path) {
  x <- read_tsv_auto(path)
  ko_col <- find_column(x, c(
    "KO", "ko", "ko_id", "kegg_ko", "KEGG_Ortholog", "orthology"
  ),
                        "core-class KO column")
  class_col <- find_column(
    x,
    c("functional_core_class", "recommended_functional_core_class",
      "functional_core_class_minloc3", "reconciled_class", "class_primary"),
    "functional core-class column"
  )
  z <- data.frame(KO = normalise_KO(x[[ko_col]]),
                  functional_core_class = clean_chr(x[[class_col]]),
                  stringsAsFactors = FALSE)
  collapse_annotation_by_KO(z, "functional_core_class",
                            "functional_core_class")
}

standardise_KEGG_snapshot <- function(path) {
  x <- read_tsv_auto(path)
  ko_col <- find_column(x, c("KO", "ko", "KEGG_KO", "orthology"),
                        "KEGG-snapshot KO column")
  module_col <- find_column(
    x,
    c("KEGG_module_id", "module_id", "Module", "module", "KEGG_MODULE"),
    "KEGG MODULE identifier column"
  )
  name_col <- find_column(
    x,
    c("KEGG_module_name", "module_name", "Module_name", "description", "name"),
    "KEGG MODULE name column",
    required = FALSE
  )
  out <- data.frame(
    KO = normalise_KO(x[[ko_col]]),
    KEGG_module_id = normalise_module(x[[module_col]]),
    KEGG_module_name = if (is.na(name_col)) NA_character_ else clean_chr(x[[name_col]]),
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$KO) & !is.na(out$KEGG_module_id), , drop = FALSE]
  invalid_KO <- !grepl("^K[0-9]{5}$", out$KO)
  invalid_module <- !grepl("^M[0-9]{5}$", out$KEGG_module_id)
  if (any(invalid_KO) || any(invalid_module)) {
    stopf("KEGG snapshot contains %d invalid KO and %d invalid MODULE identifier(s).",
          sum(invalid_KO), sum(invalid_module))
  }
  unique(out)
}

query_KEGG_modules <- function() {
  if (!requireNamespace("KEGGREST", quietly = TRUE)) {
    stop("A live query was requested, but package KEGGREST is not installed.",
         call. = FALSE)
  }
  links <- KEGGREST::keggLink("module", "ko")
  modules <- KEGGREST::keggList("module")
  out <- data.frame(
    KO = normalise_KO(names(links)),
    KEGG_module_id = normalise_module(unname(links)),
    stringsAsFactors = FALSE
  )
  module_lookup <- data.frame(
    KEGG_module_id = normalise_module(names(modules)),
    KEGG_module_name = clean_chr(unname(modules)),
    stringsAsFactors = FALSE
  )
  out$KEGG_module_name <- module_lookup$KEGG_module_name[
    match(out$KEGG_module_id, module_lookup$KEGG_module_id)
  ]
  unique(out[!is.na(out$KO) & !is.na(out$KEGG_module_id), , drop = FALSE])
}

FUNCTIONAL_CATEGORY_ORDER <- c(
  "Anaerobic redox / energy",
  "Transport / membrane systems",
  "Information processing",
  "Amino acid / cofactor biosynthesis",
  "Carbon / carbohydrate metabolism",
  "Signal / stress / metals",
  "Uncharacterized / other"
)

FUNCTIONAL_PATTERNS <- list(
  "Anaerobic redox / energy" = c(
    "ferredoxin", "oxidoreductase", "reductase", "dehydrogenase",
    "hydrogenase", "heterodisulfide", "formate", "carbon-monoxide",
    "carbon monoxide", "co dehydrogenase", "codh", "acetyl-coa synthase",
    "qmo", "hdr", "dms", "dimethyl sulfoxide", "nitrate", "nitrite",
    "sulfate", "sulfite", "sulfur", "methan", "nadh-quinone",
    "succinate dehydrogenase", "fumarate reductase", "electron transfer"
  ),
  "Transport / membrane systems" = c(
    "transport", "abc", "permease", "substrate-binding", "binding protein",
    "efflux", "membrane fusion", "transporter", "translocation",
    "nickel transport", "peptide transport", "sugar transport"
  ),
  "Information processing" = c(
    "trna", "synthetase", "ribosomal", "ribosome", "translation",
    "transcription", "polymerase", "elongation factor", "nusa",
    "uvr", "mutl", "muts", "repair", "replication", "recombin",
    "preprotein translocase", "seca"
  ),
  "Amino acid / cofactor biosynthesis" = c(
    "leucine", "isoleucine", "valine", "aspartate", "homoserine",
    "threonine", "lysine", "dap", "amino acid", "aminotransferase",
    "molybdopterin", "nad", "cofactor", "pantothenate", "biotin",
    "folate", "methylenetetrahydrofolate", "carbamoyl"
  ),
  "Carbon / carbohydrate metabolism" = c(
    "glycolysis", "gluconeogenesis", "sugar", "ribose", "galactose",
    "carbohydrate", "acetate", "pyruvate", "lactate", "acetoin",
    "coa ligase", "phosphofructokinase", "glucan", "starch",
    "cellulose", "glycosidase"
  ),
  "Signal / stress / metals" = c(
    "two-component", "response regulator", "sensor kinase", "copper",
    "silver", "metal", "selenocysteine", "stress", "chaperone",
    "heat shock", "cold shock", "detox", "resistance"
  )
)

score_functional_category <- function(function_text, module_text) {
  text <- tolower(paste(
    ifelse(is.na(function_text), "", function_text),
    ifelse(is.na(module_text), "", module_text)
  ))
  scores <- vapply(FUNCTIONAL_PATTERNS, function(patterns) {
    sum(vapply(patterns, function(pattern) {
      grepl(pattern, text, fixed = TRUE)
    }, logical(1)))
  }, numeric(1))
  max_score <- max(scores)
  if (max_score == 0L) {
    selected <- "Uncharacterized / other"
    tied <- FALSE
    candidates <- NA_character_
  } else {
    top <- names(scores)[scores == max_score]
    # FUNCTIONAL_PATTERNS order is the documented deterministic tie breaker.
    selected <- top[1L]
    tied <- length(top) > 1L
    candidates <- paste(top, collapse = "; ")
  }
  list(selected = selected, max_score = max_score, tied = tied,
       candidates = candidates, scores = scores)
}

print_help <- function() {
  cat(paste(
    "Usage:",
    "  Rscript 144_BUILD_KO_FUNCTIONAL_ANNOTATION_SNAPSHOT.R \\",
    "    --kegg_snapshot PATH [--ko_table PATH] [--ko_dictionary PATH]",
    "",
    "Key options:",
    "  --base PATH                  Project base directory",
    "  --results_dir PATH           Results directory",
    "  --run156 PATH                Corrected/legacy script-156 run",
    "  --run137 PATH                Script-137 KO dictionary run",
    "  --run153 PATH                Optional final KO core-class run",
    "  --allow_kegg_query 0|1       Default 0; explicit opt-in only",
    "  --strict_expected 0|1        Expect 5,000 network-universe KOs",
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
  run156 <- get_arg(
    "--run156",
    first_existing_latest(c(
      file.path(results_dir,
                "LATEST_156_PREPARE_KO_CLR_MATRICES_FOR_ASSOCIATION_AND_NETWORK.txt"),
      file.path(results_dir,
                "LATEST_156_build_KO_CLR_matrix_prevalence_classes.txt")
    ))
  )
  run137 <- get_arg(
    "--run137",
    read_latest(file.path(
      results_dir, "LATEST_137_plot_global_PCA_and_PC1_core_structure.txt"
    ))
  )
  run153 <- get_arg(
    "--run153",
    read_latest(file.path(
      results_dir, "LATEST_153_FIG4A_FUNCTIONAL_KO_CLASSES_RECONCILED.txt"
    ))
  )

  ko_table <- get_arg("--ko_table", NA_character_)
  if ((is.na(ko_table) || ko_table == "") && !is.na(run156)) {
    ko_table <- file.path(run156, "tables", "156_KO_network_filter_stats.tsv")
  }
  ko_dictionary <- get_arg("--ko_dictionary", NA_character_)
  if ((is.na(ko_dictionary) || ko_dictionary == "") && !is.na(run137)) {
    ko_dictionary <- file.path(
      run137, "tables", "137_KO_function_dictionary_used.tsv"
    )
  }
  core_classes <- get_arg("--core_classes", NA_character_)
  if ((is.na(core_classes) || core_classes == "") && !is.na(run153)) {
    candidate <- file.path(
      run153, "tables", "153_KO_classification_primary_minloc3.tsv"
    )
    if (file.exists(candidate)) core_classes <- candidate
  }
  kegg_snapshot <- get_arg("--kegg_snapshot", NA_character_)
  allow_kegg_query <- as_flag(get_arg("--allow_kegg_query", "0"))
  strict_expected <- as_flag(get_arg("--strict_expected", "1"))
  expected_n_KOs <- as.integer(get_arg("--expected_n_kos", "5000"))

  if (is.na(ko_table) || ko_table == "" || !file.exists(ko_table)) {
    stop("KO universe not found. Supply --ko_table or a valid --run156.",
         call. = FALSE)
  }
  if (is.na(ko_dictionary) || ko_dictionary == "" ||
      !file.exists(ko_dictionary)) {
    stop("KO function dictionary not found. Supply --ko_dictionary or a valid --run137.",
         call. = FALSE)
  }
  if ((is.na(kegg_snapshot) || kegg_snapshot == "") && !allow_kegg_query) {
    stop(paste(
      "A frozen KEGG MODULE snapshot is required for reproducibility.",
      "Supply --kegg_snapshot PATH, or deliberately create one with",
      "--allow_kegg_query 1."
    ), call. = FALSE)
  }

  script_name <- "144_BUILD_KO_FUNCTIONAL_ANNOTATION_SNAPSHOT"
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  outdir <- file.path(results_dir, paste0(script_name, "_", stamp))
  tables_dir <- file.path(outdir, "tables")
  logs_dir <- file.path(outdir, "logs")
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
  latest_files <- c(
    file.path(results_dir, paste0("LATEST_", script_name, ".txt")),
    file.path(results_dir,
              "LATEST_144_annotate_functional_network_modules.txt")
  )
  invisible(lapply(latest_files, function(path) writeLines(outdir, path)))
  writeLines("RUNNING", file.path(outdir, "STATUS.txt"))
  writeLines(c(
    "status\tRUNNING",
    paste0("started_utc\t", format(Sys.time(), tz = "UTC", usetz = TRUE))
  ), file.path(logs_dir, "144_status.tsv"))

  execute <- function() {
    cat("[1/6] Reading the KO universe and function dictionary...\n")
    universe <- standardise_KO_universe(ko_table)
    dictionary <- standardise_KO_dictionary(ko_dictionary)
    universe$KO_function <- dictionary$KO_function[
      match(universe$KO, dictionary$KO)
    ]
    universe$n_distinct_function_descriptions <- dictionary$n_distinct_values[
      match(universe$KO, dictionary$KO)
    ]
    universe$KO_function[is.na(universe$KO_function)] <- universe$KO[
      is.na(universe$KO_function)
    ]
    universe$n_distinct_function_descriptions[
      is.na(universe$n_distinct_function_descriptions)
    ] <- 0L

    cat("[2/6] Adding the optional final functional core classes...\n")
    universe$functional_core_class <- NA_character_
    class_conflicts <- data.frame(
      KO = character(), functional_core_class = character(),
      n_distinct_values = integer(), stringsAsFactors = FALSE
    )
    if (!is.na(core_classes) && core_classes != "" && file.exists(core_classes)) {
      classes <- standardise_core_classes(core_classes)
      universe$functional_core_class <- classes$functional_core_class[
        match(universe$KO, classes$KO)
      ]
      class_conflicts <- classes[classes$n_distinct_values > 1L, , drop = FALSE]
      if (nrow(class_conflicts) > 0L) {
        write_tsv(class_conflicts,
                  file.path(tables_dir,
                            "144_conflicting_functional_core_classes.tsv"))
        stopf("Functional core-class table contains %d conflicting KO assignment(s).",
              nrow(class_conflicts))
      }
    }

    cat("[3/6] Reading or deliberately querying the KEGG MODULE snapshot...\n")
    if (!is.na(kegg_snapshot) && kegg_snapshot != "") {
      kegg_map <- standardise_KEGG_snapshot(kegg_snapshot)
      kegg_source <- "local_frozen_snapshot"
      kegg_source_detail <- normalizePath(kegg_snapshot, winslash = "/",
                                          mustWork = TRUE)
      kegg_query_utc <- NA_character_
    } else {
      kegg_map <- query_KEGG_modules()
      kegg_source <- "live_KEGGREST_query_then_frozen"
      kegg_source_detail <- paste0("KEGGREST_", as.character(
        utils::packageVersion("KEGGREST")
      ))
      kegg_query_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
    }
    kegg_map <- kegg_map[kegg_map$KO %in% universe$KO, , drop = FALSE]
    kegg_map <- unique(kegg_map[order(kegg_map$KO, kegg_map$KEGG_module_id),,
                                drop = FALSE])
    if (nrow(kegg_map) == 0L) {
      stop("No KEGG MODULE links match the requested KO universe.", call. = FALSE)
    }

    # Canonical snapshot and the legacy alias consumed by script 173.
    write_tsv(kegg_map,
              file.path(tables_dir,
                        "144_KO_to_KEGG_MODULE_snapshot_long.tsv"))
    write_tsv(kegg_map,
              file.path(tables_dir, "144_KO_to_KEGG_MODULE_long.tsv"))

    cat("[4/6] Collapsing KEGG MODULE annotations by KO...\n")
    split_map <- split(kegg_map, kegg_map$KO)
    collapsed <- data.frame(
      KO = names(split_map),
      n_KEGG_modules = vapply(split_map, function(z) {
        length(unique(z$KEGG_module_id))
      }, integer(1)),
      KEGG_module_ids = vapply(split_map, function(z) {
        paste(sort(unique(z$KEGG_module_id)), collapse = "; ")
      }, character(1)),
      KEGG_module_names = vapply(split_map, function(z) {
        collapse_unique(z$KEGG_module_name)
      }, character(1)),
      stringsAsFactors = FALSE
    )
    universe$n_KEGG_modules <- collapsed$n_KEGG_modules[
      match(universe$KO, collapsed$KO)
    ]
    universe$KEGG_module_ids <- collapsed$KEGG_module_ids[
      match(universe$KO, collapsed$KO)
    ]
    universe$KEGG_module_names <- collapsed$KEGG_module_names[
      match(universe$KO, collapsed$KO)
    ]
    universe$n_KEGG_modules[is.na(universe$n_KEGG_modules)] <- 0L

    cat("[5/6] Assigning auditable functional categories...\n")
    scored <- Map(
      score_functional_category,
      universe$KO_function,
      universe$KEGG_module_names
    )
    universe$functional_category <- vapply(
      scored, function(z) z$selected, character(1)
    )
    universe$functional_category_max_score <- vapply(
      scored, function(z) z$max_score, numeric(1)
    )
    universe$functional_category_tie <- vapply(
      scored, function(z) z$tied, logical(1)
    )
    universe$functional_category_tie_candidates <- vapply(
      scored, function(z) z$candidates, character(1)
    )
    for (category in names(FUNCTIONAL_PATTERNS)) {
      safe_name <- gsub("[^A-Za-z0-9]+", "_", category)
      universe[[paste0("score_", safe_name)]] <- vapply(
        scored, function(z) unname(z$scores[category]), numeric(1)
      )
    }
    universe <- universe[order(universe$source_order), , drop = FALSE]
    universe$source_order <- NULL

    annotation_output <- universe[c(
      "KO", "KO_function", "functional_core_class", "functional_category",
      "n_KEGG_modules", "KEGG_module_ids", "KEGG_module_names",
      "n_distinct_function_descriptions"
    )]
    write_tsv(annotation_output,
              file.path(tables_dir,
                        "144_KO_functional_annotation_snapshot.tsv"))
    write_tsv(annotation_output,
              file.path(tables_dir, "144_network_nodes_annotated.tsv"))
    write_tsv(universe,
              file.path(tables_dir,
                        "144_functional_category_assignment_audit.tsv"))

    cat("[6/6] Validating and saving provenance...\n")
    category_counts <- as.data.frame(table(
      factor(annotation_output$functional_category,
             levels = FUNCTIONAL_CATEGORY_ORDER)
    ), stringsAsFactors = FALSE)
    names(category_counts) <- c("functional_category", "n_KOs")
    write_tsv(category_counts,
              file.path(tables_dir, "144_functional_category_counts.tsv"))

    validation <- data.frame(
      metric = c(
        "n_KOs", "n_duplicate_KOs", "n_missing_functional_categories",
        "n_KOs_without_dictionary_description", "n_KOs_with_KEGG_MODULE",
        "n_KEGG_MODULE_links", "n_functional_category_ties"
      ),
      observed = c(
        nrow(annotation_output), anyDuplicated(annotation_output$KO),
        sum(is.na(annotation_output$functional_category)),
        sum(annotation_output$KO_function == annotation_output$KO),
        sum(annotation_output$n_KEGG_modules > 0L), nrow(kegg_map),
        sum(universe$functional_category_tie)
      ),
      stringsAsFactors = FALSE
    )
    write_tsv(validation, file.path(tables_dir, "144_validation.tsv"))
    if (strict_expected && nrow(annotation_output) != expected_n_KOs) {
      stopf("Expected %d network-universe KOs but observed %d.",
            expected_n_KOs, nrow(annotation_output))
    }
    if (anyDuplicated(annotation_output$KO) ||
        anyNA(annotation_output$functional_category)) {
      stop("KO uniqueness or functional-category completeness validation failed.",
           call. = FALSE)
    }

    provenance <- data.frame(
      field = c(
        "script", "ko_universe_input", "ko_dictionary_input",
        "functional_core_classes_input", "kegg_source", "kegg_source_detail",
        "kegg_query_utc", "allow_kegg_query", "strict_expected",
        "expected_n_KOs", "functional_category_tie_breaker"
      ),
      value = c(
        script_name, normalizePath(ko_table, winslash = "/", mustWork = TRUE),
        normalizePath(ko_dictionary, winslash = "/", mustWork = TRUE),
        ifelse(is.na(core_classes), NA_character_, core_classes),
        kegg_source, kegg_source_detail, kegg_query_utc, allow_kegg_query,
        strict_expected, expected_n_KOs,
        "first category in documented FUNCTIONAL_PATTERNS order"
      ),
      stringsAsFactors = FALSE
    )
    write_tsv(provenance, file.path(tables_dir, "144_provenance.tsv"))
    write_tsv(data.frame(
      metric = c("n_KOs", "n_KEGG_MODULE_links", "n_KOs_with_modules",
                 "n_category_ties"),
      value = c(nrow(annotation_output), nrow(kegg_map),
                sum(annotation_output$n_KEGG_modules > 0L),
                sum(universe$functional_category_tie)),
      stringsAsFactors = FALSE
    ), file.path(tables_dir, "144_summary.tsv"))
    writeLines(c(
      "status\tPASS",
      paste0("completed_utc\t", format(Sys.time(), tz = "UTC", usetz = TRUE)),
      paste0("output_directory\t", outdir)
    ), file.path(logs_dir, "144_status.tsv"))
    writeLines("PASS", file.path(outdir, "STATUS.txt"))
    cat("PASS: ", outdir, "\n", sep = "")
  }

  tryCatch(
    execute(),
    error = function(e) {
      writeLines("FAIL", file.path(outdir, "STATUS.txt"))
      writeLines(c(
        "status\tFAIL",
        paste0("failed_utc\t", format(Sys.time(), tz = "UTC", usetz = TRUE)),
        paste0("message\t", conditionMessage(e)),
        paste0("output_directory\t", outdir)
      ), file.path(logs_dir, "144_status.tsv"))
      stop(e)
    }
  )
  invisible(outdir)
}

main()
