#!/usr/bin/env Rscript

options(warn = 1)

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Se requiere el paquete 'data.table'.")
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Se requiere el paquete 'igraph'.")
  }
})

# =============================================================================
# 171_FUNCTIONAL_KO_CONSENSUS_NETWORK_FINAL.R
#
# Objetivo
#   Inferir una red consenso de coabundancia funcional a partir de perfiles CLR
#   de KEGG orthologues (KOs), controlando el desequilibrio geografico mediante
#   submuestreo estratificado. Los nodos se anotan con la definicion funcional
#   final del manuscrito, calculada por el script 153:
#
#     absolute_global_core: detectado en 100% de 194 metagenomas.
#     near_global: prevalencia >=95% y <100%.
#     realm_core_only: presente en todas las localidades de al menos un realm
#       con >=3 localidades y abundancia relativa media realm-locality >=0.0001.
#     multi_realm_non_core: detectado en mas de un realm sin cumplir lo anterior.
#     realm_specific: detectado en un solo realm sin cumplir lo anterior.
#
# Configuracion primaria del manuscrito
#   - 194 metagenomas funcionales.
#   - 1,500 KOs candidatos, seleccionados por varianza CLR.
#   - KOs candidatos de red: prevalencia >=20 y abundancia total >=50,
#     establecidos por el script 156.
#   - 100 submuestreos balanceados por marine realm, country y locality.
#   - Maximo 4 metagenomas por estrato.
#   - 2,000 asociaciones de proporcionalidad mas fuertes por submuestreo.
#   - Arista consenso: estabilidad >=0.80, consistencia de signo =1.00 y
#     asociacion absoluta media >=0.95.
#
# Interpretacion
#   Las aristas representan asociaciones reproducibles de coabundancia en
#   espacio CLR. NO demuestran interacciones metabolicas, intercambio de
#   metabolitos, causalidad ni actividad funcional.
# =============================================================================

get_arg <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- which(args == flag)
  if (length(hit) == 0L) return(default)
  hit <- hit[length(hit)]
  if (hit >= length(args)) return(default)
  args[hit + 1L]
}

as_flag <- function(x) {
  tolower(trimws(as.character(x))) %in% c("1", "true", "t", "yes", "y")
}

BASE <- get_arg("--base", path.expand("~/Biogeografia_2026"))
RESULTS_DIR <- get_arg("--results_dir", file.path(BASE, "resultados"))

read_latest <- function(script_basename) {
  path <- file.path(RESULTS_DIR, paste0("LATEST_", script_basename, ".txt"))
  if (!file.exists(path)) return(NA_character_)
  value <- readLines(path, warn = FALSE)
  if (length(value) == 0L) return(NA_character_)
  trimws(value[1L])
}

RUN156 <- get_arg(
  "--run156",
  read_latest("156_build_KO_CLR_matrix_prevalence_classes")
)
RUN153 <- get_arg(
  "--run153",
  read_latest("153_FIG4A_FUNCTIONAL_KO_CLASSES_RECONCILED")
)
RUN137 <- get_arg(
  "--run137",
  read_latest("137_plot_global_PCA_and_PC1_core_structure")
)

KO_CLR_PATH <- get_arg(
  "--ko_clr",
  if (!is.na(RUN156)) {
    file.path(RUN156, "tables", "156_KO_CLR_network_filtered.tsv.gz")
  } else {
    NA_character_
  }
)
METADATA_PATH <- get_arg(
  "--metadata",
  if (!is.na(RUN156)) {
    file.path(RUN156, "tables", "156_metadata_aligned.tsv")
  } else {
    NA_character_
  }
)
KO_STATS_PATH <- get_arg(
  "--ko_stats",
  if (!is.na(RUN156)) {
    file.path(RUN156, "tables", "156_KO_network_filter_stats.tsv")
  } else {
    NA_character_
  }
)
KO_CORE_CLASS_PATH <- get_arg(
  "--ko_core_classes",
  if (!is.na(RUN153)) {
    file.path(RUN153, "tables", "153_KO_classification_primary_minloc3.tsv")
  } else {
    NA_character_
  }
)
KO_DICT_PATH <- get_arg(
  "--ko_dict",
  if (!is.na(RUN137)) {
    file.path(RUN137, "tables", "137_KO_function_dictionary_used.tsv")
  } else {
    NA_character_
  }
)

SEED <- as.integer(get_arg("--seed", "123"))
REUSE_EDGE_RUN <- get_arg("--reuse_edge_run", "")
N_BOOT <- as.integer(get_arg("--n_boot", "500"))
COMMUNITY_METHOD <- get_arg("--community_method", "leiden")
RESOLUTION <- as.numeric(get_arg("--resolution", "1"))
N_STARTS <- as.integer(get_arg("--n_starts", "50"))
PARTITION_SEED <- as.integer(get_arg("--partition_seed", "2001"))
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1L])
source(file.path(dirname(normalizePath(script_file)), "lib", "network_partition.R"))
MAX_KOS <- as.integer(get_arg("--max_kos", "1500"))
MAX_N_PER_STRATUM <- as.integer(get_arg("--max_n_per_stratum", "4"))
MIN_BOOTSTRAP_SAMPLES <- as.integer(
  get_arg("--min_bootstrap_samples", "60")
)
MIN_ABS_PROP <- as.numeric(get_arg("--min_abs_prop", "0"))
TARGET_EDGES_PER_BOOTSTRAP <- as.integer(
  get_arg("--target_edges_per_bootstrap", "2000")
)
CONSENSUS_STABILITY_MIN <- as.numeric(
  get_arg("--consensus_stability_min", "0.80")
)
SIGN_CONSISTENCY_MIN <- as.numeric(
  get_arg("--sign_consistency_min", "1.00")
)
CONSENSUS_MEAN_ABS_MIN <- as.numeric(
  get_arg("--consensus_mean_abs_min", "0.95")
)
MIN_MODULE_SIZE <- as.integer(get_arg("--min_module_size", "10"))
TOP_HUBS_PER_MODULE <- as.integer(
  get_arg("--top_hubs_per_module", "10")
)

REALM_COL_ARG <- get_arg("--realm_col", "auto")
COUNTRY_COL_ARG <- get_arg("--country_col", "auto")
LOCALITY_COL_ARG <- get_arg("--locality_col", "auto")

STRICT_EXPECTED <- as_flag(get_arg("--strict_expected", "1"))
ENFORCE_REFERENCE <- as_flag(get_arg("--enforce_reference", "0"))

EXPECTED_N_SAMPLES <- as.integer(get_arg("--expected_n_samples", "191"))
EXPECTED_TOTAL_KOS <- as.integer(get_arg("--expected_total_kos", "17482"))
EXPECTED_NETWORK_UNIVERSE <- as.integer(
  get_arg("--expected_network_universe", "5000")
)
EXPECTED_SELECTED_KOS <- as.integer(
  get_arg("--expected_selected_kos", "1500")
)
EXPECTED_BOOTSTRAP_MEDIAN_N <- as.numeric(
  get_arg("--expected_bootstrap_median_n", "81")
)
REFERENCE_CONSENSUS_EDGES <- as.integer(
  get_arg("--reference_consensus_edges", "1223")
)
REFERENCE_CONNECTED_NODES <- as.integer(
  get_arg("--reference_connected_nodes", "326")
)
REFERENCE_MAJOR_MODULES <- as.integer(
  get_arg("--reference_major_modules", "6")
)

SCRIPT_NUMBER <- "171"
SCRIPT_BASENAME <- "FUNCTIONAL_KO_CONSENSUS_NETWORK_FINAL"
STAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUTDIR <- file.path(
  RESULTS_DIR,
  paste0(SCRIPT_NUMBER, "_", SCRIPT_BASENAME, "_", STAMP)
)
TABLES_DIR <- file.path(OUTDIR, "tables")
LOGS_DIR <- file.path(OUTDIR, "logs")
RDS_DIR <- file.path(OUTDIR, "rds")
LATEST_FILE <- file.path(
  RESULTS_DIR,
  paste0("LATEST_", SCRIPT_NUMBER, "_", SCRIPT_BASENAME, ".txt")
)
STATUS_FILE <- file.path(OUTDIR, "STATUS.txt")

CLASS_LEVELS <- c(
  "absolute_global_core",
  "near_global",
  "realm_core_only",
  "multi_realm_non_core",
  "realm_specific"
)

dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOGS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)
writeLines(OUTDIR, LATEST_FILE)
writeLines("RUNNING", STATUS_FILE)

clean_chr <- function(x) {
  y <- trimws(as.character(x))
  y[
    is.na(y) | y == "" |
      tolower(y) %in% c("na", "n/a", "nan", "null", "none", "unknown")
  ] <- NA_character_
  y
}

safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

write_tsv <- function(x, path) {
  data.table::fwrite(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )
}

read_tsv <- function(path, required = TRUE) {
  if (is.na(path) || path == "" || !file.exists(path)) {
    if (required) stop("No existe archivo requerido: ", path)
    return(data.frame())
  }
  as.data.frame(
    data.table::fread(
      path,
      sep = "\t",
      header = TRUE,
      data.table = FALSE,
      check.names = FALSE,
      showProgress = FALSE
    )
  )
}

validate_run <- function(run_path, label) {
  if (is.na(run_path) || run_path == "" || !dir.exists(run_path)) {
    stop(label, " invalido: ", run_path)
  }
  status_path <- file.path(run_path, "STATUS.txt")
  if (!file.exists(status_path)) {
    stop(label, " no contiene STATUS.txt: ", run_path)
  }
  status <- trimws(readLines(status_path, warn = FALSE)[1L])
  if (!identical(status, "PASS")) {
    stop(label, " no tiene STATUS PASS: ", run_path)
  }
}

require_columns <- function(x, columns, label) {
  missing <- setdiff(columns, names(x))
  if (length(missing) > 0L) {
    stop("Faltan columnas en ", label, ": ", paste(missing, collapse = ", "))
  }
}

find_col <- function(x, candidates, required = TRUE, label = "column") {
  direct <- candidates[candidates %in% names(x)]
  if (length(direct) > 0L) return(direct[1L])
  lower_names <- tolower(names(x))
  for (candidate_i in candidates) {
    hit <- which(lower_names == tolower(candidate_i))
    if (length(hit) > 0L) return(names(x)[hit[1L]])
  }
  if (required) {
    stop(
      "No se encontro ", label, ". Candidatas: ",
      paste(candidates, collapse = ", ")
    )
  }
  NA_character_
}

detect_ko_col <- function(x) {
  direct <- names(x)[
    tolower(names(x)) %in% c("ko", "kegg_ko", "kegg_ko_id", "feature")
  ]
  if (length(direct) > 0L) return(direct[1L])
  candidates <- names(x)[vapply(x, function(column_i) {
    sum(grepl("^K[0-9]{5}$", clean_chr(column_i)), na.rm = TRUE) >= 10L
  }, logical(1))]
  if (length(candidates) == 0L) stop("No fue posible detectar la columna KO.")
  candidates[1L]
}

resolve_balance_column <- function(meta, argument, candidates, label) {
  if (!identical(argument, "auto")) {
    if (!argument %in% names(meta)) {
      stop("La columna solicitada para ", label, " no existe: ", argument)
    }
    return(argument)
  }
  find_col(meta, candidates, required = TRUE, label = label)
}

sanitize_stratum <- function(x) {
  y <- clean_chr(x)
  y[is.na(y)] <- "Unassigned"
  gsub("[^A-Za-z0-9_.-]+", "_", y)
}

prepare_dictionary <- function(path) {
  if (is.na(path) || path == "" || !file.exists(path)) {
    warning("No existe diccionario KO; se conservaran los IDs como etiquetas.")
    return(data.frame(KO = character(), KO_function = character()))
  }
  x <- read_tsv(path)
  ko_col <- detect_ko_col(x)
  names(x)[names(x) == ko_col] <- "KO"
  function_col <- find_col(
    x,
    c("KO_function", "function", "description", "definition", "name"),
    required = FALSE
  )
  if (is.na(function_col)) {
    x$KO_function <- NA_character_
  } else if (function_col != "KO_function") {
    names(x)[names(x) == function_col] <- "KO_function"
  }
  x$KO <- clean_chr(x$KO)
  x$KO_function <- clean_chr(x$KO_function)
  x <- x[grepl("^K[0-9]{5}$", x$KO), c("KO", "KO_function"), drop = FALSE]
  x[!duplicated(x$KO), , drop = FALSE]
}

build_sampling_bias_audit <- function(meta, columns) {
  output <- list()
  index <- 1L
  for (column_i in columns) {
    values <- clean_chr(meta[[column_i]])
    values[is.na(values)] <- "Unassigned"
    counts <- as.data.frame(table(values), stringsAsFactors = FALSE)
    names(counts) <- c("level", "n_samples")
    counts$variable <- column_i
    counts$n_levels <- length(unique(values))
    counts$max_level_n <- max(counts$n_samples)
    counts$min_level_n <- min(counts$n_samples)
    counts$imbalance_ratio_max_min <-
      counts$max_level_n / pmax(1, counts$min_level_n)
    output[[index]] <- counts[, c(
      "variable", "level", "n_samples", "n_levels", "max_level_n",
      "min_level_n", "imbalance_ratio_max_min"
    )]
    index <- index + 1L
  }
  do.call(rbind, output)
}

make_bootstrap_design <- function(
  meta,
  n_boot,
  max_n_per_stratum,
  seed
) {
  set.seed(seed)
  split_indices <- split(seq_len(nrow(meta)), meta$BALANCE_STRATUM)
  output <- vector("list", n_boot)

  for (bootstrap_i in seq_len(n_boot)) {
    selected <- unlist(lapply(split_indices, function(indices_i) {
      n_take <- min(length(indices_i), max_n_per_stratum)
      if (length(indices_i) <= n_take) {
        indices_i
      } else {
        sample(indices_i, size = n_take, replace = FALSE)
      }
    }), use.names = FALSE)
    selected <- unique(selected)
    output[[bootstrap_i]] <- data.frame(
      bootstrap_id = bootstrap_i,
      run_accession = meta$run_accession[selected],
      BALANCE_STRATUM = meta$BALANCE_STRATUM[selected],
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, output)
}

compute_proportionality_edges <- function(
  X,
  upper_indices,
  target_edges,
  min_abs_prop
) {
  n_samples <- nrow(X)
  if (n_samples < 5L) stop("Muy pocas muestras para proporcionalidad.")

  centered <- scale(X, center = TRUE, scale = FALSE)
  storage.mode(centered) <- "double"
  covariance <- crossprod(centered) / (n_samples - 1)
  variances <- diag(covariance)
  denominator <- outer(variances, variances, "+")
  proportionality <- (2 * covariance) / denominator
  proportionality[!is.finite(proportionality)] <- NA_real_
  diag(proportionality) <- NA_real_

  values <- proportionality[upper_indices]
  absolute_values <- abs(values)
  eligible <- is.finite(absolute_values) & absolute_values >= min_abs_prop

  if (!any(eligible)) {
    return(data.table::data.table(
      KO1 = character(), KO2 = character(), weight = numeric(),
      abs_weight = numeric(), sign = integer()
    ))
  }

  eligible_indices <- which(eligible)
  ordered <- order(
    absolute_values[eligible_indices],
    decreasing = TRUE,
    method = "radix"
  )
  keep_n <- min(target_edges, length(ordered))
  selected_vector_indices <- eligible_indices[ordered[seq_len(keep_n)]]
  selected_matrix_indices <- upper_indices[selected_vector_indices, , drop = FALSE]
  selected_values <- values[selected_vector_indices]

  data.table::data.table(
    KO1 = colnames(X)[selected_matrix_indices[, 1L]],
    KO2 = colnames(X)[selected_matrix_indices[, 2L]],
    weight = as.numeric(selected_values),
    abs_weight = as.numeric(abs(selected_values)),
    sign = ifelse(selected_values >= 0, 1L, -1L)
  )
}

class_count_table <- function(annotation, set_name) {
  counts <- as.data.frame(
    table(factor(annotation$functional_core_class, levels = CLASS_LEVELS)),
    stringsAsFactors = FALSE
  )
  names(counts) <- c("functional_core_class", "n_KOs")
  counts$set <- set_name
  counts$total_KOs <- sum(counts$n_KOs)
  counts$percent_KOs <- 100 * counts$n_KOs / pmax(1, counts$total_KOs)
  counts[, c(
    "set", "functional_core_class", "n_KOs", "total_KOs", "percent_KOs"
  )]
}

candidate_class_enrichment <- function(universe, selected) {
  output <- lapply(CLASS_LEVELS, function(class_i) {
    universe_in_class <- universe$functional_core_class == class_i
    universe_selected <- universe$KO %in% selected$KO
    contingency <- matrix(
      c(
        sum(universe_selected & universe_in_class),
        sum(universe_selected & !universe_in_class),
        sum(!universe_selected & universe_in_class),
        sum(!universe_selected & !universe_in_class)
      ),
      nrow = 2,
      byrow = TRUE
    )
    test <- stats::fisher.test(contingency)
    data.frame(
      functional_core_class = class_i,
      selected_in_class = contingency[1, 1],
      selected_not_class = contingency[1, 2],
      nonselected_in_class = contingency[2, 1],
      nonselected_not_class = contingency[2, 2],
      odds_ratio = unname(test$estimate),
      p_value = test$p.value,
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, output)
  result$q_value <- stats::p.adjust(result$p_value, method = "BH")
  result
}

module_class_enrichment <- function(nodes_connected) {
  modules <- sort(unique(nodes_connected$module))
  output <- list()
  index <- 1L
  for (module_i in modules) {
    in_module <- nodes_connected$module == module_i
    for (class_i in CLASS_LEVELS) {
      in_class <- nodes_connected$functional_core_class == class_i
      contingency <- matrix(
        c(
          sum(in_module & in_class),
          sum(in_module & !in_class),
          sum(!in_module & in_class),
          sum(!in_module & !in_class)
        ),
        nrow = 2,
        byrow = TRUE
      )
      test <- tryCatch(
        stats::fisher.test(contingency),
        error = function(error_i) NULL
      )
      output[[index]] <- data.frame(
        module = module_i,
        functional_core_class = class_i,
        in_module_in_class = contingency[1, 1],
        in_module_not_class = contingency[1, 2],
        outside_module_in_class = contingency[2, 1],
        outside_module_not_class = contingency[2, 2],
        odds_ratio = if (is.null(test)) NA_real_ else unname(test$estimate),
        p_value = if (is.null(test)) NA_real_ else test$p.value,
        stringsAsFactors = FALSE
      )
      index <- index + 1L
    }
  }
  result <- do.call(rbind, output)
  result$q_value <- stats::p.adjust(result$p_value, method = "BH")
  result
}

build_module_eigengenes <- function(Y, nodes, min_module_size) {
  module_sizes <- table(nodes$module[!is.na(nodes$module)])
  major_modules <- as.integer(names(module_sizes[module_sizes >= min_module_size]))
  scores <- data.frame(run_accession = rownames(Y), stringsAsFactors = FALSE)
  summaries <- list()
  index <- 1L

  for (module_i in sort(major_modules)) {
    kos <- nodes$KO[nodes$module == module_i]
    kos <- intersect(kos, colnames(Y))
    X <- Y[, kos, drop = FALSE]
    variable <- apply(X, 2, stats::sd) > 0
    X <- X[, variable, drop = FALSE]
    if (ncol(X) < min_module_size) next

    pca <- stats::prcomp(X, center = TRUE, scale. = TRUE)
    eigengene <- pca$x[, 1L]
    standardized_mean <- rowMeans(scale(X), na.rm = TRUE)
    if (stats::cor(eigengene, standardized_mean) < 0) eigengene <- -eigengene

    eigengene_name <- paste0("module_", module_i, "_ME1")
    scores[[eigengene_name]] <- as.numeric(eigengene)
    summaries[[index]] <- data.frame(
      module = module_i,
      eigengene = eigengene_name,
      n_KOs = ncol(X),
      variance_explained_PC1 = pca$sdev[1L]^2 / sum(pca$sdev^2),
      stringsAsFactors = FALSE
    )
    index <- index + 1L
  }

  list(
    scores = scores,
    summary = if (length(summaries) == 0L) data.frame() else do.call(rbind, summaries)
  )
}

main <- function() {
  validate_run(RUN156, "RUN156")
  # Explicit class tables supersede discovery through a historical run folder.
  if (is.na(KO_CORE_CLASS_PATH) || !file.exists(KO_CORE_CLASS_PATH)) {
    stop("No existe la tabla de clases funcionales: ", KO_CORE_CLASS_PATH)
  }

  if (
    N_BOOT < 1L || MAX_KOS < 2L || MAX_N_PER_STRATUM < 1L ||
      TARGET_EDGES_PER_BOOTSTRAP < 1L ||
      CONSENSUS_STABILITY_MIN < 0 || CONSENSUS_STABILITY_MIN > 1 ||
      SIGN_CONSISTENCY_MIN < 0 || SIGN_CONSISTENCY_MIN > 1 ||
      CONSENSUS_MEAN_ABS_MIN < 0 || CONSENSUS_MEAN_ABS_MIN > 1
  ) {
    stop("Parametros de red invalidos.")
  }

  log_path <- file.path(LOGS_DIR, "171_log.txt")
  sink(log_path, split = TRUE)
  on.exit({
    while (sink.number() > 0L) sink()
  }, add = TRUE)

  cat("============================================================\n")
  cat("171_FUNCTIONAL_KO_CONSENSUS_NETWORK_FINAL\n")
  cat("Started: ", as.character(Sys.time()), "\n", sep = "")
  cat("RUN156: ", RUN156, "\n", sep = "")
  cat("RUN153: ", RUN153, "\n", sep = "")
  cat("RUN137: ", RUN137, "\n", sep = "")
  cat("KO_CLR_PATH: ", KO_CLR_PATH, "\n", sep = "")
  cat("METADATA_PATH: ", METADATA_PATH, "\n", sep = "")
  cat("KO_STATS_PATH: ", KO_STATS_PATH, "\n", sep = "")
  cat("KO_CORE_CLASS_PATH: ", KO_CORE_CLASS_PATH, "\n", sep = "")
  cat("KO_DICT_PATH: ", KO_DICT_PATH, "\n", sep = "")
  cat("N_BOOT: ", N_BOOT, "\n", sep = "")
  cat("MAX_KOS: ", MAX_KOS, "\n", sep = "")
  cat("MAX_N_PER_STRATUM: ", MAX_N_PER_STRATUM, "\n", sep = "")
  cat("TARGET_EDGES_PER_BOOTSTRAP: ", TARGET_EDGES_PER_BOOTSTRAP, "\n", sep = "")
  cat("CONSENSUS_STABILITY_MIN: ", CONSENSUS_STABILITY_MIN, "\n", sep = "")
  cat("SIGN_CONSISTENCY_MIN: ", SIGN_CONSISTENCY_MIN, "\n", sep = "")
  cat("CONSENSUS_MEAN_ABS_MIN: ", CONSENSUS_MEAN_ABS_MIN, "\n", sep = "")
  cat("STRICT_EXPECTED: ", STRICT_EXPECTED, "\n", sep = "")
  cat("ENFORCE_REFERENCE: ", ENFORCE_REFERENCE, "\n", sep = "")
  cat("OUTDIR: ", OUTDIR, "\n", sep = "")
  cat("============================================================\n\n")

  cat("[1] Leyendo matriz CLR, metadata y clases funcionales canonicas...\n")
  clr_table <- read_tsv(KO_CLR_PATH)
  ko_col <- detect_ko_col(clr_table)
  names(clr_table)[names(clr_table) == ko_col] <- "KO"
  clr_table$KO <- clean_chr(clr_table$KO)
  clr_table <- clr_table[
    grepl("^K[0-9]{5}$", clr_table$KO) & !duplicated(clr_table$KO),
    ,
    drop = FALSE
  ]

  sample_columns <- setdiff(names(clr_table), "KO")
  matrix_KO_sample <- as.matrix(clr_table[, sample_columns, drop = FALSE])
  storage.mode(matrix_KO_sample) <- "double"
  rownames(matrix_KO_sample) <- clr_table$KO
  Y_universe <- t(matrix_KO_sample)
  rm(clr_table, matrix_KO_sample)
  gc()

  metadata <- read_tsv(METADATA_PATH)
  metadata_id_col <- find_col(
    metadata,
    c("run_accession", ".sample_name", "sample_id"),
    required = TRUE,
    label = "ID de metadata"
  )
  names(metadata)[names(metadata) == metadata_id_col] <- "run_accession"
  metadata$run_accession <- clean_chr(metadata$run_accession)
  if (anyNA(metadata$run_accession) || anyDuplicated(metadata$run_accession)) {
    stop("Metadata con run_accession faltante o duplicado.")
  }

  common_samples <- intersect(rownames(Y_universe), metadata$run_accession)
  Y_universe <- Y_universe[common_samples, , drop = FALSE]
  metadata <- metadata[
    match(common_samples, metadata$run_accession),
    ,
    drop = FALSE
  ]
  if (!identical(rownames(Y_universe), metadata$run_accession)) {
    stop("Fallo la alineacion por run_accession.")
  }
  if (any(!is.finite(Y_universe))) stop("La matriz CLR contiene valores no finitos.")

  KO_stats <- read_tsv(KO_STATS_PATH)
  ko_stats_col <- detect_ko_col(KO_stats)
  names(KO_stats)[names(KO_stats) == ko_stats_col] <- "KO"
  require_columns(
    KO_stats,
    c(
      "KO", "prevalence_n", "total_count",
      "CLR_variance", "network_rank_by_CLR_variance"
    ),
    "156_KO_network_filter_stats.tsv"
  )
  KO_stats$KO <- clean_chr(KO_stats$KO)
  KO_stats$prevalence_n <- safe_num(KO_stats$prevalence_n)
  KO_stats$total_count <- safe_num(KO_stats$total_count)
  KO_stats$CLR_variance <- safe_num(KO_stats$CLR_variance)
  KO_stats$network_rank_by_CLR_variance <- safe_num(
    KO_stats$network_rank_by_CLR_variance
  )
  if (anyDuplicated(KO_stats$KO)) stop("KOs duplicados en KO_stats.")

  core_classes <- read_tsv(KO_CORE_CLASS_PATH)
  require_columns(
    core_classes,
    c(
      "KO", "functional_core_class", "global_prevalence_n",
      "global_prevalence_prop", "total_abundance",
      "min_realm_localities", "min_mean_relative_abundance_realm_core"
    ),
    "153_KO_classification_primary_minloc3.tsv"
  )
  core_classes$KO <- clean_chr(core_classes$KO)
  core_classes$functional_core_class <- clean_chr(
    core_classes$functional_core_class
  )
  if (anyDuplicated(core_classes$KO)) stop("KOs duplicados en clases 153.")
  if (STRICT_EXPECTED && nrow(core_classes) != EXPECTED_TOTAL_KOS) {
    stop(
      "La clasificacion 153 debe contener ", EXPECTED_TOTAL_KOS,
      " KOs; contiene ", nrow(core_classes), "."
    )
  }
  if ("n_functional_samples" %in% names(core_classes)) {
    n_functional_samples <- unique(safe_num(core_classes$n_functional_samples))
    n_functional_samples <- n_functional_samples[is.finite(n_functional_samples)]
    if (length(n_functional_samples) != 1L ||
        n_functional_samples != EXPECTED_N_SAMPLES) {
      stop("La clasificacion 153 no corresponde a ", EXPECTED_N_SAMPLES,
           " metagenomas.")
    }
  }
  if (any(!core_classes$functional_core_class %in% CLASS_LEVELS)) {
    stop("Clases funcionales inesperadas en el output 153.")
  }
  if (data.table::uniqueN(core_classes$min_realm_localities) != 1L ||
      unique(core_classes$min_realm_localities) != 3) {
    stop("La clasificacion 153 no corresponde al minimo primario de 3 localidades.")
  }
  near_threshold <- unique(safe_num(core_classes$near_global_threshold))
  near_threshold <- near_threshold[is.finite(near_threshold)]
  realm_abundance_threshold <- unique(safe_num(
    core_classes$min_mean_relative_abundance_realm_core
  ))
  realm_abundance_threshold <- realm_abundance_threshold[
    is.finite(realm_abundance_threshold)
  ]
  if (length(near_threshold) != 1L || abs(near_threshold - 0.95) > 1e-12 ||
      length(realm_abundance_threshold) != 1L ||
      abs(realm_abundance_threshold - 0.0001) > 1e-12) {
    stop("Los umbrales funcionales de 153 no coinciden con el manuscrito.")
  }

  legacy_class_cols <- intersect(
    c("functional_core_class", "functional_prevalence_class"),
    names(KO_stats)
  )
  if (length(legacy_class_cols) > 0L) {
    KO_stats <- KO_stats[, setdiff(names(KO_stats), legacy_class_cols), drop = FALSE]
  }
  core_keep <- unique(c(
    "KO", "functional_core_class", "global_prevalence_n",
    "global_prevalence_prop", "n_functional_samples", "total_abundance",
    "n_realms_detected", "realms_detected", "realm_core_realms",
    "is_absolute_global_core", "is_near_global", "is_realm_core_only",
    "is_multi_realm_non_core", "is_realm_specific",
    "near_global_threshold", "min_realm_localities",
    "min_mean_relative_abundance_realm_core", "presence_threshold",
    "abundance_weighting"
  ))
  core_keep <- intersect(core_keep, names(core_classes))
  KO_stats <- merge(
    KO_stats,
    core_classes[, core_keep, drop = FALSE],
    by = "KO",
    all.x = TRUE,
    sort = FALSE
  )
  if (anyNA(KO_stats$functional_core_class)) {
    stop("Hay KOs del universo de red sin clase funcional canonica de 153.")
  }

  missing_stats <- setdiff(colnames(Y_universe), KO_stats$KO)
  missing_matrix <- setdiff(KO_stats$KO, colnames(Y_universe))
  if (length(missing_stats) > 0L || length(missing_matrix) > 0L) {
    stop("El universo KO no coincide exactamente entre matriz y anotacion.")
  }
  KO_stats <- KO_stats[match(colnames(Y_universe), KO_stats$KO), , drop = FALSE]

  cat("    samples: ", nrow(Y_universe), "\n", sep = "")
  cat("    network-universe KOs: ", ncol(Y_universe), "\n", sep = "")

  cat("[2] Seleccionando 1,500 KOs por varianza CLR...\n")
  observed_variance <- apply(Y_universe, 2, stats::var)
  observed_variance[!is.finite(observed_variance)] <- 0
  selection_order <- order(
    -observed_variance,
    names(observed_variance)
  )
  selected_n <- min(MAX_KOS, length(selection_order))
  selected_KOs <- names(observed_variance)[selection_order[seq_len(selected_n)]]
  Y <- Y_universe[, selected_KOs, drop = FALSE]

  universe_annotation <- KO_stats
  candidate_annotation <- KO_stats[
    match(selected_KOs, KO_stats$KO),
    ,
    drop = FALSE
  ]
  candidate_annotation$selection_rank_by_recalculated_CLR_variance <-
    seq_len(nrow(candidate_annotation))
  candidate_annotation$recalculated_CLR_variance <- observed_variance[selected_KOs]
  candidate_annotation$CLR_variance_absolute_difference <- abs(
    candidate_annotation$recalculated_CLR_variance - candidate_annotation$CLR_variance
  )

  dictionary <- prepare_dictionary(KO_DICT_PATH)
  candidate_annotation$KO_function <- dictionary$KO_function[
    match(candidate_annotation$KO, dictionary$KO)
  ]
  candidate_annotation$KO_function[
    is.na(candidate_annotation$KO_function)
  ] <- candidate_annotation$KO[
    is.na(candidate_annotation$KO_function)
  ]

  class_counts <- rbind(
    class_count_table(core_classes, paste0("all_detected_KOs_", nrow(core_classes))),
    class_count_table(universe_annotation, "network_universe_5000"),
    class_count_table(candidate_annotation, "selected_candidates_1500")
  )
  write_tsv(
    universe_annotation,
    file.path(TABLES_DIR, "171_network_universe_KO_annotation.tsv")
  )
  selection_enrichment <- candidate_class_enrichment(
    universe_annotation,
    candidate_annotation
  )
  write_tsv(
    candidate_annotation,
    file.path(TABLES_DIR, "171_candidate_KO_selection.tsv")
  )
  write_tsv(
    selection_enrichment,
    file.path(TABLES_DIR, "171_candidate_class_selection_enrichment.tsv")
  )

  cat("    selected KOs: ", ncol(Y), "\n", sep = "")
  print(class_counts[class_counts$set == "selected_candidates_1500", ])

  cat("[3] Auditando y construyendo el diseno balanceado...\n")
  realm_col <- resolve_balance_column(
    metadata,
    REALM_COL_ARG,
    c("REALM_min5", "REALM", "realm", "MarineRealm", "marine_realm"),
    "marine realm"
  )
  country_col <- resolve_balance_column(
    metadata,
    COUNTRY_COL_ARG,
    c("Country", "country"),
    "country"
  )
  locality_col <- resolve_balance_column(
    metadata,
    LOCALITY_COL_ARG,
    c("Locality", "locality"),
    "locality"
  )
  balance_columns <- c(realm_col, country_col, locality_col)
  metadata$BALANCE_STRATUM <- paste(
    sanitize_stratum(metadata[[realm_col]]),
    sanitize_stratum(metadata[[country_col]]),
    sanitize_stratum(metadata[[locality_col]]),
    sep = "__"
  )

  sampling_audit <- build_sampling_bias_audit(metadata, balance_columns)
  stratum_counts <- as.data.frame(
    table(metadata$BALANCE_STRATUM),
    stringsAsFactors = FALSE
  )
  names(stratum_counts) <- c("BALANCE_STRATUM", "n_samples")
  stratum_counts$n_retained_per_bootstrap <- pmin(
    stratum_counts$n_samples,
    MAX_N_PER_STRATUM
  )
  write_tsv(
    sampling_audit,
    file.path(TABLES_DIR, "171_sampling_bias_audit.tsv")
  )
  write_tsv(
    stratum_counts,
    file.path(TABLES_DIR, "171_balance_stratum_counts.tsv")
  )

  bootstrap_design <- make_bootstrap_design(
    metadata,
    N_BOOT,
    MAX_N_PER_STRATUM,
    SEED
  )
  bootstrap_n <- aggregate(
    run_accession ~ bootstrap_id,
    data = bootstrap_design,
    FUN = length
  )
  names(bootstrap_n)[2L] <- "n_samples"
  data.table::fwrite(
    bootstrap_design,
    file.path(TABLES_DIR, "171_bootstrap_design.tsv.gz"),
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )
  write_tsv(
    bootstrap_n,
    file.path(TABLES_DIR, "171_bootstrap_sample_sizes.tsv")
  )

  if (min(bootstrap_n$n_samples) < MIN_BOOTSTRAP_SAMPLES) {
    stop("El diseno balanceado produjo menos muestras que el minimo permitido.")
  }
  cat(
    "    balance columns: ", paste(balance_columns, collapse = ";"), "\n",
    sep = ""
  )
  cat("    balance strata: ", nrow(stratum_counts), "\n", sep = "")
  cat(
    "    bootstrap n median [range]: ",
    stats::median(bootstrap_n$n_samples), " [",
    min(bootstrap_n$n_samples), "-", max(bootstrap_n$n_samples), "]\n",
    sep = ""
  )

  cat("[4] Infiriendo aristas de proporcionalidad en ", N_BOOT, " submuestreos...\n",sep="")
  upper_indices <- which(
    upper.tri(matrix(FALSE, nrow = ncol(Y), ncol = ncol(Y)), diag = FALSE),
    arr.ind = TRUE
  )
  if(nzchar(REUSE_EDGE_RUN)) {
    cached_tables <- file.path(REUSE_EDGE_RUN,"tables")
    cached_design <- data.table::fread(file.path(cached_tables,"171_bootstrap_design.tsv.gz"))
    stopifnot(identical(as.data.frame(cached_design),bootstrap_design))
    cached_candidates <- read_tsv(file.path(cached_tables,"171_candidate_KO_selection.tsv"))
    stopifnot(identical(cached_candidates$KO,candidate_annotation$KO))
    cached_params <- read_tsv(file.path(cached_tables,"171_parameters.tsv"))
    for(key in c("ko_clr_path","metadata_path","ko_stats_path")) {
      cached_path <- cached_params$value[cached_params$parameter==key]
      current_path <- switch(key,ko_clr_path=KO_CLR_PATH,metadata_path=METADATA_PATH,ko_stats_path=KO_STATS_PATH)
      stopifnot(length(cached_path)==1L, file.exists(cached_path),
        unname(tools::md5sum(cached_path))==unname(tools::md5sum(current_path)))
    }
    all_edges <- data.table::fread(file.path(cached_tables,"171_bootstrap_edge_membership.tsv.gz"))
    stopifnot(data.table::uniqueN(all_edges$bootstrap_id)==N_BOOT)
    edge_counts <- as.integer(table(factor(all_edges$bootstrap_id,levels=seq_len(N_BOOT))))
    stopifnot(all(edge_counts==TARGET_EDGES_PER_BOOTSTRAP))
    writeLines(REUSE_EDGE_RUN,file.path(OUTDIR,"REUSED_EDGE_RUN.txt"))
  } else {
  edge_list <- vector("list", N_BOOT)
  edge_counts <- integer(N_BOOT)

  for (bootstrap_i in seq_len(N_BOOT)) {
    sample_ids <- bootstrap_design$run_accession[
      bootstrap_design$bootstrap_id == bootstrap_i
    ]
    X <- Y[sample_ids, , drop = FALSE]
    edges_i <- compute_proportionality_edges(
      X,
      upper_indices,
      TARGET_EDGES_PER_BOOTSTRAP,
      MIN_ABS_PROP
    )
    edges_i$bootstrap_id <- bootstrap_i
    edge_list[[bootstrap_i]] <- edges_i
    edge_counts[bootstrap_i] <- nrow(edges_i)

    if (bootstrap_i == 1L || bootstrap_i %% 5L == 0L || bootstrap_i == N_BOOT) {
      cat(
        "    bootstrap ", bootstrap_i, "/", N_BOOT,
        ": n=", nrow(X), ", edges=", nrow(edges_i), "\n",
        sep = ""
      )
    }
    rm(X, edges_i)
    if (bootstrap_i %% 10L == 0L) gc()
  }

  all_edges <- data.table::rbindlist(edge_list, use.names = TRUE, fill = TRUE)
  rm(edge_list)
  gc()
  }
  if (nrow(all_edges) == 0L) stop("No se recuperaron aristas.")

  data.table::fwrite(
    all_edges,
    file.path(TABLES_DIR, "171_bootstrap_edge_membership.tsv.gz"),
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )
  edge_count_summary <- data.frame(
    bootstrap_id = seq_len(N_BOOT),
    n_edges = edge_counts,
    stringsAsFactors = FALSE
  )
  write_tsv(
    edge_count_summary,
    file.path(TABLES_DIR, "171_bootstrap_edge_counts.tsv")
  )

  cat("[5] Construyendo la red consenso estricta...\n")
  edge_stability <- all_edges[
    ,
    .(
      n_boot_present = data.table::uniqueN(bootstrap_id),
      edge_stability = data.table::uniqueN(bootstrap_id) / N_BOOT,
      mean_weight = mean(weight),
      median_weight = stats::median(weight),
      mean_abs_weight = mean(abs_weight),
      median_abs_weight = stats::median(abs_weight),
      positive_fraction = mean(sign > 0),
      negative_fraction = mean(sign < 0)
    ),
    by = .(KO1, KO2)
  ]
  edge_stability[, sign_consistency := pmax(positive_fraction, negative_fraction)]
  edge_stability[, consensus_sign := ifelse(
    positive_fraction >= negative_fraction, 1L, -1L
  )]
  edge_stability[, consensus_weight := consensus_sign * mean_abs_weight]
  data.table::setorder(
    edge_stability,
    -edge_stability,
    -sign_consistency,
    -mean_abs_weight,
    KO1,
    KO2
  )

  consensus_edges <- edge_stability[
    edge_stability >= CONSENSUS_STABILITY_MIN &
      sign_consistency >= SIGN_CONSISTENCY_MIN &
      mean_abs_weight >= CONSENSUS_MEAN_ABS_MIN
  ]
  if (nrow(consensus_edges) == 0L) {
    stop("La configuracion estricta produjo cero aristas consenso.")
  }

  data.table::fwrite(
    edge_stability,
    file.path(TABLES_DIR, "171_edge_stability_all.tsv.gz"),
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )
  data.table::fwrite(
    consensus_edges,
    file.path(TABLES_DIR, "171_consensus_edges.tsv.gz"),
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )

  cat("    candidate unique edges: ", nrow(edge_stability), "\n", sep = "")
  cat("    consensus edges: ", nrow(consensus_edges), "\n", sep = "")

  cat("[6] Construyendo grafo, modulos y arquitectura de prevalencia...\n")
  all_candidate_KOs <- candidate_annotation$KO
  edge_graph_table <- as.data.frame(
    consensus_edges[, .(
      KO1, KO2, consensus_weight, mean_abs_weight,
      edge_stability, sign_consistency
    )]
  )
  graph <- igraph::graph_from_data_frame(
    edge_graph_table,
    directed = FALSE,
    vertices = data.frame(name = all_candidate_KOs, stringsAsFactors = FALSE)
  )
  igraph::E(graph)$weight_abs <- igraph::E(graph)$mean_abs_weight
  components_full <- igraph::components(graph)
  partition <- fit_network_partition(graph, method=COMMUNITY_METHOD, resolution=RESOLUTION,
    seeds=PARTITION_SEED + seq_len(N_STARTS)-1L, weights=igraph::E(graph)$weight_abs)
  community_fit <- partition$community
  write_tsv(partition$fits, file.path(TABLES_DIR,"171_partition_multistart.tsv"))
  saveRDS(partition, file.path(RDS_DIR,"171_partition_multistart.rds"))

  all_node_metrics <- data.frame(
    KO = igraph::V(graph)$name,
    degree = igraph::degree(graph),
    strength_abs = igraph::strength(
      graph,
      weights = igraph::E(graph)$weight_abs
    ),
    component = as.integer(components_full$membership),
    module = as.integer(igraph::membership(community_fit)),
    stringsAsFactors = FALSE
  )

  nodes <- candidate_annotation
  metric_index <- match(nodes$KO, all_node_metrics$KO)
  nodes$degree <- all_node_metrics$degree[metric_index]
  nodes$strength_abs <- all_node_metrics$strength_abs[metric_index]
  nodes$component <- all_node_metrics$component[metric_index]
  nodes$module <- all_node_metrics$module[metric_index]
  nodes$degree[is.na(nodes$degree)] <- 0
  nodes$strength_abs[is.na(nodes$strength_abs)] <- 0
  nodes$connected <- nodes$degree > 0

  connected_graph <- igraph::induced_subgraph(
    graph,
    vids = igraph::V(graph)[igraph::degree(graph) > 0]
  )
  connected_components <- igraph::components(connected_graph)

  module_sizes <- as.data.frame(
    table(nodes$module[!is.na(nodes$module)]),
    stringsAsFactors = FALSE
  )
  names(module_sizes) <- c("module", "n_KOs")
  module_sizes$module <- as.integer(as.character(module_sizes$module))
  major_modules <- module_sizes$module[module_sizes$n_KOs >= MIN_MODULE_SIZE]
  nodes$major_module <- !is.na(nodes$module) & nodes$module %in% major_modules

  connected_annotation <- nodes[nodes$connected, , drop = FALSE]
  major_annotation <- nodes[nodes$major_module, , drop = FALSE]
  class_counts <- rbind(
    class_counts,
    class_count_table(connected_annotation, "consensus_connected_nodes"),
    class_count_table(major_annotation, "major_modules_nodes")
  )

  module_summaries <- lapply(sort(unique(connected_annotation$module)), function(module_i) {
    module_nodes <- connected_annotation[
      connected_annotation$module == module_i,
      ,
      drop = FALSE
    ]
    class_n <- table(factor(
      module_nodes$functional_core_class,
      levels = CLASS_LEVELS
    ))
    data.frame(
      module = module_i,
      n_KOs = nrow(module_nodes),
      major_module = nrow(module_nodes) >= MIN_MODULE_SIZE,
      n_edges_internal = sum(
        consensus_edges$KO1 %in% module_nodes$KO &
          consensus_edges$KO2 %in% module_nodes$KO
      ),
      mean_degree = mean(module_nodes$degree),
      mean_strength_abs = mean(module_nodes$strength_abs),
      n_absolute_global_core = unname(class_n["absolute_global_core"]),
      n_near_global = unname(class_n["near_global"]),
      n_realm_core_only = unname(class_n["realm_core_only"]),
      n_multi_realm_non_core = unname(class_n["multi_realm_non_core"]),
      n_realm_specific = unname(class_n["realm_specific"]),
      percent_absolute_global_core = 100 * unname(class_n["absolute_global_core"]) / nrow(module_nodes),
      percent_near_global = 100 * unname(class_n["near_global"]) / nrow(module_nodes),
      percent_realm_core_only = 100 * unname(class_n["realm_core_only"]) / nrow(module_nodes),
      percent_multi_realm_non_core = 100 * unname(class_n["multi_realm_non_core"]) / nrow(module_nodes),
      percent_realm_specific = 100 * unname(class_n["realm_specific"]) / nrow(module_nodes),
      dominant_core_class = names(which.max(class_n)),
      stringsAsFactors = FALSE
    )
  })
  module_summary <- do.call(rbind, module_summaries)
  module_summary <- module_summary[order(-module_summary$n_KOs), , drop = FALSE]
  module_enrichment <- module_class_enrichment(connected_annotation)

  node_class <- setNames(
    nodes$functional_core_class,
    nodes$KO
  )
  edge_architecture <- as.data.frame(consensus_edges)
  edge_architecture$class_KO1 <- node_class[edge_architecture$KO1]
  edge_architecture$class_KO2 <- node_class[edge_architecture$KO2]
  edge_architecture$class_index_KO1 <- match(
    edge_architecture$class_KO1,
    CLASS_LEVELS
  )
  edge_architecture$class_index_KO2 <- match(
    edge_architecture$class_KO2,
    CLASS_LEVELS
  )
  swap <- edge_architecture$class_index_KO1 > edge_architecture$class_index_KO2
  class_left <- ifelse(
    swap,
    edge_architecture$class_KO2,
    edge_architecture$class_KO1
  )
  class_right <- ifelse(
    swap,
    edge_architecture$class_KO1,
    edge_architecture$class_KO2
  )
  edge_architecture$class_pair <- paste(class_left, class_right, sep = "--")

  edge_class_mixing <- do.call(
    rbind,
    lapply(split(edge_architecture, edge_architecture$class_pair), function(x) {
      data.frame(
        class_pair = x$class_pair[1L],
        n_edges = nrow(x),
        percent_edges = 100 * nrow(x) / nrow(edge_architecture),
        n_positive = sum(x$consensus_weight > 0),
        n_negative = sum(x$consensus_weight < 0),
        mean_abs_weight = mean(x$mean_abs_weight),
        mean_stability = mean(x$edge_stability),
        stringsAsFactors = FALSE
      )
    })
  )
  edge_class_mixing <- edge_class_mixing[
    order(-edge_class_mixing$n_edges),
    ,
    drop = FALSE
  ]

  class_codes <- match(
    node_class[igraph::V(graph)$name],
    CLASS_LEVELS
  )
  class_assortativity <- igraph::assortativity_nominal(
    graph,
    types = class_codes,
    directed = FALSE
  )
  candidate_graph_density <- nrow(consensus_edges) /
    choose(nrow(candidate_annotation), 2L)
  connected_graph_density <- igraph::edge_density(
    connected_graph,
    loops = FALSE
  )

  global_metrics <- data.frame(
    metric = c(
      "n_candidate_KOs",
      "n_connected_KOs",
      "n_isolated_candidate_KOs",
      "n_consensus_edges",
      "edge_density_all_candidate_KOs",
      "edge_density_connected_graph",
      "n_components_among_connected_KOs",
      "largest_connected_component_size",
      "n_community_modules_connected",
      "n_major_modules_ge_min_size",
      "n_KOs_in_major_modules",
      "modularity_primary",
      "core_class_assortativity",
      "mean_degree_connected",
      "mean_strength_abs_connected"
    ),
    value = c(
      nrow(nodes),
      nrow(connected_annotation),
      sum(!nodes$connected),
      nrow(consensus_edges),
      candidate_graph_density,
      connected_graph_density,
      connected_components$no,
      max(connected_components$csize),
      length(unique(connected_annotation$module)),
      length(major_modules),
      nrow(major_annotation),
      partition$Q,
      class_assortativity,
      mean(connected_annotation$degree),
      mean(connected_annotation$strength_abs)
    ),
    stringsAsFactors = FALSE
  )

  module_eigengenes <- build_module_eigengenes(
    Y,
    nodes,
    MIN_MODULE_SIZE
  )
  hub_table <- nodes[nodes$major_module, , drop = FALSE]
  hub_table <- do.call(
    rbind,
    lapply(split(hub_table, hub_table$module), function(x) {
      x <- x[order(-x$strength_abs, -x$degree, x$KO), , drop = FALSE]
      x$module_hub_rank <- seq_len(nrow(x))
      x
    })
  )
  top_hubs <- hub_table[
    hub_table$module_hub_rank <= TOP_HUBS_PER_MODULE,
    ,
    drop = FALSE
  ]

  write_tsv(nodes, file.path(TABLES_DIR, "171_consensus_nodes.tsv"))
  write_tsv(class_counts, file.path(TABLES_DIR, "171_core_class_counts.tsv"))
  write_tsv(global_metrics, file.path(TABLES_DIR, "171_network_global_metrics.tsv"))
  write_tsv(module_summary, file.path(TABLES_DIR, "171_module_summary.tsv"))
  write_tsv(
    module_enrichment,
    file.path(TABLES_DIR, "171_module_core_class_enrichment.tsv")
  )
  write_tsv(
    edge_class_mixing,
    file.path(TABLES_DIR, "171_edge_core_class_mixing.tsv")
  )
  write_tsv(
    edge_architecture,
    file.path(TABLES_DIR, "171_consensus_edges_functional_core_annotated.tsv")
  )
  write_tsv(hub_table, file.path(TABLES_DIR, "171_major_module_hubs_all.tsv"))
  write_tsv(top_hubs, file.path(TABLES_DIR, "171_top_hubs_annotated.tsv"))
  write_tsv(
    module_eigengenes$scores,
    file.path(TABLES_DIR, "171_major_module_eigengenes.tsv")
  )
  write_tsv(
    module_eigengenes$summary,
    file.path(TABLES_DIR, "171_major_module_eigengene_summary.tsv")
  )

  cat("[7] Validaciones estrictas y comparacion con referencia...\n")
  consensus_duplicate_keys <- anyDuplicated(
    paste(consensus_edges$KO1, consensus_edges$KO2, sep = "||")
  )
  missing_consensus_endpoints <- length(setdiff(
    unique(c(consensus_edges$KO1, consensus_edges$KO2)),
    candidate_annotation$KO
  ))
  consensus_filter_violations <- sum(
    consensus_edges$edge_stability < CONSENSUS_STABILITY_MIN |
      consensus_edges$sign_consistency < SIGN_CONSISTENCY_MIN |
      consensus_edges$mean_abs_weight < CONSENSUS_MEAN_ABS_MIN
  )
  max_variance_difference <- max(
    candidate_annotation$CLR_variance_absolute_difference,
    na.rm = TRUE
  )

  validation <- data.frame(
    metric = c(
      "n_samples",
      "n_network_universe_KOs",
      "n_selected_candidate_KOs",
      "n_bootstraps",
      "median_bootstrap_n_samples",
      "min_bootstrap_n_samples",
      "min_edges_retained_per_bootstrap",
      "max_edges_retained_per_bootstrap",
      "n_missing_core_classes",
      "n_unexpected_core_classes",
      "n_missing_consensus_endpoints",
      "n_duplicate_consensus_edge_keys",
      "n_consensus_filter_violations",
      "n_nonfinite_CLR_values",
      "n_nonfinite_consensus_weights",
      "max_CLR_variance_difference_vs_script156"
    ),
    observed = c(
      nrow(Y),
      ncol(Y_universe),
      ncol(Y),
      N_BOOT,
      stats::median(bootstrap_n$n_samples),
      min(bootstrap_n$n_samples),
      min(edge_counts),
      max(edge_counts),
      sum(is.na(candidate_annotation$functional_core_class)),
      sum(!candidate_annotation$functional_core_class %in% CLASS_LEVELS),
      missing_consensus_endpoints,
      as.integer(consensus_duplicate_keys > 0),
      consensus_filter_violations,
      sum(!is.finite(Y)),
      sum(!is.finite(consensus_edges$consensus_weight)),
      max_variance_difference
    ),
    expected = c(
      EXPECTED_N_SAMPLES,
      EXPECTED_NETWORK_UNIVERSE,
      EXPECTED_SELECTED_KOS,
      N_BOOT,
      EXPECTED_BOOTSTRAP_MEDIAN_N,
      EXPECTED_BOOTSTRAP_MEDIAN_N,
      TARGET_EDGES_PER_BOOTSTRAP,
      TARGET_EDGES_PER_BOOTSTRAP,
      0, 0, 0, 0, 0, 0, 0, 0
    ),
    tolerance = c(
      rep(0, 15),
      1e-10
    ),
    stringsAsFactors = FALSE
  )
  validation$pass <- abs(validation$observed - validation$expected) <=
    validation$tolerance
  write_tsv(validation, file.path(TABLES_DIR, "171_expected_validation.tsv"))

  reference_check <- data.frame(
    metric = c(
      "n_consensus_edges",
      "n_connected_nodes",
      "n_major_modules"
    ),
    observed = c(
      nrow(consensus_edges),
      nrow(connected_annotation),
      length(major_modules)
    ),
    reference = c(
      REFERENCE_CONSENSUS_EDGES,
      REFERENCE_CONNECTED_NODES,
      REFERENCE_MAJOR_MODULES
    ),
    stringsAsFactors = FALSE
  )
  reference_check$matches_reference <-
    reference_check$observed == reference_check$reference
  write_tsv(
    reference_check,
    file.path(TABLES_DIR, "171_reference_reproduction_check.tsv")
  )

  if (STRICT_EXPECTED && any(!validation$pass)) {
    stop(
      "Fallo la validacion estricta: ",
      paste(validation$metric[!validation$pass], collapse = ", ")
    )
  }
  if (ENFORCE_REFERENCE && any(!reference_check$matches_reference)) {
    stop("La red no reprodujo los conteos de referencia solicitados.")
  }

  parameters <- data.frame(
    parameter = c(
      "method",
      "edge_interpretation",
      "run156",
      "run153",
      "run137",
      "ko_clr_path",
      "metadata_path",
      "ko_stats_path",
      "ko_core_class_path",
      "ko_dictionary_path",
      "canonical_sample_id",
      "core_classes",
      "network_KO_filter",
      "candidate_selection",
      "seed",
      "n_boot",
      "max_kos",
      "max_n_per_stratum",
      "balance_columns",
      "min_bootstrap_samples",
      "min_abs_prop",
      "target_edges_per_bootstrap",
      "consensus_stability_min",
      "sign_consistency_min",
      "consensus_mean_abs_min",
      "min_module_size",
      "top_hubs_per_module"
    ),
    value = c(
      "CLR proportionality consensus network with balanced geographic subsampling",
      "Robust co-abundance association; not biochemical interaction or causality",
      RUN156,
      RUN153,
      RUN137,
      KO_CLR_PATH,
      METADATA_PATH,
      KO_STATS_PATH,
      KO_CORE_CLASS_PATH,
      KO_DICT_PATH,
      "run_accession",
      paste(
        "absolute_global_core=100%; near_global>=95% and <100%;",
        "realm_core_only=all localities in >=1 realm with >=3 localities",
        "and mean realm-locality relative abundance>=0.0001;",
        "multi_realm_non_core; realm_specific"
      ),
      "prevalence_n>=20 and total_count>=50 inherited from script 156",
      "top CLR variance; classes not used for selection",
      SEED,
      N_BOOT,
      MAX_KOS,
      MAX_N_PER_STRATUM,
      paste(balance_columns, collapse = ";"),
      MIN_BOOTSTRAP_SAMPLES,
      MIN_ABS_PROP,
      TARGET_EDGES_PER_BOOTSTRAP,
      CONSENSUS_STABILITY_MIN,
      SIGN_CONSISTENCY_MIN,
      CONSENSUS_MEAN_ABS_MIN,
      MIN_MODULE_SIZE,
      TOP_HUBS_PER_MODULE
    ),
    stringsAsFactors = FALSE
  )
  parameters <- rbind(parameters, data.frame(parameter=c("community_method","resolution","n_starts","partition_seed_first","best_partition_seed","best_modularity"),
    value=c(COMMUNITY_METHOD,RESOLUTION,N_STARTS,PARTITION_SEED,partition$best_seed,partition$Q)))
  write_tsv(parameters, file.path(TABLES_DIR, "171_parameters.tsv"))

  methods_text <- c(
    "Functional co-abundance network inference",
    "",
    paste0(
      "We inferred a composition-aware functional co-abundance network from ",
      "CLR-transformed KEGG orthologue profiles across ", nrow(Y),
      " mangrove sediment metagenomes."
    ),
    paste0(
      "To limit geographic sampling imbalance, ", N_BOOT,
      " repeated datasets were generated by stratified subsampling without ",
      "replacement across ",
      "marine realm, country and locality, retaining at most ",
      MAX_N_PER_STRATUM, " samples per stratum."
    ),
    paste0(
      "The repeated subsamples contained a median of ",
      stats::median(bootstrap_n$n_samples), " metagenomes."
    ),
    paste0(
      "Pairwise associations among ", ncol(Y),
      " candidate KOs were estimated using proportionality in CLR space, ",
      "rho = 2cov(i,j)/(var(i)+var(j))."
    ),
    paste0(
      "For each bootstrap replicate, the ", TARGET_EDGES_PER_BOOTSTRAP,
      " strongest associations were retained."
    ),
    paste0(
      "Consensus edges required bootstrap stability >=",
      CONSENSUS_STABILITY_MIN, ", sign consistency >=",
      SIGN_CONSISTENCY_MIN, " and mean absolute association >=",
      CONSENSUS_MEAN_ABS_MIN, "."
    ),
    paste0(
      "Modules were identified by ", COMMUNITY_METHOD, " community detection at resolution ", RESOLUTION, "; the highest modularity among ",N_STARTS," starts was retained using absolute ",
      "consensus weights; biological interpretation was restricted to modules ",
      "containing at least ", MIN_MODULE_SIZE, " KOs."
    ),
    paste0(
      "KOs were annotated using the final hierarchical functional-persistence ",
      "classification: absolute global core, near-global, realm-core only, ",
      "multi-realm non-core and realm-specific."
    ),
    paste0(
      "Realm-core support required occurrence in every locality of at least ",
      "one marine realm represented by at least three localities and a mean ",
      "realm-locality relative abundance >=0.0001."
    ),
    "Edges represent reproducible co-abundance associations and not direct metabolic interactions."
  )
  writeLines(methods_text, file.path(TABLES_DIR, "171_methods_draft.txt"))

  summary_df <- data.frame(
    item = c(
      "status",
      "run156",
      "run153",
      "run137",
      "canonical_sample_id",
      "n_samples",
      "n_network_universe_KOs",
      "n_selected_candidate_KOs",
      "n_boot",
      "median_bootstrap_n_samples",
      "balance_columns",
      "n_balance_strata",
      "n_unique_candidate_edges",
      "n_consensus_edges",
      "n_connected_nodes",
      "edge_density_all_candidate_KOs",
      "edge_density_connected_graph",
      "n_major_modules",
      "n_KOs_in_major_modules",
      "modularity_primary",
      "core_class_assortativity",
      "reference_counts_all_match",
      "candidate_KO_selection",
      "bootstrap_design",
      "edge_stability_all",
      "consensus_edges",
      "consensus_nodes",
      "global_metrics",
      "module_summary",
      "module_class_enrichment",
      "edge_class_mixing",
      "module_eigengenes",
      "top_hubs"
    ),
    value = c(
      "PASS",
      RUN156,
      RUN153,
      RUN137,
      "run_accession",
      nrow(Y),
      ncol(Y_universe),
      ncol(Y),
      N_BOOT,
      stats::median(bootstrap_n$n_samples),
      paste(balance_columns, collapse = ";"),
      nrow(stratum_counts),
      nrow(edge_stability),
      nrow(consensus_edges),
      nrow(connected_annotation),
      candidate_graph_density,
      connected_graph_density,
      length(major_modules),
      nrow(major_annotation),
      partition$Q,
      class_assortativity,
      all(reference_check$matches_reference),
      file.path(TABLES_DIR, "171_candidate_KO_selection.tsv"),
      file.path(TABLES_DIR, "171_bootstrap_design.tsv.gz"),
      file.path(TABLES_DIR, "171_edge_stability_all.tsv.gz"),
      file.path(TABLES_DIR, "171_consensus_edges.tsv.gz"),
      file.path(TABLES_DIR, "171_consensus_nodes.tsv"),
      file.path(TABLES_DIR, "171_network_global_metrics.tsv"),
      file.path(TABLES_DIR, "171_module_summary.tsv"),
      file.path(TABLES_DIR, "171_module_core_class_enrichment.tsv"),
      file.path(TABLES_DIR, "171_edge_core_class_mixing.tsv"),
      file.path(TABLES_DIR, "171_major_module_eigengenes.tsv"),
      file.path(TABLES_DIR, "171_top_hubs_annotated.tsv")
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(summary_df, file.path(TABLES_DIR, "171_summary.tsv"))

  saveRDS(
    list(
      graph = graph,
      community = community_fit,
      partition_multistart = partition,
      nodes = nodes,
      consensus_edges = as.data.frame(consensus_edges),
      candidate_annotation = candidate_annotation,
      module_summary = module_summary,
      module_eigengenes = module_eigengenes,
      parameters = parameters
    ),
    file.path(RDS_DIR, "171_consensus_network_objects.rds")
  )

  writeLines("PASS", STATUS_FILE)

  cat("\n============================================================\n")
  cat("STATUS: PASS\n")
  cat("Output dir: ", OUTDIR, "\n", sep = "")
  cat("Samples: ", nrow(Y), "\n", sep = "")
  cat("Candidate KOs: ", ncol(Y), "\n", sep = "")
  cat("Median bootstrap n: ", stats::median(bootstrap_n$n_samples), "\n", sep = "")
  cat("Consensus edges: ", nrow(consensus_edges), "\n", sep = "")
  cat("Connected nodes: ", nrow(connected_annotation), "\n", sep = "")
  cat("Major modules: ", length(major_modules), "\n", sep = "")
  cat("KOs in major modules: ", nrow(major_annotation), "\n", sep = "")
  cat("Modularity: ", signif(partition$Q, 5), "\n", sep = "")
  cat("Reference counts all match: ", all(reference_check$matches_reference), "\n", sep = "")
  cat("Validation:\n")
  print(validation, row.names = FALSE)
  cat("Finished: ", as.character(Sys.time()), "\n", sep = "")
  cat("============================================================\n")
}

status <- tryCatch(
  {
    main()
    0L
  },
  error = function(error_i) {
    message("ERROR: ", conditionMessage(error_i))
    writeLines("FAIL", STATUS_FILE)
    writeLines(
      c(
        paste0("time\t", as.character(Sys.time())),
        paste0("error\t", conditionMessage(error_i))
      ),
      file.path(LOGS_DIR, "171_error.txt")
    )
    1L
  }
)

quit(save = "no", status = status)
