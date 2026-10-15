#!/usr/bin/env Rscript

# =============================================================================
# 174_VALIDATE_FUNCTIONAL_NETWORK_FINAL.R
#
# Statistical validation of the final functional persistence network.
#
# Node classes are inherited unchanged from script 153: absolute global core,
# near-global, realm-core only, multi-realm non-core and realm-specific.
#
# Important scope:
#   * The network selected the 1,500 highest-variance KOs from the 5,000-KO
#     network matrix without using persistence-class labels.
#   * Post-network inferential units are KOs, edges, or graph partitions; the
#     194 metagenomes are not treated as independent replicates here.
#   * Associations are co-abundance patterns, not direct interactions, causal
#     effects, expression, metabolic activity, or biogeochemical fluxes.
# =============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Se requiere el paquete 'data.table'.")
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Se requiere el paquete 'igraph'.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Se requiere el paquete 'ggplot2'.")
  }
})

get_arg <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- which(args == flag)
  if (length(hit) == 0L || hit[length(hit)] == length(args)) return(default)
  args[hit[length(hit)] + 1L]
}

as_flag <- function(x) {
  tolower(trimws(as.character(x))) %in% c("1", "true", "t", "yes", "y")
}

BASE <- normalizePath(
  get_arg("--base", "/home/fjbalvino/Biogeografia_2026"),
  mustWork = FALSE
)
RESULTS_DIR <- file.path(BASE, "resultados")
SCRIPT_BASENAME <- "174_VALIDATE_FUNCTIONAL_NETWORK_FINAL"
STAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUTDIR <- file.path(RESULTS_DIR, paste0(SCRIPT_BASENAME, "_", STAMP))
TABLES_DIR <- file.path(OUTDIR, "tables")
FIGURES_DIR <- file.path(OUTDIR, "figures")
LOGS_DIR <- file.path(OUTDIR, "logs")
RDS_DIR <- file.path(OUTDIR, "rds")

for (d in c(TABLES_DIR, FIGURES_DIR, LOGS_DIR, RDS_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

LATEST_FILE <- file.path(RESULTS_DIR, paste0("LATEST_", SCRIPT_BASENAME, ".txt"))
writeLines(OUTDIR, LATEST_FILE)
writeLines("RUNNING", file.path(OUTDIR, "STATUS.txt"))

read_latest <- function(prefix) {
  path <- file.path(RESULTS_DIR, paste0("LATEST_", prefix, ".txt"))
  if (!file.exists(path)) return(NA_character_)
  x <- readLines(path, warn = FALSE)
  if (length(x) == 0L) return(NA_character_)
  trimws(x[1])
}

RUN156 <- get_arg(
  "--run156",
  read_latest("156_build_KO_CLR_matrix_prevalence_classes")
)
RUN171 <- get_arg(
  "--run171",
  read_latest("171_FUNCTIONAL_KO_CONSENSUS_NETWORK_FINAL")
)
RUN173 <- get_arg(
  "--run173",
  read_latest("173_ANNOTATE_FUNCTIONAL_NETWORK_MODULES_FINAL")
)

N_PERM <- as.integer(get_arg("--n_perm", "4999"))
N_REWIRE <- as.integer(get_arg("--n_rewire", "999"))
N_EDGE_BOOT <- as.integer(get_arg("--n_edge_boot", "500"))
MIN_MODULE_SIZE <- as.integer(get_arg("--min_module_size", "10"))
PROCESS_MIN_TOTAL <- as.integer(get_arg("--process_min_total", "3"))
KEGG_MAPPING_COMPLETE <- as_flag(get_arg("--kegg_mapping_complete", "0"))
BLOCK_TARGET_SIZE <- as.integer(get_arg("--block_target_size", "35"))
BLOCK_MIN_SIZE <- as.integer(get_arg("--block_min_size", "12"))
SEED <- as.integer(get_arg("--seed", "123"))
DPI <- as.integer(get_arg("--dpi", "600"))
STRICT_EXPECTED <- as_flag(get_arg("--strict_expected", "1"))
ENFORCE_MANUSCRIPT_REFERENCE <- as_flag(
  get_arg("--enforce_manuscript_reference", "1")
)

EXPECTED <- list(
  network_universe = 5000L,
  candidates = 1500L,
  consensus_edges = 1223L,
  connected = 326L,
  major_modules = 6L,
  major_nodes = 237L,
  major_absolute_global_core = 9L,
  major_near_global = 195L,
  major_realm_core_only = 33L,
  major_multi_realm_non_core = 0L,
  major_realm_specific = 0L
)

CLASS_LEVELS <- c(
  "absolute_global_core", "near_global", "realm_core_only",
  "multi_realm_non_core", "realm_specific"
)

read_tsv <- function(path, required = TRUE) {
  if (is.na(path) || path == "" || !file.exists(path)) {
    if (required) stop("No existe archivo requerido: ", path)
    return(data.table::data.table())
  }
  data.table::fread(
    path,
    sep = "\t",
    header = TRUE,
    data.table = TRUE,
    check.names = FALSE,
    showProgress = FALSE
  )
}

write_tsv <- function(x, path) {
  data.table::fwrite(x, path, sep = "\t", quote = FALSE, na = "NA")
}

clean_chr <- function(x) {
  y <- trimws(as.character(x))
  y[is.na(y) | y == "" | tolower(y) %in% c(
    "na", "n/a", "nan", "null", "none", "unknown"
  )] <- NA_character_
  y
}

safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

pick_col <- function(nms, candidates, required = FALSE, label = NULL) {
  hit <- candidates[tolower(candidates) %in% tolower(nms)]
  if (length(hit) == 0L) {
    if (required) {
      stop(
        "No se identificó columna ",
        ifelse(is.null(label), paste(candidates, collapse = "/"), label),
        ". Disponibles: ", paste(nms, collapse = "; ")
      )
    }
    return(NA_character_)
  }
  nms[match(tolower(hit[1]), tolower(nms))]
}

require_cols <- function(x, cols, label) {
  missing <- setdiff(cols, names(x))
  if (length(missing) > 0L) {
    stop("Faltan columnas en ", label, ": ", paste(missing, collapse = ", "))
  }
}

normalize_core_class <- function(x) {
  y <- tolower(gsub("[ _-]+", "", clean_chr(x)))
  out <- rep(NA_character_, length(y))
  out[y %in% c("absoluteglobalcore", "absolutecore")] <- "absolute_global_core"
  out[y %in% c("nearglobal", "nearglobalcore")] <- "near_global"
  out[y %in% c("realmcoreonly", "realmcore")] <- "realm_core_only"
  out[y %in% c("multirealmnoncore", "multirealm") ] <- "multi_realm_non_core"
  out[y %in% c("realmspecific", "singlerealm") ] <- "realm_specific"
  out
}

canonical_pair <- function(a, b) {
  paste(pmin(as.character(a), as.character(b)),
        pmax(as.character(a), as.character(b)), sep = "__")
}

empirical_p_upper <- function(observed, null) {
  null <- null[is.finite(null)]
  if (!is.finite(observed) || length(null) == 0L) return(NA_real_)
  (1 + sum(null >= observed)) / (length(null) + 1)
}

empirical_p_lower <- function(observed, null) {
  null <- null[is.finite(null)]
  if (!is.finite(observed) || length(null) == 0L) return(NA_real_)
  (1 + sum(null <= observed)) / (length(null) + 1)
}

empirical_p_two <- function(observed, null) {
  pu <- empirical_p_upper(observed, null)
  pl <- empirical_p_lower(observed, null)
  if (!is.finite(pu) || !is.finite(pl)) return(NA_real_)
  min(1, 2 * min(pu, pl))
}

comb2 <- function(x) x * (x - 1) / 2

adjusted_rand_index <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  x <- as.character(x[ok])
  y <- as.character(y[ok])
  n <- length(x)
  if (n < 2L) return(NA_real_)
  tab <- table(x, y)
  a <- rowSums(tab)
  b <- colSums(tab)
  index <- sum(comb2(tab))
  expected <- sum(comb2(a)) * sum(comb2(b)) / comb2(n)
  maximum <- 0.5 * (sum(comb2(a)) + sum(comb2(b)))
  den <- maximum - expected
  if (!is.finite(den) || den == 0) return(NA_real_)
  (index - expected) / den
}

progress_tick <- function(i, n, label, every = 250L) {
  if (i == 1L || i == n || i %% every == 0L) {
    cat("    ", label, " ", i, "/", n, "\n", sep = "")
    flush.console()
  }
}

save_gg <- function(plot, stem, width, height) {
  ggplot2::ggsave(
    paste0(stem, ".pdf"), plot, width = width, height = height,
    units = "in", device = grDevices::cairo_pdf, bg = "white"
  )
  ggplot2::ggsave(
    paste0(stem, ".png"), plot, width = width, height = height,
    units = "in", dpi = DPI, bg = "white"
  )
  ggplot2::ggsave(
    paste0(stem, ".tiff"), plot, width = width, height = height,
    units = "in", dpi = DPI, compression = "lzw", bg = "white"
  )
}

standardize_KO_table <- function(dt, label) {
  nms <- names(dt)
  ko_col <- pick_col(
    nms,
    c("KO", "ko", "KEGG_ko", "kegg_ko", "orthology", "name", "node"),
    required = TRUE,
    label = paste0("KO en ", label)
  )
  out <- data.table::data.table(KO = clean_chr(dt[[ko_col]]))

  aliases <- list(
    functional_core_class = c(
      "functional_core_class", "core_class", "class"
    ),
    prevalence_n = c(
      "prevalence_n", "global_prevalence_n", "n_detected", "detected_n"
    ),
    prevalence_prop = c(
      "prevalence_prop", "global_prevalence_prop", "prevalence"
    ),
    total_count = c(
      "total_count", "total_abundance", "total_count_script154"
    ),
    mean_relative_abundance = c(
      "mean_relative_abundance", "mean_relative_abundance_global",
      "mean_rel_abund"
    ),
    CLR_variance = c(
      "CLR_variance", "recalculated_CLR_variance", "clr_variance"
    ),
    selection_rank = c(
      "selection_rank_by_recalculated_CLR_variance",
      "network_rank_by_CLR_variance", "selection_rank"
    ),
    connected = c("connected", "is_connected"),
    degree = c("degree", "source_degree", "plot_degree"),
    strength_abs = c(
      "strength_abs", "source_strength_abs", "strength", "plot_strength"
    ),
    module = c("module"),
    module_label = c("module_label"),
    module_display = c("module_display"),
    module_program = c("module_program"),
    major_module = c("major_module", "is_major_module"),
    functional_category = c("functional_category")
  )

  for (nm in names(aliases)) {
    col <- pick_col(nms, aliases[[nm]])
    if (!is.na(col)) out[[nm]] <- dt[[col]]
  }

  out <- out[!is.na(KO) & KO != ""]
  if (anyDuplicated(out$KO)) stop("KOs duplicados en ", label, ".")

  if ("functional_core_class" %in% names(out)) {
    out[, functional_core_class := normalize_core_class(
      functional_core_class
    )]
  }

  for (nm in intersect(
    c(
      "prevalence_n", "prevalence_prop", "total_count",
      "mean_relative_abundance", "CLR_variance", "selection_rank",
      "degree", "strength_abs", "module"
    ),
    names(out)
  )) {
    out[, (nm) := safe_num(get(nm))]
  }

  for (nm in intersect(c("connected", "major_module"), names(out))) {
    if (!is.logical(out[[nm]])) {
      out[, (nm) := tolower(as.character(get(nm))) %in% c("1", "true", "t", "yes")]
    }
  }

  out
}

standardize_edges <- function(dt, label) {
  nms <- names(dt)
  c1 <- pick_col(nms, c("KO1", "from", "node1", "source", "Var1"), TRUE,
                 paste0("primer extremo en ", label))
  c2 <- pick_col(nms, c("KO2", "to", "node2", "target", "Var2"), TRUE,
                 paste0("segundo extremo en ", label))
  out <- data.table::data.table(
    KO1 = clean_chr(dt[[c1]]),
    KO2 = clean_chr(dt[[c2]])
  )

  stability_col <- pick_col(nms, c(
    "edge_stability", "stability", "bootstrap_stability", "selection_frequency"
  ))
  sign_col <- pick_col(nms, c(
    "sign_consistency", "sign_stability", "sign_agreement"
  ))
  weight_col <- pick_col(nms, c(
    "mean_abs_weight", "mean_abs_prop", "mean_abs_association",
    "abs_weight", "weight"
  ))

  out[, edge_stability := if (is.na(stability_col)) NA_real_ else safe_num(dt[[stability_col]])]
  out[, sign_consistency := if (is.na(sign_col)) NA_real_ else safe_num(dt[[sign_col]])]
  out[, mean_abs_weight := if (is.na(weight_col)) 1 else abs(safe_num(dt[[weight_col]]))]
  out <- out[!is.na(KO1) & !is.na(KO2) & KO1 != KO2]
  out[, pair_id := canonical_pair(KO1, KO2)]
  data.table::setorder(out, pair_id, -edge_stability, -mean_abs_weight)
  out <- out[!duplicated(pair_id)]
  out
}

merge_fill <- function(base, add) {
  add <- data.table::copy(add)
  add <- add[!duplicated(KO)]
  common <- intersect(setdiff(names(add), "KO"), names(base))
  if (length(common) > 0L) {
    data.table::setnames(add, common, paste0(common, "__add"))
  }
  out <- merge(base, add, by = "KO", all.x = TRUE, sort = FALSE)
  for (nm in common) {
    add_nm <- paste0(nm, "__add")
    old <- out[[nm]]
    new <- out[[add_nm]]
    if (is.numeric(old) || is.integer(old)) {
      use <- !is.finite(safe_num(old)) & is.finite(safe_num(new))
    } else if (is.logical(old)) {
      use <- is.na(old) & !is.na(new)
    } else {
      use <- is.na(clean_chr(old)) & !is.na(clean_chr(new))
    }
    out[[nm]][use] <- new[use]
    out[[add_nm]] <- NULL
  }
  out
}

robust_scale_vector <- function(x) {
  x <- safe_num(x)
  med <- stats::median(x, na.rm = TRUE)
  if (!is.finite(med)) med <- 0
  x[!is.finite(x)] <- med
  s <- stats::mad(x, center = med, constant = 1, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) return(rep(0, length(x)))
  (x - med) / s
}

make_covariate_blocks <- function(dt, vars, target_size, min_size, seed) {
  vars <- vars[vars %in% names(dt)]
  if (length(vars) == 0L || nrow(dt) < 2L * min_size) {
    return(rep(1L, nrow(dt)))
  }
  X <- do.call(cbind, lapply(vars, function(v) robust_scale_vector(dt[[v]])))
  X <- as.matrix(X)
  keep <- apply(X, 2, function(z) stats::sd(z) > 0)
  X <- X[, keep, drop = FALSE]
  if (ncol(X) == 0L) return(rep(1L, nrow(dt)))

  k <- max(1L, floor(nrow(dt) / target_size))
  k <- min(k, floor(nrow(dt) / min_size))
  if (k <= 1L) return(rep(1L, nrow(dt)))

  set.seed(seed)
  km <- stats::kmeans(X, centers = k, nstart = 50, iter.max = 200)
  cl <- km$cluster

  repeat {
    sizes <- table(cl)
    small <- as.integer(names(sizes)[sizes < min_size])
    if (length(small) == 0L || length(unique(cl)) <= 1L) break
    s <- small[which.min(sizes[as.character(small)])]
    idx_s <- which(cl == s)
    other <- setdiff(unique(cl), s)
    if (length(other) == 0L) break
    center_s <- colMeans(X[idx_s, , drop = FALSE])
    distances <- vapply(other, function(o) {
      center_o <- colMeans(X[cl == o, , drop = FALSE])
      sqrt(sum((center_s - center_o)^2))
    }, numeric(1))
    cl[idx_s] <- other[which.min(distances)]
  }

  as.integer(match(cl, sort(unique(cl))))
}

permute_within_blocks <- function(x, blocks) {
  out <- x
  for (b in unique(blocks)) {
    idx <- which(blocks == b)
    if (length(idx) > 1L) out[idx] <- sample(x[idx], length(idx), replace = FALSE)
  }
  out
}

block_audit <- function(dt, blocks, analysis, covariates) {
  ans <- data.table::data.table(block = blocks)[, .(n_KOs = .N), by = block]
  ans[, analysis := analysis]
  ans[, covariates := paste(covariates, collapse = ";")]
  data.table::setcolorder(ans, c("analysis", "covariates", "block", "n_KOs"))
  ans
}

selection_test <- function(dt, selected, analysis, covariates, n_perm, seed) {
  selected <- as.logical(selected)
  if (length(selected) != nrow(dt) || anyNA(selected)) {
    stop("Indicador de selección inválido en: ", analysis)
  }
  blocks <- make_covariate_blocks(
    dt, covariates, BLOCK_TARGET_SIZE, BLOCK_MIN_SIZE, seed
  )
  classes <- CLASS_LEVELS
  observed <- vapply(classes, function(z) {
    mean(dt$functional_core_class[selected] == z)
  }, numeric(1))
  null <- matrix(NA_real_, nrow = n_perm, ncol = length(classes),
                 dimnames = list(NULL, classes))
  set.seed(seed + 1L)
  for (p in seq_len(n_perm)) {
    selected_p <- permute_within_blocks(selected, blocks)
    null[p, ] <- vapply(classes, function(z) {
      mean(dt$functional_core_class[selected_p] == z)
    }, numeric(1))
    progress_tick(p, n_perm, paste0(analysis, " permutations"))
  }
  result <- data.table::data.table(
    analysis = analysis,
    functional_core_class = classes,
    observed_selected_proportion = observed,
    null_expected_proportion = colMeans(null, na.rm = TRUE),
    null_sd = apply(null, 2, stats::sd, na.rm = TRUE)
  )
  result[, difference := observed_selected_proportion - null_expected_proportion]
  result[, fold_enrichment := (observed_selected_proportion + 1e-8) /
           (null_expected_proportion + 1e-8)]
  result[, p_upper := vapply(seq_len(.N), function(i) {
    empirical_p_upper(observed_selected_proportion[i], null[, i])
  }, numeric(1))]
  result[, p_two_sided := vapply(seq_len(.N), function(i) {
    empirical_p_two(observed_selected_proportion[i], null[, i])
  }, numeric(1))]
  result[, q_BY := stats::p.adjust(p_upper, method = "BY")]
  result[, selected_n := sum(selected)]
  result[, background_n := nrow(dt)]
  result[, n_blocks := data.table::uniqueN(blocks)]
  result[, covariates := paste(covariates, collapse = ";")]
  list(
    result = result,
    blocks = block_audit(dt, blocks, analysis, covariates)
  )
}

enrichment_from_null <- function(observed_matrix, null_matrix,
                                 module_levels, category_levels,
                                 category_name) {
  observed <- as.vector(observed_matrix)
  null_mean <- colMeans(null_matrix, na.rm = TRUE)
  null_sd <- apply(null_matrix, 2, stats::sd, na.rm = TRUE)
  z_observed <- (observed - null_mean) / null_sd
  z_observed[!is.finite(z_observed)] <- NA_real_
  z_null <- sweep(sweep(null_matrix, 2, null_mean, "-"), 2, null_sd, "/")
  z_null[!is.finite(z_null)] <- NA_real_
  max_z <- apply(z_null, 1, function(z) {
    if (all(is.na(z))) NA_real_ else max(z, na.rm = TRUE)
  })
  grid <- data.table::as.data.table(expand.grid(
    module_label = module_levels,
    category = category_levels,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ))
  grid[, observed_n := observed]
  grid[, expected_n := null_mean]
  grid[, null_sd := null_sd]
  grid[, fold_enrichment := (observed_n + 0.5) / (expected_n + 0.5)]
  grid[, log2_fold_enrichment := log2(fold_enrichment)]
  grid[, z_permutation := z_observed]
  grid[, p_empirical := vapply(seq_along(observed), function(i) {
    empirical_p_upper(observed[i], null_matrix[, i])
  }, numeric(1))]
  grid[, q_BY := stats::p.adjust(p_empirical, method = "BY")]
  grid[, p_maxT_FWER := vapply(z_observed, function(z) {
    empirical_p_upper(z, max_z)
  }, numeric(1))]
  grid[, enrichment_type := category_name]
  grid
}

run_label_enrichment <- function(module_idx, category_idx,
                                 module_levels, category_levels,
                                 perm_module_matrix, category_name) {
  observed <- table(
    factor(module_idx, levels = seq_along(module_levels)),
    factor(category_idx, levels = seq_along(category_levels))
  )
  null <- matrix(
    0,
    nrow = ncol(perm_module_matrix),
    ncol = length(module_levels) * length(category_levels)
  )
  for (p in seq_len(ncol(perm_module_matrix))) {
    tab <- table(
      factor(perm_module_matrix[, p], levels = seq_along(module_levels)),
      factor(category_idx, levels = seq_along(category_levels))
    )
    null[p, ] <- as.vector(tab)
    progress_tick(p, ncol(perm_module_matrix),
                  paste0(category_name, " enrichment"))
  }
  enrichment_from_null(
    observed, null, module_levels, category_levels, category_name
  )
}

build_graph <- function(edges, all_vertices = NULL) {
  ed <- data.table::copy(edges)
  ed <- ed[!is.na(KO1) & !is.na(KO2) & KO1 != KO2]
  if (nrow(ed) == 0L) {
    g <- igraph::make_empty_graph(
      n = ifelse(is.null(all_vertices), 0L, length(all_vertices)),
      directed = FALSE
    )
    if (!is.null(all_vertices)) igraph::V(g)$name <- all_vertices
    return(g)
  }
  data.table::setorder(ed, pair_id, -mean_abs_weight)
  ed <- ed[!duplicated(pair_id)]
  edge_df <- data.frame(
    from = ed$KO1,
    to = ed$KO2,
    weight = pmax(ed$mean_abs_weight, 1e-8),
    stringsAsFactors = FALSE
  )
  vertices <- if (is.null(all_vertices)) {
    sort(unique(c(ed$KO1, ed$KO2)))
  } else {
    unique(all_vertices)
  }
  g <- igraph::graph_from_data_frame(
    edge_df,
    directed = FALSE,
    vertices = data.frame(name = vertices, stringsAsFactors = FALSE)
  )
  igraph::simplify(
    g, remove.multiple = TRUE, remove.loops = TRUE,
    edge.attr.comb = list(weight = "max")
  )
}

community_partition <- function(g, method = "louvain") {
  # Match inference in 10C: optimize on the complete graph, not separately
  # within components (which changes the modularity normalization).
  vertex_names <- igraph::V(g)$name
  if (igraph::ecount(g) == 0L) {
    return(setNames(seq_along(vertex_names), vertex_names))
  }
  weights <- igraph::E(g)$weight
  cl <- tryCatch({
    if (method == "louvain") {
      igraph::cluster_louvain(g, weights = weights)
    } else if (method == "infomap") {
      igraph::cluster_infomap(g, e.weights = weights)
    } else if (method == "walktrap") {
      igraph::cluster_walktrap(g, weights = weights)
    } else if (method == "leiden") {
      igraph::cluster_leiden(g, objective_function = "modularity",
                            weights = weights, resolution = 1, n_iterations = 10)
    } else {
      stop("Algoritmo no reconocido: ", method)
    }
  }, error = function(e) {
    stop("Community detection failed for ", method, ": ", conditionMessage(e))
  })
  membership <- igraph::membership(cl)
  stopifnot(setequal(names(membership), vertex_names))
  membership[vertex_names]
}

graph_metrics <- function(g, membership) {
  degree <- igraph::degree(g)
  connected <- names(degree)[degree > 0]
  if (length(connected) == 0L) {
    return(list(
      n_nodes = igraph::vcount(g), n_connected_nodes = 0L,
      n_edges = igraph::ecount(g), largest_component = 0L,
      modularity = NA_real_, n_modules = NA_integer_, n_large_modules = 0L
    ))
  }
  comp <- igraph::components(g)
  modularity <- tryCatch(
    igraph::modularity(
      g, membership = membership[igraph::V(g)$name],
      weights = igraph::E(g)$weight
    ),
    error = function(e) NA_real_
  )
  sizes <- table(membership[connected])
  list(
    n_nodes = igraph::vcount(g),
    n_connected_nodes = length(connected),
    n_edges = igraph::ecount(g),
    largest_component = max(comp$csize),
    modularity = modularity,
    n_modules = length(sizes),
    n_large_modules = sum(sizes >= MIN_MODULE_SIZE)
  )
}

membership_for_reference <- function(membership, reference_kos) {
  x <- membership[reference_kos]
  missing <- is.na(x)
  if (any(missing)) {
    start <- ifelse(all(is.na(x)), 0L, max(x, na.rm = TRUE))
    x[missing] <- start + seq_len(sum(missing))
  }
  as.character(x)
}

main <- function() {
  log_file <- file.path(LOGS_DIR, "174_log.txt")
  log_con <- file(log_file, open = "wt")
  sink(log_con, split = TRUE)
  sink(log_con, type = "message")
  on.exit({
    # Para type="message", sink.number() devuelve el número de conexión.
    # El valor 2 corresponde al stderr normal; no es un contador de sinks.
    if (sink.number(type = "message") != 2L) sink(type = "message")
    while (sink.number(type = "output") > 0L) sink(type = "output")
    try(close(log_con), silent = TRUE)
  }, add = TRUE)

  cat("============================================================\n")
  cat("174_VALIDATE_FUNCTIONAL_NETWORK_FINAL\n")
  cat("Started: ", as.character(Sys.time()), "\n", sep = "")
  cat("RUN156: ", RUN156, "\n", sep = "")
  cat("RUN171: ", RUN171, "\n", sep = "")
  cat("RUN173: ", RUN173, "\n", sep = "")
  cat("N_PERM: ", N_PERM, "\n", sep = "")
  cat("N_REWIRE: ", N_REWIRE, "\n", sep = "")
  cat("N_EDGE_BOOT: ", N_EDGE_BOOT, "\n", sep = "")
  cat("STRICT_EXPECTED: ", STRICT_EXPECTED, "\n", sep = "")
  cat("ENFORCE_MANUSCRIPT_REFERENCE: ",
      ENFORCE_MANUSCRIPT_REFERENCE, "\n", sep = "")
  cat("OUTDIR: ", OUTDIR, "\n", sep = "")
  cat("============================================================\n\n")

  for (run in c(RUN156, RUN171, RUN173)) {
    if (is.na(run) || run == "" || !dir.exists(run)) stop("Run inválido: ", run)
  }
  if (N_PERM < 99L || N_REWIRE < 99L || N_EDGE_BOOT < 0L) {
    stop("Los números de permutaciones/rewires/perturbaciones son inválidos.")
  }

  paths <- list(
    universe = file.path(RUN171, "tables", "171_network_universe_KO_annotation.tsv"),
    candidates = file.path(RUN171, "tables", "171_candidate_KO_selection.tsv"),
    nodes163 = file.path(RUN171, "tables", "171_consensus_nodes.tsv"),
    edges163 = file.path(RUN171, "tables", "171_consensus_edges.tsv.gz"),
    edge_all = file.path(RUN171, "tables", "171_edge_stability_all.tsv.gz"),
    nodes165 = file.path(RUN173, "tables", "173_network_nodes_annotated.tsv"),
    kegg_map = file.path(RUN173, "tables", "173_KO_to_KEGG_MODULE_snapshot_long.tsv")
  )

  cat("[1] Leyendo y alineando universos, nodos, aristas y anotaciones...\n")
  universe <- standardize_KO_table(read_tsv(paths$universe), "universo 156")
  candidate_selection <- standardize_KO_table(
    read_tsv(paths$candidates), "candidatos 163"
  )
  candidates <- standardize_KO_table(read_tsv(paths$nodes163), "nodos 163")
  candidates <- merge_fill(candidates, candidate_selection)
  major <- standardize_KO_table(read_tsv(paths$nodes165), "nodos 165")
  edges <- standardize_edges(read_tsv(paths$edges163), "aristas consenso 163")
  edge_all <- standardize_edges(read_tsv(paths$edge_all), "aristas candidatas 163")
  kegg_map <- read_tsv(paths$kegg_map, required = FALSE)

  for (nm in c("universe", "candidates", "major")) {
    dt <- get(nm)
    require_cols(dt, "functional_core_class", nm)
    if (anyNA(dt$functional_core_class)) {
      stop("Clases de prevalencia faltantes o inválidas en ", nm, ".")
    }
    if (any(!dt$functional_core_class %in% CLASS_LEVELS)) {
      stop("Clase funcional inesperada en ", nm, ".")
    }
    assign(nm, dt)
  }

  if (!setequal(candidates$KO, candidate_selection$KO)) {
    stop("Los KOs de 171_consensus_nodes y 171_candidate_KO_selection no coinciden.")
  }
  if (!all(candidates$KO %in% universe$KO)) {
    stop("Hay candidatos 163 ausentes del universo 156.")
  }
  if (!all(major$KO %in% candidates$KO)) {
    stop("Hay nodos 165 ausentes de los candidatos 163.")
  }

  connected_kos <- sort(unique(c(edges$KO1, edges$KO2)))
  candidates[, connected := KO %in% connected_kos]
  candidates[, major_module := KO %in% major$KO]

  degree_dt <- data.table::rbindlist(list(
    edges[, .(KO = KO1, degree_add = 1, strength_add = mean_abs_weight)],
    edges[, .(KO = KO2, degree_add = 1, strength_add = mean_abs_weight)]
  ))[, .(
    degree_calc = sum(degree_add),
    strength_calc = sum(strength_add, na.rm = TRUE)
  ), by = KO]
  candidates <- merge(candidates, degree_dt, by = "KO", all.x = TRUE, sort = FALSE)
  if (!"degree" %in% names(candidates)) candidates[, degree := NA_real_]
  if (!"strength_abs" %in% names(candidates)) candidates[, strength_abs := NA_real_]
  candidates[!is.finite(degree), degree := degree_calc]
  candidates[!is.finite(strength_abs), strength_abs := strength_calc]
  candidates[!is.finite(degree), degree := 0]
  candidates[!is.finite(strength_abs), strength_abs := 0]
  candidates[, c("degree_calc", "strength_calc") := NULL]

  major <- merge_fill(major, candidates)
  if (!"module_label" %in% names(major)) {
    if (!"module" %in% names(major)) stop("No hay módulo en nodos 165.")
    major[, module_label := paste0("M", sprintf("%02d", as.integer(module)))]
  }
  major[, module_label := clean_chr(module_label)]
  if (anyNA(major$module_label)) stop("Hay module_label faltantes en nodos 165.")
  module_levels <- sort(unique(major$module_label))

  for (dt_name in c("universe", "candidates", "major")) {
    dt <- data.table::copy(get(dt_name))
    if (!"prevalence_n" %in% names(dt)) dt[, prevalence_n := NA_real_]
    if (!"prevalence_prop" %in% names(dt)) dt[, prevalence_prop := NA_real_]
    if (!"total_count" %in% names(dt)) dt[, total_count := NA_real_]
    if (!"mean_relative_abundance" %in% names(dt)) {
      dt[, mean_relative_abundance := NA_real_]
    }
    if (!"CLR_variance" %in% names(dt)) dt[, CLR_variance := NA_real_]
    dt[, log10_total_count := log10(pmax(total_count, 0) + 1)]
    dt[, log10_mean_relative_abundance := log10(
      pmax(mean_relative_abundance, 0) + 1e-12
    )]
    assign(dt_name, dt)
  }

  cat("    Network-universe KOs: ", nrow(universe), "\n", sep = "")
  cat("    Candidate KOs: ", nrow(candidates), "\n", sep = "")
  cat("    Consensus edges: ", nrow(edges), "\n", sep = "")
  cat("    Connected KOs: ", length(connected_kos), "\n", sep = "")
  cat("    Major-module KOs: ", nrow(major), "\n", sep = "")
  cat("    Modules: ", paste(module_levels, collapse = "; "), "\n", sep = "")

  scope <- c(
    "STATISTICAL SCOPE AND ASSUMPTIONS",
    "",
    "1. Functional persistence classes are inherited unchanged from final script 153.",
    "2. Realm-core support requires all localities in at least one realm with >=3 localities and mean realm-locality relative abundance >=0.0001.",
    "3. Candidate selection was deterministic by CLR variance. Its class composition is a methodological selection audit, not independent biological evidence.",
    "5. Connected-node enrichment is tested by blocked permutations controlling prevalence, total abundance and CLR variance.",
    "6. Major-module and module-specific enrichment additionally control graph degree and absolute strength.",
    "7. KOs, edges and graph partitions—not metagenomes—are the inferential units in post-network tests.",
    "8. Class assortativity is evaluated with degree-preserving graph rewiring.",
    "9. Robustness is evaluated across edge thresholds, community algorithms and edge-stability perturbations, reclustering the full candidate graph as in primary inference.",
    "10. Functional and KEGG process enrichment is conditional/post-selection evidence because annotation and module discovery use the same network.",
    "11. KEGG categories can overlap; process results report BY FDR and permutation maxT family-wise correction.",
    "12. Co-abundance does not demonstrate direct interaction, causality, expression, activity, or biogeochemical flux."
  )
  writeLines(scope, file.path(TABLES_DIR, "174_statistical_scope_and_assumptions.txt"))

  cat("[2] Resumiendo clases a través de universos anidados...\n")
  summarize_classes <- function(dt, universe_name) {
    dt[, .(n_KOs = .N), by = functional_core_class][
      , `:=`(
        universe = universe_name,
        total_KOs = nrow(dt),
        proportion = n_KOs / nrow(dt)
      )
    ]
  }
  nested <- data.table::rbindlist(list(
    summarize_classes(universe, "Network universe (5,000)"),
    summarize_classes(candidates, "Selected candidates (1,500)"),
    summarize_classes(candidates[connected == TRUE], "Consensus-connected"),
    summarize_classes(major, "Major modules")
  ), use.names = TRUE, fill = TRUE)
  nested[, functional_core_class := factor(
    functional_core_class, levels = CLASS_LEVELS
  )]
  nested[, universe := factor(
    universe,
    levels = c(
      "Network universe (5,000)", "Selected candidates (1,500)",
      "Consensus-connected", "Major modules"
    )
  )]
  data.table::setorder(nested, universe, functional_core_class)
  write_tsv(nested, file.path(TABLES_DIR, "174_core_classes_nested_universes.tsv"))

  # Unconditioned selection audit only; do not interpret as biological validation.
  candidate_audit <- data.table::rbindlist(lapply(CLASS_LEVELS, function(cl) {
    selected <- candidates$functional_core_class == cl
    background <- universe$functional_core_class == cl
    a <- sum(selected)
    b <- nrow(candidates) - a
    c0 <- sum(background) - a
    d <- nrow(universe) - nrow(candidates) - c0
    ft <- stats::fisher.test(matrix(c(a, b, c0, d), nrow = 2, byrow = TRUE))
    data.table::data.table(
      functional_core_class = cl,
      selected_n = a,
      selected_total = nrow(candidates),
      universe_n = sum(background),
      universe_total = nrow(universe),
      selected_proportion = a / nrow(candidates),
      universe_proportion = sum(background) / nrow(universe),
      odds_ratio_unconditioned = unname(ft$estimate),
      p_fisher_unconditioned = ft$p.value,
      inferential_scope = "Methodological selection audit; candidate selection is deterministic by CLR variance"
    )
  }))
  candidate_audit[, q_BY_unconditioned := stats::p.adjust(
    p_fisher_unconditioned, method = "BY"
  )]
  write_tsv(candidate_audit, file.path(TABLES_DIR, "174_candidate_selection_method_audit.tsv"))

  cat("[3] Probando enriquecimiento condicional del backbone...\n")
  connected_test <- selection_test(
    candidates,
    candidates$connected,
    "Connected versus 1,500 candidates",
    c("prevalence_prop", "log10_total_count", "CLR_variance"),
    N_PERM,
    SEED + 100L
  )
  major_connected_test <- selection_test(
    candidates[connected == TRUE],
    candidates[connected == TRUE]$major_module,
    "Major modules versus connected nodes",
    c(
      "prevalence_prop", "log10_total_count", "CLR_variance",
      "degree", "strength_abs"
    ),
    N_PERM,
    SEED + 200L
  )
  major_candidate_test <- selection_test(
    candidates,
    candidates$major_module,
    "Major modules versus 1,500 candidates",
    c("prevalence_prop", "log10_total_count", "CLR_variance"),
    N_PERM,
    SEED + 300L
  )
  selection_results <- data.table::rbindlist(list(
    connected_test$result,
    major_connected_test$result,
    major_candidate_test$result
  ))
  selection_blocks <- data.table::rbindlist(list(
    connected_test$blocks,
    major_connected_test$blocks,
    major_candidate_test$blocks
  ))
  write_tsv(selection_results, file.path(
    TABLES_DIR, "174_backbone_selection_matched_permutation.tsv"
  ))
  write_tsv(selection_blocks, file.path(
    TABLES_DIR, "174_covariate_permutation_blocks.tsv"
  ))

  cat("[4] Construyendo permutaciones comunes de módulos...\n")
  major[, module_idx := match(module_label, module_levels)]
  module_covariates <- c(
    "prevalence_prop", "log10_total_count", "CLR_variance",
    "degree", "strength_abs"
  )
  module_blocks <- make_covariate_blocks(
    major, module_covariates, BLOCK_TARGET_SIZE, BLOCK_MIN_SIZE, SEED + 400L
  )
  perm_module_matrix <- matrix(
    NA_integer_, nrow = nrow(major), ncol = N_PERM
  )
  set.seed(SEED + 401L)
  for (p in seq_len(N_PERM)) {
    perm_module_matrix[, p] <- permute_within_blocks(major$module_idx, module_blocks)
    progress_tick(p, N_PERM, "module-label permutations")
  }
  write_tsv(
    block_audit(major, module_blocks, "Module enrichment", module_covariates),
    file.path(TABLES_DIR, "174_module_permutation_blocks.tsv")
  )

  cat("[5] Probando enriquecimiento de clases funcionales y categorías...\n")
  network_class_levels <- CLASS_LEVELS
  major[, class_idx := match(functional_core_class, network_class_levels)]
  if (anyNA(major$class_idx)) {
    stop("Los módulos mayores contienen clases funcionales no canónicas.")
  }
  prevalence_enrichment <- run_label_enrichment(
    major$module_idx,
    major$class_idx,
    module_levels,
    network_class_levels,
    perm_module_matrix,
    "functional_core_class"
  )

  functional_enrichment <- data.table::data.table()
  if ("functional_category" %in% names(major)) {
    major[, functional_category := clean_chr(functional_category)]
    major[is.na(functional_category), functional_category := "Unannotated"]
    functional_levels <- sort(unique(major$functional_category))
    major[, functional_idx := match(functional_category, functional_levels)]
    functional_enrichment <- run_label_enrichment(
      major$module_idx,
      major$functional_idx,
      module_levels,
      functional_levels,
      perm_module_matrix,
      "functional_category"
    )
  }
  write_tsv(prevalence_enrichment, file.path(
    TABLES_DIR, "174_module_core_class_enrichment_permutation.tsv"
  ))
  write_tsv(functional_enrichment, file.path(
    TABLES_DIR, "174_module_functional_category_enrichment_permutation.tsv"
  ))
  supported_functional <- functional_enrichment[
    is.finite(q_BY) & q_BY <= 0.05 &
      is.finite(p_maxT_FWER) & p_maxT_FWER <= 0.05
  ]
  write_tsv(supported_functional, file.path(
    TABLES_DIR, "174_supported_module_functional_specialization.tsv"
  ))

  expected_supported <- data.table::data.table(
    module_label = c("M02", "M03", "M06"),
    category = c(
      "Information processing",
      "Anaerobic redox / energy",
      "Transport / membrane systems"
    )
  )
  observed_supported_keys <- paste(
    supported_functional$module_label,
    supported_functional$category,
    sep = "||"
  )
  expected_supported_keys <- paste(
    expected_supported$module_label,
    expected_supported$category,
    sep = "||"
  )
  specialization_reference <- data.table::data.table(
    check = c(
      "exact_supported_module_category_pairs",
      "fold_enrichment_within_reported_range_2.22_3.36",
      "BY_q_values_compatible_with_reported_0.0061"
    ),
    pass = c(
      setequal(observed_supported_keys, expected_supported_keys),
      nrow(supported_functional) == 3L &&
        all(supported_functional$fold_enrichment >= 2.20 &
              supported_functional$fold_enrichment <= 3.40),
      nrow(supported_functional) == 3L &&
        all(abs(supported_functional$q_BY - 0.0061) <= 0.001)
    )
  )
  write_tsv(specialization_reference, file.path(
    TABLES_DIR, "174_manuscript_specialization_reference_check.tsv"
  ))
  if (ENFORCE_MANUSCRIPT_REFERENCE && any(!specialization_reference$pass)) {
    stop(
      "La especializacion funcional no reprodujo exactamente la referencia ",
      "del manuscrito; revisar 174_manuscript_specialization_reference_check.tsv."
    )
  }

  cat("[6] Probando enriquecimiento y solapamiento de procesos KEGG...\n")
  process_enrichment <- data.table::data.table()
  process_overlap <- data.table::data.table()
  writeLines(if (KEGG_MAPPING_COMPLETE) {
    "KEGG process tests enabled: user asserted complete mapping coverage for the analysis universe."
  } else {
    "KEGG process tests omitted: the supplied snapshot is restricted to historical network nodes; missing links cannot be treated as confirmed absence."
  }, file.path(TABLES_DIR, "174_KEGG_mapping_scope.txt"))
  if (KEGG_MAPPING_COMPLETE && nrow(kegg_map) > 0L) {
    ko_col <- pick_col(names(kegg_map), c("KO", "ko"), TRUE, "KO en mapa KEGG")
    process_col <- pick_col(
      names(kegg_map),
      c("KEGG_module_name", "module_name", "KEGG_module_id"),
      TRUE,
      "proceso KEGG"
    )
    process_map <- unique(data.table::data.table(
      KO = clean_chr(kegg_map[[ko_col]]),
      process = clean_chr(kegg_map[[process_col]])
    ))
    process_map <- process_map[
      !is.na(KO) & !is.na(process) & KO %in% major$KO
    ]
    process_totals <- process_map[, .(total_KOs = data.table::uniqueN(KO)), by = process]
    keep_process <- process_totals[total_KOs >= PROCESS_MIN_TOTAL]$process
    process_map <- process_map[process %in% keep_process]

    if (nrow(process_map) > 0L) {
      process_levels <- sort(unique(process_map$process))
      A <- matrix(
        0L, nrow = nrow(major), ncol = length(process_levels),
        dimnames = list(major$KO, process_levels)
      )
      A[cbind(match(process_map$KO, major$KO),
              match(process_map$process, process_levels))] <- 1L
      observed_process_counts <- rowsum(
        A,
        group = factor(major$module_idx, levels = seq_along(module_levels)),
        reorder = FALSE
      )
      null_process_counts <- matrix(
        0,
        nrow = N_PERM,
        ncol = length(module_levels) * length(process_levels)
      )
      module_pairs <- t(utils::combn(seq_along(module_levels), 2))
      null_overlap <- matrix(NA_real_, nrow = N_PERM, ncol = nrow(module_pairs))

      for (p in seq_len(N_PERM)) {
        count_p <- rowsum(
          A,
          group = factor(
            perm_module_matrix[, p], levels = seq_along(module_levels)
          ),
          reorder = FALSE
        )
        null_process_counts[p, ] <- as.vector(count_p)
        present_p <- count_p > 0
        for (q in seq_len(nrow(module_pairs))) {
          i <- module_pairs[q, 1]
          j <- module_pairs[q, 2]
          inter <- sum(present_p[i, ] & present_p[j, ])
          union_n <- sum(present_p[i, ] | present_p[j, ])
          null_overlap[p, q] <- ifelse(union_n == 0L, NA_real_, inter / union_n)
        }
        progress_tick(p, N_PERM, "KEGG-process permutations")
      }

      process_enrichment <- enrichment_from_null(
        observed_process_counts,
        null_process_counts,
        module_levels,
        process_levels,
        "KEGG_process"
      )
      process_enrichment <- merge(
        process_enrichment,
        process_totals,
        by.x = "category",
        by.y = "process",
        all.x = TRUE,
        sort = FALSE
      )

      observed_present <- observed_process_counts > 0
      overlap_rows <- vector("list", nrow(module_pairs))
      for (q in seq_len(nrow(module_pairs))) {
        i <- module_pairs[q, 1]
        j <- module_pairs[q, 2]
        shared <- process_levels[observed_present[i, ] & observed_present[j, ]]
        union_set <- process_levels[observed_present[i, ] | observed_present[j, ]]
        observed_jaccard <- ifelse(
          length(union_set) == 0L, NA_real_, length(shared) / length(union_set)
        )
        overlap_rows[[q]] <- data.table::data.table(
          module_i = module_levels[i],
          module_j = module_levels[j],
          shared_n = length(shared),
          union_n = length(union_set),
          observed_jaccard = observed_jaccard,
          null_mean_jaccard = mean(null_overlap[, q], na.rm = TRUE),
          null_sd_jaccard = stats::sd(null_overlap[, q], na.rm = TRUE),
          p_upper = empirical_p_upper(observed_jaccard, null_overlap[, q]),
          p_two_sided = empirical_p_two(observed_jaccard, null_overlap[, q]),
          shared_processes = ifelse(
            length(shared) == 0L, NA_character_, paste(shared, collapse = " | ")
          )
        )
      }
      process_overlap <- data.table::rbindlist(overlap_rows)
      process_overlap[, q_BY := stats::p.adjust(p_upper, method = "BY")]
    }
  }
  write_tsv(process_enrichment, file.path(
    TABLES_DIR, "174_module_KEGG_process_enrichment_permutation.tsv"
  ))
  write_tsv(process_overlap, file.path(
    TABLES_DIR, "174_module_process_overlap_significance.tsv"
  ))

  cat("[7] Evaluando mezcla de clases con rewiring degree-preserving...\n")
  g_primary <- build_graph(edges, all_vertices = candidates$KO)
  vertex_classes <- candidates[
    match(igraph::V(g_primary)$name, KO)
  ]$functional_core_class
  if (anyNA(vertex_classes) || any(!vertex_classes %in% network_class_levels)) {
    stop("La red conectada contiene clases funcionales no canónicas.")
  }
  igraph::V(g_primary)$core_class <- vertex_classes

  class_assortativity <- function(g) {
    types <- as.integer(factor(
      igraph::V(g)$core_class,
      levels = network_class_levels
    ))
    # igraph >=2 requires non-negative nominal types; 1/2 are valid.
    tryCatch(
      igraph::assortativity_nominal(g, types = types, directed = FALSE),
      error = function(e) NA_real_
    )
  }
  near_global_edge_fraction <- function(g) {
    endpoint_names <- igraph::ends(g, igraph::E(g), names = TRUE)
    c1 <- igraph::V(g)$core_class[
      match(endpoint_names[, 1], igraph::V(g)$name)
    ]
    c2 <- igraph::V(g)$core_class[
      match(endpoint_names[, 2], igraph::V(g)$name)
    ]
    mean(c1 == "near_global" & c2 == "near_global")
  }

  observed_assortativity <- class_assortativity(g_primary)
  observed_near_global_edges <- near_global_edge_fraction(g_primary)
  null_assortativity <- rep(NA_real_, N_REWIRE)
  null_near_global_edges <- rep(NA_real_, N_REWIRE)
  set.seed(SEED + 500L)
  for (b in seq_len(N_REWIRE)) {
    g_null <- tryCatch(
      igraph::rewire(
        g_primary,
        with = igraph::keeping_degseq(
          niter = max(100L, 10L * igraph::ecount(g_primary)),
          loops = FALSE
        )
      ),
      error = function(e) NULL
    )
    if (!is.null(g_null)) {
      null_assortativity[b] <- class_assortativity(g_null)
      null_near_global_edges[b] <- near_global_edge_fraction(g_null)
    }
    progress_tick(b, N_REWIRE, "degree-preserving rewires", every = 100L)
  }
  assortativity_results <- data.table::data.table(
    statistic = c(
      "Nominal assortativity among five functional persistence classes",
      "Fraction of edges joining two near-global KOs"
    ),
    observed = c(observed_assortativity, observed_near_global_edges),
    null_mean = c(
      mean(null_assortativity, na.rm = TRUE),
      mean(null_near_global_edges, na.rm = TRUE)
    ),
    null_sd = c(
      stats::sd(null_assortativity, na.rm = TRUE),
      stats::sd(null_near_global_edges, na.rm = TRUE)
    ),
    p_upper = c(
      empirical_p_upper(observed_assortativity, null_assortativity),
      empirical_p_upper(observed_near_global_edges, null_near_global_edges)
    ),
    p_two_sided = c(
      empirical_p_two(observed_assortativity, null_assortativity),
      empirical_p_two(observed_near_global_edges, null_near_global_edges)
    ),
    n_valid_rewires = c(
      sum(is.finite(null_assortativity)),
      sum(is.finite(null_near_global_edges))
    )
  )
  endpoint_names <- igraph::ends(g_primary, igraph::E(g_primary), names = TRUE)
  edge_mixing <- data.table::data.table(
    class1 = igraph::V(g_primary)$core_class[
      match(endpoint_names[, 1], igraph::V(g_primary)$name)
    ],
    class2 = igraph::V(g_primary)$core_class[
      match(endpoint_names[, 2], igraph::V(g_primary)$name)
    ]
  )
  edge_mixing[, edge_type := ifelse(
    class1 <= class2,
    paste(class1, class2, sep = "__"),
    paste(class2, class1, sep = "__")
  )]
  edge_mixing_summary <- edge_mixing[, .(n_edges = .N), by = edge_type]
  edge_mixing_summary[, proportion := n_edges / sum(n_edges)]
  write_tsv(assortativity_results, file.path(
    TABLES_DIR, "174_prevalence_degree_preserving_assortativity.tsv"
  ))
  write_tsv(edge_mixing_summary, file.path(
    TABLES_DIR, "174_edge_core_class_mixing.tsv"
  ))

  cat("[8] Evaluando sensibilidad a umbrales...\n")
  threshold_sensitivity <- data.table::data.table()
  if (nrow(edge_all) > 0L && any(is.finite(edge_all$edge_stability)) &&
      any(is.finite(edge_all$mean_abs_weight))) {
    grid <- data.table::CJ(
      stability_min = c(0.70, 0.75, 0.80, 0.85, 0.90),
      mean_abs_min = c(0.90, 0.95, 0.97),
      sign_consistency_min = c(0.95, 1.00)
    )
    primary_edge_set <- unique(edges$pair_id)
    rows <- vector("list", nrow(grid))
    for (i in seq_len(nrow(grid))) {
      threshold <- grid[i]
      ed <- edge_all[
        edge_stability >= threshold$stability_min &
          mean_abs_weight >= threshold$mean_abs_min &
          (is.na(sign_consistency) |
             sign_consistency >= threshold$sign_consistency_min)
      ]
      g <- build_graph(ed, all_vertices = candidates$KO)
      membership <- community_partition(g, "louvain")
      metrics <- graph_metrics(g, membership)
      predicted_reference <- membership_for_reference(membership, major$KO)
      edge_set <- unique(ed$pair_id)
      rows[[i]] <- data.table::data.table(
        stability_min = threshold$stability_min,
        mean_abs_min = threshold$mean_abs_min,
        sign_consistency_min = threshold$sign_consistency_min,
        n_edges = metrics$n_edges,
        n_connected_nodes = metrics$n_connected_nodes,
        largest_component = metrics$largest_component,
        modularity = metrics$modularity,
        n_modules = metrics$n_modules,
        n_large_modules = metrics$n_large_modules,
        edge_jaccard_vs_primary = ifelse(
          length(union(edge_set, primary_edge_set)) == 0L,
          NA_real_,
          length(intersect(edge_set, primary_edge_set)) /
            length(union(edge_set, primary_edge_set))
        ),
        ARI_vs_primary_major_modules = adjusted_rand_index(
          major$module_label, predicted_reference
        ),
        major_node_retention = mean(
          major$KO %in% unique(c(ed$KO1, ed$KO2))
        )
      )
      progress_tick(i, nrow(grid), "threshold combinations", every = 5L)
    }
    threshold_sensitivity <- data.table::rbindlist(rows)
  }
  write_tsv(threshold_sensitivity, file.path(
    TABLES_DIR, "174_network_threshold_sensitivity.tsv"
  ))

  cat("[9] Comparando algoritmos comunitarios...\n")
  methods <- c("louvain", "infomap", "walktrap")
  if ("cluster_leiden" %in% getNamespaceExports("igraph")) {
    methods <- c(methods, "leiden")
  }
  algorithm_memberships <- list()
  algorithm_rows <- vector("list", length(methods))
  for (i in seq_along(methods)) {
    method <- methods[i]
    membership <- community_partition(g_primary, method)
    algorithm_memberships[[method]] <- membership
    metrics <- graph_metrics(g_primary, membership)
    predicted_reference <- membership_for_reference(membership, major$KO)
    algorithm_rows[[i]] <- data.table::data.table(
      method = method,
      ARI_vs_primary_major_modules = adjusted_rand_index(
        major$module_label, predicted_reference
      ),
      modularity = metrics$modularity,
      n_modules = metrics$n_modules,
      n_large_modules = metrics$n_large_modules,
      largest_component = metrics$largest_component
    )
  }
  algorithm_sensitivity <- data.table::rbindlist(algorithm_rows)
  algorithm_pairwise <- data.table::data.table()
  if (length(methods) >= 2L) {
    pairs <- t(utils::combn(methods, 2))
    algorithm_pairwise <- data.table::rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
      data.table::data.table(
        method1 = pairs[i, 1],
        method2 = pairs[i, 2],
        ARI_pairwise = adjusted_rand_index(
          algorithm_memberships[[pairs[i, 1]]],
          algorithm_memberships[[pairs[i, 2]]]
        )
      )
    }))
  }
  write_tsv(algorithm_sensitivity, file.path(
    TABLES_DIR, "174_community_algorithm_sensitivity.tsv"
  ))
  write_tsv(algorithm_pairwise, file.path(
    TABLES_DIR, "174_community_algorithm_pairwise_ARI.tsv"
  ))

  cat("[10] Perturbando aristas según su estabilidad bootstrap...\n")
  edge_boot_summary <- data.table::data.table()
  node_stability <- data.table::data.table()
  module_stability <- data.table::data.table()
  if (N_EDGE_BOOT > 0L && nrow(edges) > 0L) {
    edges_boot <- data.table::copy(edges)
    edges_boot[!is.finite(edge_stability), edge_stability := 1]
    edges_boot[, edge_stability := pmin(1, pmax(0, edge_stability))]
    node_accumulator <- setNames(rep(0, nrow(major)), major$KO)
    boot_rows <- vector("list", N_EDGE_BOOT)
    set.seed(SEED + 600L)

    for (b in seq_len(N_EDGE_BOOT)) {
      keep <- stats::runif(nrow(edges_boot)) <= edges_boot$edge_stability
      ed_b <- edges_boot[keep]
      g_b <- build_graph(ed_b, all_vertices = candidates$KO)
      membership_b <- community_partition(g_b, "louvain")
      predicted_reference <- membership_for_reference(membership_b, major$KO)
      ari <- adjusted_rand_index(major$module_label, predicted_reference)

      for (module_i in module_levels) {
        kos <- major[module_label == module_i]$KO
        if (length(kos) <= 1L) {
          node_accumulator[kos] <- node_accumulator[kos] + 1
          next
        }
        memberships <- predicted_reference[match(kos, major$KO)]
        cluster_sizes <- table(memberships)
        scores <- (as.numeric(cluster_sizes[memberships]) - 1) / (length(kos) - 1)
        node_accumulator[kos] <- node_accumulator[kos] + scores
      }
      boot_rows[[b]] <- data.table::data.table(
        perturbation = b,
        retained_edges = nrow(ed_b),
        retained_edge_fraction = nrow(ed_b) / nrow(edges_boot),
        ARI_vs_primary_modules = ari
      )
      progress_tick(b, N_EDGE_BOOT, "edge perturbations", every = 50L)
    }
    edge_boot_summary <- data.table::rbindlist(boot_rows)
    node_stability <- data.table::data.table(
      KO = names(node_accumulator),
      mean_within_module_coassignment = as.numeric(node_accumulator / N_EDGE_BOOT)
    )
    node_stability <- merge(
      node_stability,
      major[, .(KO, module_label, functional_core_class)],
      by = "KO",
      all.x = TRUE,
      sort = FALSE
    )
    module_stability <- node_stability[, .(
      n_KOs = .N,
      mean_stability = mean(mean_within_module_coassignment, na.rm = TRUE),
      median_stability = stats::median(
        mean_within_module_coassignment, na.rm = TRUE
      ),
      min_stability = min(mean_within_module_coassignment, na.rm = TRUE),
      q10_stability = as.numeric(stats::quantile(
        mean_within_module_coassignment, 0.10, na.rm = TRUE
      ))
    ), by = module_label]
  }
  write_tsv(edge_boot_summary, file.path(
    TABLES_DIR, "174_edge_perturbation_summary.tsv"
  ))
  write_tsv(node_stability, file.path(
    TABLES_DIR, "174_node_module_membership_stability.tsv"
  ))
  write_tsv(module_stability, file.path(
    TABLES_DIR, "174_module_membership_stability.tsv"
  ))

  cat("[11] Generando figuras diagnósticas...\n")
  class_colors <- c(
    absolute_global_core = "#7F0000", near_global = "#C62828",
    realm_core_only = "#EF5350", multi_realm_non_core = "#6B6B6B",
    realm_specific = "#BDBDBD"
  )
  nested_plot <- data.table::copy(nested)
  nested_plot[, proportion_percent := 100 * proportion]
  p_nested <- ggplot2::ggplot(
    nested_plot,
    ggplot2::aes(
      x = universe,
      y = proportion_percent,
      fill = functional_core_class
    )
  ) +
    ggplot2::geom_col(width = 0.72) +
    ggplot2::scale_fill_manual(values = class_colors, drop = FALSE) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = NULL,
      y = "KOs (%)",
      fill = "Functional persistence class",
      title = "Functional persistence architecture across nested network universes"
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
  save_gg(p_nested, file.path(FIGURES_DIR, "174_nested_functional_core_composition"), 8.5, 4.8)

  selection_plot_data <- selection_results[
    functional_core_class == "near_global"
  ]
  p_selection <- ggplot2::ggplot(
    selection_plot_data,
    ggplot2::aes(
      x = reorder(analysis, fold_enrichment),
      y = fold_enrichment
    )
  ) +
    ggplot2::geom_hline(yintercept = 1, color = "grey65", linetype = 2) +
    ggplot2::geom_point(size = 3.5, color = class_colors[["near_global"]]) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = NULL,
      y = "Conditional fold enrichment",
      title = "Conditional enrichment of near-global KOs",
      subtitle = "Blocked permutations control detectability, CLR variance and topology as applicable"
    ) +
    ggplot2::theme_classic(base_size = 11)
  save_gg(p_selection, file.path(FIGURES_DIR, "174_conditional_near_global_enrichment"), 8.5, 4.8)

  module_plot_data <- prevalence_enrichment[
    category == "near_global"
  ]
  module_plot_data[, significant := p_maxT_FWER <= 0.05]
  p_module <- ggplot2::ggplot(
    module_plot_data,
    ggplot2::aes(
      x = module_label,
      y = log2_fold_enrichment,
      fill = significant
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "grey65", linetype = 2) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#BDBDBD")) +
    ggplot2::labs(
      x = "Network module",
      y = "log2 conditional fold enrichment",
      fill = "maxT P <= 0.05",
      title = "Module-specific enrichment of near-global KOs"
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
  save_gg(p_module, file.path(FIGURES_DIR, "174_module_near_global_enrichment"), 8.5, 4.8)

  if (nrow(threshold_sensitivity) > 0L) {
    threshold_plot_data <- threshold_sensitivity[
      sign_consistency_min == 1
    ]
    p_threshold <- ggplot2::ggplot(
      threshold_plot_data,
      ggplot2::aes(
        x = stability_min,
        y = ARI_vs_primary_major_modules,
        color = factor(mean_abs_min),
        group = factor(mean_abs_min)
      )
    ) +
      ggplot2::geom_line(linewidth = 0.7) +
      ggplot2::geom_point(size = 2.4) +
      ggplot2::scale_color_brewer(palette = "Dark2") +
      ggplot2::labs(
        x = "Minimum bootstrap edge stability",
        y = "ARI versus primary modules",
        color = "Minimum |association|",
        title = "Network-partition sensitivity to consensus thresholds"
      ) +
      ggplot2::theme_classic(base_size = 11) +
      ggplot2::theme(legend.position = "bottom")
    save_gg(p_threshold, file.path(FIGURES_DIR, "174_threshold_partition_sensitivity"), 8.5, 4.8)
  }

  cat("[12] Validaciones finales y resumen...\n")
  class_major <- table(factor(
    major$functional_core_class,
    levels = CLASS_LEVELS
  ))
  validations <- data.table::data.table(
    metric = c(
      "n_network_universe_KOs",
      "n_candidate_KOs",
      "n_consensus_edges",
      "n_connected_KOs",
      "n_major_modules",
      "n_major_module_KOs",
      "n_major_absolute_global_core",
      "n_major_near_global",
      "n_major_realm_core_only",
      "n_major_multi_realm_non_core",
      "n_major_realm_specific",
      "n_missing_candidate_classes",
      "n_missing_major_classes",
      "n_missing_edge_endpoints",
      "n_duplicate_consensus_edges",
      "n_nonfinite_consensus_weights",
      "n_prevalence_enrichment_rows",
      "n_valid_assortativity_rewires",
      "n_threshold_grid_rows"
    ),
    observed = c(
      nrow(universe),
      nrow(candidates),
      nrow(edges),
      length(connected_kos),
      length(module_levels),
      nrow(major),
      unname(class_major["absolute_global_core"]),
      unname(class_major["near_global"]),
      unname(class_major["realm_core_only"]),
      unname(class_major["multi_realm_non_core"]),
      unname(class_major["realm_specific"]),
      sum(is.na(candidates$functional_core_class)),
      sum(is.na(major$functional_core_class)),
      sum(!edges$KO1 %in% candidates$KO | !edges$KO2 %in% candidates$KO),
      anyDuplicated(edges$pair_id),
      sum(!is.finite(edges$mean_abs_weight)),
      nrow(prevalence_enrichment),
      min(assortativity_results$n_valid_rewires),
      nrow(threshold_sensitivity)
    ),
    expected = c(
      EXPECTED$network_universe,
      EXPECTED$candidates,
      EXPECTED$consensus_edges,
      EXPECTED$connected,
      EXPECTED$major_modules,
      EXPECTED$major_nodes,
      EXPECTED$major_absolute_global_core,
      EXPECTED$major_near_global,
      EXPECTED$major_realm_core_only,
      EXPECTED$major_multi_realm_non_core,
      EXPECTED$major_realm_specific,
      0L,
      0L,
      0L,
      0L,
      0L,
      length(module_levels) * length(CLASS_LEVELS),
      N_REWIRE,
      30L
    )
  )
  validations[, pass := observed == expected]
  write_tsv(validations, file.path(TABLES_DIR, "174_expected_validation.tsv"))

  if (STRICT_EXPECTED && any(!validations$pass)) {
    stop(
      "Falló al menos una validación estricta. Revisa 174_expected_validation.tsv"
    )
  }

  summary_table <- data.table::data.table(
    item = c(
      "status", "run156", "run171", "run173", "n_perm", "n_rewire",
      "n_edge_boot", "network_universe_KOs", "candidate_KOs",
      "connected_KOs", "major_module_KOs", "major_modules",
      "consensus_edges", "major_absolute_global_core", "major_near_global",
      "major_realm_core_only", "major_multi_realm_non_core",
      "major_realm_specific",
      "network_scope", "output_dir"
    ),
    value = c(
      "PASS", RUN156, RUN171, RUN173, N_PERM, N_REWIRE, N_EDGE_BOOT,
      nrow(universe), nrow(candidates), length(connected_kos), nrow(major),
      length(module_levels), nrow(edges),
      unname(class_major["absolute_global_core"]),
      unname(class_major["near_global"]),
      unname(class_major["realm_core_only"]),
      unname(class_major["multi_realm_non_core"]),
      unname(class_major["realm_specific"]),
      "Five final persistence classes inherited from script 153",
      OUTDIR
    )
  )
  write_tsv(summary_table, file.path(TABLES_DIR, "174_summary.tsv"))

  saveRDS(list(
    nested = nested,
    candidate_selection_audit = candidate_audit,
    selection_results = selection_results,
    prevalence_enrichment = prevalence_enrichment,
    functional_enrichment = functional_enrichment,
    supported_functional_specialization = supported_functional,
    manuscript_specialization_reference = specialization_reference,
    process_enrichment = process_enrichment,
    process_overlap = process_overlap,
    assortativity = assortativity_results,
    threshold_sensitivity = threshold_sensitivity,
    algorithm_sensitivity = algorithm_sensitivity,
    edge_perturbation = edge_boot_summary,
    node_stability = node_stability,
    module_stability = module_stability
  ), file.path(RDS_DIR, "174_validation_results.rds"))

  writeLines("PASS", file.path(OUTDIR, "STATUS.txt"))
  cat("\n============================================================\n")
  cat("STATUS: PASS\n")
  cat("Output dir: ", OUTDIR, "\n", sep = "")
  cat("Network-universe KOs: ", nrow(universe), "\n", sep = "")
  cat("Candidate KOs: ", nrow(candidates), "\n", sep = "")
  cat("Connected KOs: ", length(connected_kos), "\n", sep = "")
  cat("Major modules: ", length(module_levels), "\n", sep = "")
  cat("Major-module KOs: ", nrow(major), "\n", sep = "")
  cat("Major absolute/near-global/realm-core/multi-realm/realm-specific: ", paste(as.integer(class_major), collapse = "/"), "\n", sep = "")
  cat("Permutations: ", N_PERM, "\n", sep = "")
  cat("Degree-preserving rewires: ", N_REWIRE, "\n", sep = "")
  cat("Edge perturbations: ", N_EDGE_BOOT, "\n", sep = "")
  cat("Finished: ", as.character(Sys.time()), "\n", sep = "")
  cat("============================================================\n")
}

status <- tryCatch(
  {
    main()
    "PASS"
  },
  error = function(e) {
    writeLines("FAIL", file.path(OUTDIR, "STATUS.txt"))
    writeLines(conditionMessage(e), file.path(LOGS_DIR, "174_error.txt"))
    message("ERROR: ", conditionMessage(e))
    "FAIL"
  }
)

if (status != "PASS") quit(save = "no", status = 1L)
