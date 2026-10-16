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
# 172_CURATE_FUNCTIONAL_NETWORK_LAYOUT_FINAL.R
#
# Objetivo
#   Curar y visualizar la red consenso validada por el script 163 sin volver a
#   inferir aristas ni comunidades. Se conservan solamente los modulos con al
#   menos 10 KOs para la interpretacion biologica y se fija un layout
#   reproducible para los analisis de anotacion y validacion posteriores.
#
# Arquitectura funcional
#   Cinco clases jerarquicas heredadas de 153: absolute global core,
#   near-global, realm-core only, multi-realm non-core y realm-specific.
#
# Alcance
#   Las aristas son asociaciones reproducibles de coabundancia en espacio CLR.
#   No representan interacciones bioquimicas directas, causalidad ni actividad.
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

RUN171 <- get_arg(
  "--run171",
  read_latest("171_FUNCTIONAL_KO_CONSENSUS_NETWORK_FINAL")
)

MIN_MODULE_SIZE <- as.integer(get_arg("--min_module_size", "10"))
TOP_LABELS_PER_MODULE <- as.integer(get_arg("--top_labels_per_module", "3"))
SEED <- as.integer(get_arg("--seed", "123"))
NITER_LAYOUT <- as.integer(get_arg("--niter_layout", "2000"))
PLOT_WIDTH <- as.numeric(get_arg("--plot_width", "10"))
PLOT_HEIGHT <- as.numeric(get_arg("--plot_height", "8"))
DPI <- as.integer(get_arg("--dpi", "600"))
STRICT_EXPECTED <- as_flag(get_arg("--strict_expected", "0"))

EXPECTED_SOURCE_EDGES <- as.integer(
  get_arg("--expected_source_edges", "1223")
)
EXPECTED_CANDIDATE_KOS <- as.integer(
  get_arg("--expected_candidate_kos", "1500")
)
EXPECTED_CONNECTED_KOS <- as.integer(
  get_arg("--expected_connected_kos", "326")
)
EXPECTED_MAJOR_MODULES <- as.integer(
  get_arg("--expected_major_modules", "6")
)
EXPECTED_MAJOR_MODULE_KOS <- as.integer(
  get_arg("--expected_major_module_kos", "237")
)
SCRIPT_NUMBER <- "172"
SCRIPT_BASENAME <- "CURATE_FUNCTIONAL_NETWORK_LAYOUT_FINAL"
STAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUTDIR <- file.path(
  RESULTS_DIR,
  paste0(SCRIPT_NUMBER, "_", SCRIPT_BASENAME, "_", STAMP)
)
TABLES_DIR <- file.path(OUTDIR, "tables")
FIGURES_DIR <- file.path(OUTDIR, "figures")
LOGS_DIR <- file.path(OUTDIR, "logs")
RDS_DIR <- file.path(OUTDIR, "rds")
LATEST_FILE <- file.path(
  RESULTS_DIR,
  paste0("LATEST_", SCRIPT_NUMBER, "_", SCRIPT_BASENAME, ".txt")
)
STATUS_FILE <- file.path(OUTDIR, "STATUS.txt")

CLASS_LEVELS <- c(
  "absolute_global_core", "near_global", "realm_core_only",
  "multi_realm_non_core", "realm_specific"
)
CLASS_COLORS <- c(
  "absolute_global_core" = "#7F0000",
  "near_global" = "#C62828",
  "realm_core_only" = "#EF5350",
  "multi_realm_non_core" = "#6B6B6B",
  "realm_specific" = "#BDBDBD"
)
CLASS_LABELS <- c(
  "absolute_global_core" = "Absolute global core",
  "near_global" = "Near-global",
  "realm_core_only" = "Realm-core only",
  "multi_realm_non_core" = "Multi-realm non-core",
  "realm_specific" = "Realm-specific"
)

dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
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

as_logical_safe <- function(x) {
  if (is.logical(x)) return(x)
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

read_tsv <- function(path, as_data_table = FALSE) {
  if (is.na(path) || path == "" || !file.exists(path)) {
    stop("No existe archivo requerido: ", path)
  }
  data.table::fread(
    path,
    sep = "\t",
    header = TRUE,
    data.table = as_data_table,
    check.names = FALSE,
    showProgress = FALSE
  )
}

write_tsv <- function(x, path) {
  data.table::fwrite(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )
}

require_columns <- function(x, columns, label) {
  missing <- setdiff(columns, names(x))
  if (length(missing) > 0L) {
    stop("Faltan columnas en ", label, ": ", paste(missing, collapse = ", "))
  }
}

validate_run <- function(run_path, label) {
  if (is.na(run_path) || run_path == "" || !dir.exists(run_path)) {
    stop(label, " invalido: ", run_path)
  }
  status_path <- file.path(run_path, "STATUS.txt")
  if (!file.exists(status_path)) stop(label, " no contiene STATUS.txt.")
  status <- trimws(readLines(status_path, warn = FALSE)[1L])
  if (!identical(status, "PASS")) stop(label, " no tiene STATUS PASS.")
}

shorten_label <- function(x, max_chars = 58L) {
  y <- clean_chr(x)
  y[is.na(y)] <- ""
  y <- gsub("\\s*\\[EC:[^\\]]+\\]", "", y)
  y <- gsub("\\s+", " ", y)
  y <- trimws(y)
  ifelse(
    nchar(y) > max_chars,
    paste0(substr(y, 1L, max_chars - 3L), "..."),
    y
  )
}

save_base_plot <- function(file_base, width, height, plot_fun, dpi) {
  grDevices::cairo_pdf(
    filename = paste0(file_base, ".pdf"),
    width = width,
    height = height
  )
  plot_fun()
  grDevices::dev.off()

  grDevices::png(
    filename = paste0(file_base, ".png"),
    width = round(width * dpi),
    height = round(height * dpi),
    res = dpi
  )
  plot_fun()
  grDevices::dev.off()

  grDevices::tiff(
    filename = paste0(file_base, ".tiff"),
    width = round(width * dpi),
    height = round(height * dpi),
    res = dpi,
    compression = "lzw"
  )
  plot_fun()
  grDevices::dev.off()
}

main <- function() {
  validate_run(RUN171, "RUN171")
  if (
    MIN_MODULE_SIZE < 2L || TOP_LABELS_PER_MODULE < 0L ||
      NITER_LAYOUT < 100L || DPI < 72L ||
      !is.finite(PLOT_WIDTH) || !is.finite(PLOT_HEIGHT) ||
      PLOT_WIDTH <= 0 || PLOT_HEIGHT <= 0
  ) {
    stop("Parametros de visualizacion invalidos.")
  }

  log_path <- file.path(LOGS_DIR, "172_log.txt")
  log_connection <- file(log_path, open = "wt")
  sink(log_connection, split = TRUE)
  sink(log_connection, type = "message")
  on.exit({
    # Message sinks report a connection number; 2 is the normal stderr.
    if (sink.number(type = "message") != 2L) sink(type = "message")
    while (sink.number(type = "output") > 0L) sink(type = "output")
    close(log_connection)
  }, add = TRUE)

  cat("============================================================\n")
  cat("172_CURATE_FUNCTIONAL_NETWORK_LAYOUT_FINAL\n")
  cat("Started: ", as.character(Sys.time()), "\n", sep = "")
  cat("RUN171: ", RUN171, "\n", sep = "")
  cat("MIN_MODULE_SIZE: ", MIN_MODULE_SIZE, "\n", sep = "")
  cat("TOP_LABELS_PER_MODULE: ", TOP_LABELS_PER_MODULE, "\n", sep = "")
  cat("SEED: ", SEED, "\n", sep = "")
  cat("NITER_LAYOUT: ", NITER_LAYOUT, "\n", sep = "")
  cat("DPI: ", DPI, "\n", sep = "")
  cat("STRICT_EXPECTED: ", STRICT_EXPECTED, "\n", sep = "")
  cat("OUTDIR: ", OUTDIR, "\n", sep = "")
  cat("============================================================\n\n")

  edges_path <- file.path(RUN171, "tables", "171_consensus_edges.tsv.gz")
  nodes_path <- file.path(RUN171, "tables", "171_consensus_nodes.tsv")
  modules_path <- file.path(RUN171, "tables", "171_module_summary.tsv")
  metrics_path <- file.path(RUN171, "tables", "171_network_global_metrics.tsv")

  cat("[1] Leyendo y validando outputs canonicos del script 163...\n")
  edges <- read_tsv(edges_path, as_data_table = TRUE)
  nodes <- read_tsv(nodes_path, as_data_table = TRUE)
  modules <- read_tsv(modules_path, as_data_table = TRUE)
  global_metrics_163 <- read_tsv(metrics_path, as_data_table = TRUE)

  require_columns(
    edges,
    c(
      "KO1", "KO2", "consensus_weight", "mean_abs_weight",
      "edge_stability", "sign_consistency"
    ),
    "171_consensus_edges.tsv.gz"
  )
  require_columns(
    nodes,
    c(
      "KO", "degree", "strength_abs", "module", "connected",
      "major_module", "functional_core_class", "KO_function"
    ),
    "171_consensus_nodes.tsv"
  )
  require_columns(
    modules,
    c("module", "n_KOs", "major_module", "n_edges_internal"),
    "171_module_summary.tsv"
  )
  require_columns(global_metrics_163, c("metric", "value"), "163 metrics")

  raw_edge_n <- nrow(edges)
  raw_node_n <- nrow(nodes)

  cat("    source edges: ", raw_edge_n, "\n", sep = "")
  cat("    source candidate nodes: ", raw_node_n, "\n", sep = "")

  cat("[2] Curando modulos, nodos y aristas...\n")
  edges[, KO1 := clean_chr(KO1)]
  edges[, KO2 := clean_chr(KO2)]
  edges[, consensus_weight := safe_num(consensus_weight)]
  edges[, mean_abs_weight := safe_num(mean_abs_weight)]
  edges[, edge_stability := safe_num(edge_stability)]
  edges[, sign_consistency := safe_num(sign_consistency)]

  nodes[, KO := clean_chr(KO)]
  nodes[, degree := safe_num(degree)]
  nodes[, strength_abs := safe_num(strength_abs)]
  nodes[, module := as.integer(safe_num(module))]
  nodes[, connected := as_logical_safe(connected)]
  nodes[, major_module := as_logical_safe(major_module)]
  nodes[, functional_core_class := clean_chr(functional_core_class)]
  nodes[, KO_function := clean_chr(KO_function)]
  nodes[is.na(KO_function), KO_function := KO]

  modules[, module := as.integer(safe_num(module))]
  modules[, n_KOs := as.integer(safe_num(n_KOs))]
  modules[, n_edges_internal := as.integer(safe_num(n_edges_internal))]
  modules[, major_module := as_logical_safe(major_module)]

  data.table::setorder(nodes, KO, -degree, -strength_abs)
  nodes <- nodes[!is.na(KO) & KO != "" & !duplicated(KO)]

  curated_modules <- modules[
    major_module & n_KOs >= MIN_MODULE_SIZE & n_edges_internal > 0
  ]
  data.table::setorder(curated_modules, -n_KOs, -n_edges_internal, module)
  curated_modules[, module_rank := seq_len(.N)]
  curated_modules[, module_label := sprintf("M%02d", module_rank)]
  curated_modules[, module_display := paste0(module_label, " (n=", n_KOs, ")")]

  if (nrow(curated_modules) == 0L) {
    stop("No se recuperaron modulos interpretables.")
  }

  nodes_plot <- nodes[
    connected & module %in% curated_modules$module
  ]
  nodes_plot <- merge(
    nodes_plot,
    curated_modules[, .(
      module, module_rank, module_label, module_display
    )],
    by = "module",
    all.x = TRUE,
    sort = FALSE
  )

  edges <- edges[
    !is.na(KO1) & !is.na(KO2) & KO1 != "" & KO2 != "" & KO1 != KO2 &
      is.finite(consensus_weight) & is.finite(mean_abs_weight) &
      is.finite(edge_stability) & is.finite(sign_consistency)
  ]
  edges[, node_a := pmin(KO1, KO2)]
  edges[, node_b := pmax(KO1, KO2)]
  edges[, edge_key := paste(node_a, node_b, sep = "||")]
  edges <- edges[
    ,
    .(
      KO1 = node_a[1L],
      KO2 = node_b[1L],
      consensus_weight = mean(consensus_weight),
      mean_abs_weight = mean(mean_abs_weight),
      edge_stability = max(edge_stability),
      sign_consistency = min(sign_consistency),
      n_source_rows = .N
    ),
    by = edge_key
  ]

  edges_plot <- edges[
    KO1 %in% nodes_plot$KO & KO2 %in% nodes_plot$KO
  ]
  nodes_in_graph <- unique(c(edges_plot$KO1, edges_plot$KO2))
  nodes_dropped_no_major_edge <- setdiff(nodes_plot$KO, nodes_in_graph)
  nodes_plot <- nodes_plot[KO %in% nodes_in_graph]

  module_lookup <- setNames(nodes_plot$module, nodes_plot$KO)
  module_label_lookup <- setNames(nodes_plot$module_label, nodes_plot$KO)
  edges_plot[, module_KO1 := module_lookup[KO1]]
  edges_plot[, module_KO2 := module_lookup[KO2]]
  edges_plot[, module_label_KO1 := module_label_lookup[KO1]]
  edges_plot[, module_label_KO2 := module_label_lookup[KO2]]
  edges_plot[, is_internal_module_edge := module_KO1 == module_KO2]

  cat("    curated modules: ", nrow(curated_modules), "\n", sep = "")
  cat("    curated nodes: ", nrow(nodes_plot), "\n", sep = "")
  cat("    curated edges: ", nrow(edges_plot), "\n", sep = "")
  cat(
    "    nodes dropped without major-module edges: ",
    length(nodes_dropped_no_major_edge), "\n",
    sep = ""
  )

  if (nrow(edges_plot) == 0L || nrow(nodes_plot) == 0L) {
    stop("La curacion produjo una red vacia.")
  }

  cat("[3] Construyendo el grafo de visualizacion...\n")
  edge_matrix <- as.matrix(edges_plot[, .(KO1, KO2)])
  graph <- igraph::graph_from_edgelist(edge_matrix, directed = FALSE)

  igraph::E(graph)$consensus_weight <- edges_plot$consensus_weight
  igraph::E(graph)$mean_abs_weight <- edges_plot$mean_abs_weight
  igraph::E(graph)$edge_stability <- edges_plot$edge_stability
  igraph::E(graph)$sign_consistency <- edges_plot$sign_consistency
  igraph::E(graph)$is_internal_module_edge <-
    edges_plot$is_internal_module_edge

  vertex_names <- igraph::V(graph)$name
  node_index <- match(vertex_names, nodes_plot$KO)
  if (anyNA(node_index)) stop("Faltan nodos anotados en el grafo.")

  igraph::V(graph)$KO <- vertex_names
  igraph::V(graph)$source_module <- nodes_plot$module[node_index]
  igraph::V(graph)$module_rank <- nodes_plot$module_rank[node_index]
  igraph::V(graph)$module_label <- nodes_plot$module_label[node_index]
  igraph::V(graph)$module_display <- nodes_plot$module_display[node_index]
  igraph::V(graph)$functional_core_class <-
    nodes_plot$functional_core_class[node_index]
  igraph::V(graph)$KO_function <- nodes_plot$KO_function[node_index]
  igraph::V(graph)$source_degree <- nodes_plot$degree[node_index]
  igraph::V(graph)$source_strength_abs <- nodes_plot$strength_abs[node_index]

  plot_degree <- igraph::degree(graph)
  plot_strength <- igraph::strength(
    graph,
    weights = igraph::E(graph)$mean_abs_weight
  )
  igraph::V(graph)$plot_degree <- plot_degree
  igraph::V(graph)$plot_strength <- plot_strength

  cat("    graph nodes: ", igraph::vcount(graph), "\n", sep = "")
  cat("    graph edges: ", igraph::ecount(graph), "\n", sep = "")
  cat(
    "    graph density: ",
    signif(igraph::edge_density(graph, loops = FALSE), 6), "\n",
    sep = ""
  )

  cat("[4] Calculando layout reproducible y estetica...\n")
  module_levels <- curated_modules$module_display
  module_levels <- module_levels[
    module_levels %in% unique(igraph::V(graph)$module_display)
  ]
  module_colors <- grDevices::hcl.colors(
    length(module_levels),
    palette = "Dark 3"
  )
  names(module_colors) <- module_levels

  observed_classes <- CLASS_LEVELS[
    CLASS_LEVELS %in% unique(igraph::V(graph)$functional_core_class)
  ]
  observed_class_colors <- CLASS_COLORS[observed_classes]

  module_vertex_colors <- module_colors[
    as.character(igraph::V(graph)$module_display)
  ]
  module_vertex_colors[is.na(module_vertex_colors)] <- "#BDBDBD"
  class_vertex_colors <- CLASS_COLORS[
    as.character(igraph::V(graph)$functional_core_class)
  ]
  class_vertex_colors[is.na(class_vertex_colors)] <- "#BDBDBD"

  strength_values <- igraph::V(graph)$plot_strength
  strength_values[!is.finite(strength_values)] <- 0
  if (max(strength_values) > 0) {
    vertex_sizes <- 3 + 8 * sqrt(strength_values / max(strength_values))
  } else {
    vertex_sizes <- rep(4, length(strength_values))
  }

  nodes_for_labels <- data.table::data.table(
    KO = igraph::V(graph)$KO,
    module_display = igraph::V(graph)$module_display,
    degree = igraph::V(graph)$source_degree,
    strength_abs = igraph::V(graph)$source_strength_abs,
    functional_core_class =
      igraph::V(graph)$functional_core_class,
    KO_function = igraph::V(graph)$KO_function
  )
  data.table::setorder(
    nodes_for_labels,
    module_display,
    -strength_abs,
    -degree,
    KO
  )
  top_labels <- nodes_for_labels[
    ,
    head(.SD, TOP_LABELS_PER_MODULE),
    by = module_display
  ]
  top_labels[, label := paste0(KO, " - ", shorten_label(KO_function))]

  edge_values <- igraph::E(graph)$mean_abs_weight
  if (max(edge_values) > min(edge_values)) {
    edge_widths <- 0.25 + 1.7 * (
      (edge_values - min(edge_values)) /
        (max(edge_values) - min(edge_values))
    )
  } else {
    edge_widths <- rep(0.6, length(edge_values))
  }
  edge_colors <- ifelse(
    igraph::E(graph)$is_internal_module_edge,
    grDevices::adjustcolor("#4D4D4D", alpha.f = 0.28),
    grDevices::adjustcolor("#999999", alpha.f = 0.16)
  )

  igraph::V(graph)$plot_color_module <- module_vertex_colors
  igraph::V(graph)$plot_color_class <- class_vertex_colors
  igraph::V(graph)$plot_size <- vertex_sizes
  igraph::V(graph)$plot_label <- ifelse(
    igraph::V(graph)$KO %in% top_labels$KO,
    igraph::V(graph)$KO,
    NA_character_
  )
  igraph::E(graph)$plot_width <- edge_widths
  igraph::E(graph)$plot_color <- edge_colors

  set.seed(SEED)
  layout_xy <- igraph::layout_with_fr(
    graph,
    weights = igraph::E(graph)$mean_abs_weight,
    niter = NITER_LAYOUT
  )

  layout_table <- data.table::data.table(
    KO = igraph::V(graph)$KO,
    x = layout_xy[, 1L],
    y = layout_xy[, 2L],
    source_module = igraph::V(graph)$source_module,
    module_rank = igraph::V(graph)$module_rank,
    module_label = igraph::V(graph)$module_label,
    module_display = igraph::V(graph)$module_display,
    functional_core_class =
      igraph::V(graph)$functional_core_class,
    functional_compartment = CLASS_LABELS[
      igraph::V(graph)$functional_core_class
    ],
    source_degree = igraph::V(graph)$source_degree,
    source_strength_abs = igraph::V(graph)$source_strength_abs,
    plot_degree = igraph::V(graph)$plot_degree,
    plot_strength = igraph::V(graph)$plot_strength,
    KO_function = igraph::V(graph)$KO_function
  )

  cat("[5] Exportando source data y archivos para Gephi/Cytoscape...\n")
  nodes_out <- data.table::copy(layout_table)
  nodes_out[, module_color := module_colors[module_display]]
  nodes_out[, core_class_color := CLASS_COLORS[functional_core_class]]
  edges_out <- data.table::copy(edges_plot)

  write_tsv(
    curated_modules,
    file.path(TABLES_DIR, "172_curated_modules.tsv")
  )
  write_tsv(
    nodes_out,
    file.path(TABLES_DIR, "172_network_nodes_for_plot.tsv")
  )
  data.table::fwrite(
    edges_out,
    file.path(TABLES_DIR, "172_network_edges_for_plot.tsv.gz"),
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )
  write_tsv(
    layout_table,
    file.path(TABLES_DIR, "172_network_layout_coordinates.tsv")
  )
  write_tsv(
    top_labels,
    file.path(TABLES_DIR, "172_network_top_labeled_nodes.tsv")
  )

  module_palette <- data.frame(
    module_display = names(module_colors),
    color = unname(module_colors),
    stringsAsFactors = FALSE
  )
  class_palette <- data.frame(
    functional_core_class = CLASS_LEVELS,
    functional_compartment = unname(CLASS_LABELS[CLASS_LEVELS]),
    color = unname(CLASS_COLORS[CLASS_LEVELS]),
    observed_in_graph = CLASS_LEVELS %in% observed_classes,
    stringsAsFactors = FALSE
  )
  write_tsv(
    module_palette,
    file.path(TABLES_DIR, "172_module_palette_used.tsv")
  )
  write_tsv(
    class_palette,
    file.path(TABLES_DIR, "172_core_class_palette_used.tsv")
  )

  gephi_nodes <- data.table::copy(nodes_out)
  data.table::setnames(gephi_nodes, "KO", "Id")
  gephi_nodes[, Label := Id]
  gephi_edges <- data.table::copy(edges_out)
  gephi_edges[, Source := KO1]
  gephi_edges[, Target := KO2]
  gephi_edges[, Type := "Undirected"]
  gephi_edges[, Weight := mean_abs_weight]

  data.table::fwrite(
    gephi_nodes,
    file.path(TABLES_DIR, "172_gephi_nodes.csv"),
    sep = ",",
    quote = TRUE,
    na = ""
  )
  data.table::fwrite(
    gephi_edges[, .(
      Source, Target, Type, Weight, consensus_weight,
      edge_stability, sign_consistency, is_internal_module_edge,
      module_label_KO1, module_label_KO2
    )],
    file.path(TABLES_DIR, "172_gephi_edges.csv"),
    sep = ",",
    quote = TRUE,
    na = ""
  )

  caption_numbers <- data.frame(
    metric = c(
      "source_candidate_KOs",
      "source_consensus_edges",
      "source_connected_KOs",
      "curated_major_modules",
      "KOs_in_curated_major_modules",
      "edges_between_curated_major_module_KOs",
      "internal_edges",
      "between_module_edges",
      "absolute_global_core_KOs_in_curated_graph",
      "near_global_KOs_in_curated_graph",
      "realm_core_only_KOs_in_curated_graph",
      "multi_realm_non_core_KOs_in_curated_graph",
      "realm_specific_KOs_in_curated_graph"
    ),
    value = c(
      nrow(nodes),
      raw_edge_n,
      sum(nodes$connected),
      nrow(curated_modules),
      nrow(nodes_plot),
      nrow(edges_plot),
      sum(edges_plot$is_internal_module_edge),
      sum(!edges_plot$is_internal_module_edge),
      sum(nodes_plot$functional_core_class == "absolute_global_core"),
      sum(nodes_plot$functional_core_class == "near_global"),
      sum(nodes_plot$functional_core_class == "realm_core_only"),
      sum(nodes_plot$functional_core_class == "multi_realm_non_core"),
      sum(nodes_plot$functional_core_class == "realm_specific")
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(
    caption_numbers,
    file.path(TABLES_DIR, "172_network_caption_numbers.tsv")
  )

  cat("[6] Generando figuras editoriales...\n")
  plot_network <- function(
    color_mode = c("module", "prevalence"),
    show_labels = FALSE,
    internal_edges_only = FALSE
  ) {
    color_mode <- match.arg(color_mode)
    old_parameters <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_parameters), add = TRUE)
    graphics::par(mar = c(0, 0, 0, 0), xpd = NA)

    graph_plot <- graph
    layout_plot <- layout_xy
    if (internal_edges_only) {
      keep_edges <- which(igraph::E(graph_plot)$is_internal_module_edge)
      graph_plot <- igraph::subgraph.edges(
        graph_plot,
        eids = keep_edges,
        delete.vertices = FALSE
      )
      graph_plot <- igraph::induced_subgraph(
        graph_plot,
        vids = igraph::V(graph_plot)[igraph::degree(graph_plot) > 0]
      )
      layout_plot <- layout_xy[
        match(igraph::V(graph_plot)$name, igraph::V(graph)$name),
        ,
        drop = FALSE
      ]
    }

    vertex_colors <- if (color_mode == "module") {
      igraph::V(graph_plot)$plot_color_module
    } else {
      igraph::V(graph_plot)$plot_color_class
    }
    vertex_labels <- if (show_labels) {
      igraph::V(graph_plot)$plot_label
    } else {
      NA_character_
    }

    igraph::plot.igraph(
      graph_plot,
      layout = layout_plot,
      vertex.color = vertex_colors,
      vertex.frame.color = "#252525",
      vertex.frame.width = 0.25,
      vertex.size = igraph::V(graph_plot)$plot_size,
      vertex.label = vertex_labels,
      vertex.label.cex = 0.62,
      vertex.label.color = "#111111",
      edge.color = igraph::E(graph_plot)$plot_color,
      edge.width = igraph::E(graph_plot)$plot_width,
      edge.curved = 0.03,
      asp = 0,
      rescale = TRUE
    )

    if (color_mode == "module") {
      graphics::legend(
        "bottomleft",
        legend = names(module_colors),
        col = module_colors,
        pch = 19,
        pt.cex = 1.2,
        bty = "n",
        cex = 0.72
      )
    } else {
      graphics::legend(
        "bottomleft",
        legend = unname(CLASS_LABELS[observed_classes]),
        col = observed_class_colors,
        pch = 19,
        pt.cex = 1.2,
        bty = "n",
        cex = 0.72
      )
    }
  }

  save_base_plot(
    file.path(FIGURES_DIR, "172_network_by_module_unlabeled"),
    PLOT_WIDTH,
    PLOT_HEIGHT,
    function() plot_network("module", FALSE, FALSE),
    DPI
  )
  save_base_plot(
    file.path(FIGURES_DIR, "172_network_by_module_labeled_hubs"),
    PLOT_WIDTH,
    PLOT_HEIGHT,
    function() plot_network("module", TRUE, FALSE),
    DPI
  )
  save_base_plot(
    file.path(FIGURES_DIR, "172_network_by_core_class"),
    PLOT_WIDTH,
    PLOT_HEIGHT,
    function() plot_network("prevalence", TRUE, FALSE),
    DPI
  )
  save_base_plot(
    file.path(FIGURES_DIR, "172_network_internal_module_edges"),
    PLOT_WIDTH,
    PLOT_HEIGHT,
    function() plot_network("module", TRUE, TRUE),
    DPI
  )

  cat("[7] Ejecutando validaciones estrictas...\n")
  candidate_class_counts <- table(factor(
    nodes$functional_core_class,
    levels = CLASS_LEVELS
  ))
  endpoint_missing <- length(setdiff(
    unique(c(edges_plot$KO1, edges_plot$KO2)),
    nodes_plot$KO
  ))
  duplicate_edge_keys <- anyDuplicated(edges_plot$edge_key)
  duplicate_nodes <- anyDuplicated(nodes_plot$KO)
  unexpected_classes <- sum(
    !nodes$functional_core_class %in% CLASS_LEVELS
  )
  nonfinite_layout <- sum(!is.finite(layout_xy))

  validation <- data.frame(
    metric = c(
      "n_source_consensus_edges",
      "n_source_candidate_KOs",
      "n_source_connected_KOs",
      "n_source_major_modules",
      "n_source_major_module_KOs",
      "n_curated_modules",
      "n_plot_nodes",
      "n_nodes_dropped_without_major_edge",
      "n_duplicate_plot_nodes",
      "n_duplicate_plot_edge_keys",
      "n_missing_edge_endpoints",
      "n_unexpected_core_classes",
      "n_nonfinite_layout_coordinates"
    ),
    observed = c(
      raw_edge_n,
      raw_node_n,
      sum(nodes$connected),
      sum(modules$major_module),
      sum(modules$n_KOs[modules$major_module]),
      nrow(curated_modules),
      nrow(nodes_plot),
      length(nodes_dropped_no_major_edge),
      as.integer(duplicate_nodes > 0L),
      as.integer(duplicate_edge_keys > 0L),
      endpoint_missing,
      unexpected_classes,
      nonfinite_layout
    ),
    expected = c(
      EXPECTED_SOURCE_EDGES,
      EXPECTED_CANDIDATE_KOS,
      EXPECTED_CONNECTED_KOS,
      EXPECTED_MAJOR_MODULES,
      EXPECTED_MAJOR_MODULE_KOS,
      EXPECTED_MAJOR_MODULES,
      EXPECTED_MAJOR_MODULE_KOS,
      0, 0, 0, 0, 0, 0
    ),
    stringsAsFactors = FALSE
  )
  validation$pass <- validation$observed == validation$expected
  write_tsv(
    validation,
    file.path(TABLES_DIR, "172_expected_validation.tsv")
  )

  if (STRICT_EXPECTED && any(!validation$pass)) {
    stop(
      "Fallo la validacion estricta: ",
      paste(validation$metric[!validation$pass], collapse = ", ")
    )
  }

  figure_paths <- c(
    file.path(FIGURES_DIR, "172_network_by_module_unlabeled.pdf"),
    file.path(FIGURES_DIR, "172_network_by_module_unlabeled.png"),
    file.path(FIGURES_DIR, "172_network_by_module_unlabeled.tiff"),
    file.path(FIGURES_DIR, "172_network_by_module_labeled_hubs.pdf"),
    file.path(FIGURES_DIR, "172_network_by_module_labeled_hubs.png"),
    file.path(FIGURES_DIR, "172_network_by_module_labeled_hubs.tiff"),
    file.path(FIGURES_DIR, "172_network_by_core_class.pdf"),
    file.path(FIGURES_DIR, "172_network_by_core_class.png"),
    file.path(FIGURES_DIR, "172_network_by_core_class.tiff"),
    file.path(FIGURES_DIR, "172_network_internal_module_edges.pdf"),
    file.path(FIGURES_DIR, "172_network_internal_module_edges.png"),
    file.path(FIGURES_DIR, "172_network_internal_module_edges.tiff")
  )
  figure_validation <- data.frame(
    file = figure_paths,
    exists = file.exists(figure_paths),
    size_bytes = ifelse(
      file.exists(figure_paths),
      file.info(figure_paths)$size,
      0
    ),
    stringsAsFactors = FALSE
  )
  figure_validation$pass <-
    figure_validation$exists & figure_validation$size_bytes > 0
  write_tsv(
    figure_validation,
    file.path(TABLES_DIR, "172_figure_validation.tsv")
  )
  if (any(!figure_validation$pass)) stop("Fallo la escritura de figuras.")

  summary_table <- data.frame(
    item = c(
      "status",
      "run171",
      "method",
      "edge_interpretation",
      "n_source_candidate_KOs",
      "n_source_connected_KOs",
      "n_source_consensus_edges",
      "n_curated_modules",
      "n_plot_nodes",
      "n_plot_edges",
      "n_internal_edges",
      "n_between_module_edges",
      "graph_density_curated",
      "n_absolute_global_core_plot_nodes",
      "n_near_global_plot_nodes",
      "n_realm_core_only_plot_nodes",
      "n_multi_realm_non_core_plot_nodes",
      "n_realm_specific_plot_nodes",
      "layout",
      "seed",
      "niter_layout",
      "output_dir"
    ),
    value = c(
      "PASS",
      RUN171,
      "Fixed 171 Louvain partition; deterministic Fruchterman-Reingold layout",
      "Robust CLR co-abundance association; not biochemical interaction",
      nrow(nodes),
      sum(nodes$connected),
      raw_edge_n,
      nrow(curated_modules),
      nrow(nodes_plot),
      nrow(edges_plot),
      sum(edges_plot$is_internal_module_edge),
      sum(!edges_plot$is_internal_module_edge),
      igraph::edge_density(graph, loops = FALSE),
      sum(nodes_plot$functional_core_class == "absolute_global_core"),
      sum(nodes_plot$functional_core_class == "near_global"),
      sum(nodes_plot$functional_core_class == "realm_core_only"),
      sum(nodes_plot$functional_core_class == "multi_realm_non_core"),
      sum(nodes_plot$functional_core_class == "realm_specific"),
      "Fruchterman-Reingold weighted by mean absolute proportionality",
      SEED,
      NITER_LAYOUT,
      OUTDIR
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(summary_table, file.path(TABLES_DIR, "172_summary.tsv"))

  saveRDS(
    list(
      graph = graph,
      layout = layout_xy,
      nodes = as.data.frame(nodes_out),
      edges = as.data.frame(edges_out),
      curated_modules = as.data.frame(curated_modules),
      module_colors = module_colors,
      core_class_colors = CLASS_COLORS,
      parameters = list(
        run171 = RUN171,
        min_module_size = MIN_MODULE_SIZE,
        top_labels_per_module = TOP_LABELS_PER_MODULE,
        seed = SEED,
        niter_layout = NITER_LAYOUT
      )
    ),
    file.path(RDS_DIR, "172_curated_network_layout_objects.rds")
  )

  writeLines("PASS", STATUS_FILE)

  cat("\n============================================================\n")
  cat("STATUS: PASS\n")
  cat("Output dir: ", OUTDIR, "\n", sep = "")
  cat("Source consensus edges: ", raw_edge_n, "\n", sep = "")
  cat("Curated modules: ", nrow(curated_modules), "\n", sep = "")
  cat("Plot nodes: ", nrow(nodes_plot), "\n", sep = "")
  cat("Plot edges: ", nrow(edges_plot), "\n", sep = "")
  cat(
    "Absolute/near-global/realm-core/multi-realm/realm-specific plot nodes: ",
    paste(as.integer(table(factor(
      nodes_plot$functional_core_class,
      levels = CLASS_LEVELS
    ))), collapse = "/"), "\n",
    sep = ""
  )
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
      file.path(LOGS_DIR, "172_error.txt")
    )
    1L
  }
)

quit(save = "no", status = status)
