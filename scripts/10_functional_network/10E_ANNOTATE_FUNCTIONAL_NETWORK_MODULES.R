#!/usr/bin/env Rscript

options(warn = 1)

suppressPackageStartupMessages({
  for (package_i in c("data.table", "ggplot2")) {
    if (!requireNamespace(package_i, quietly = TRUE)) {
      stop("Se requiere el paquete '", package_i, "'.")
    }
  }
})

# =============================================================================
# 173_ANNOTATE_FUNCTIONAL_NETWORK_MODULES_FINAL.R
#
# Objetivo
#   Anotar e interpretar los seis modulos de la red funcional validada por los
#   scripts 171-172 bajo la clasificacion funcional final del script 153.
#
# Reproducibilidad
#   Las anotaciones KO -> KEGG MODULE y las categorias funcionales se recuperan
#   del snapshot ya usado por el script 144. El script NO consulta KEGGREST ni
#   modifica la particion de red. Del snapshot 144 solo se usan anotaciones
#   funcionales; las clases de persistencia se conservan desde 153/171.
#
# Interpretacion
#   Los modulos son programas de coabundancia del potencial funcional codificado.
#   No representan actividad metabolica, expresion, causalidad ni interacciones
#   bioquimicas directas.
# =============================================================================

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1L])
source(file.path(dirname(normalizePath(script_file)), "lib", "functional_category_rules.R"))

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

RUN172 <- get_arg(
  "--run172",
  read_latest("172_CURATE_FUNCTIONAL_NETWORK_LAYOUT_FINAL")
)
RUN144_REFERENCE <- get_arg(
  "--run144_reference",
  read_latest("144_annotate_functional_network_modules")
)

TOP_HUBS_PER_MODULE <- as.integer(get_arg("--top_hubs_per_module", "15"))
TOP_HUBS_FIGURE <- as.integer(get_arg("--top_hubs_figure", "5"))
TOP_KEGG_PER_MODULE <- as.integer(get_arg("--top_kegg_per_module", "8"))
DPI <- as.integer(get_arg("--dpi", "600"))
MAKE_FIGURES <- as_flag(get_arg("--make_figures", "1"))
STRICT_EXPECTED <- as_flag(get_arg("--strict_expected", "0"))

EXPECTED_NODES <- as.integer(get_arg("--expected_nodes", "237"))
EXPECTED_EDGES <- as.integer(get_arg("--expected_edges", "1162"))
EXPECTED_MODULES <- as.integer(get_arg("--expected_modules", "6"))
EXPECTED_ABSOLUTE <- as.integer(get_arg("--expected_absolute", "9"))
EXPECTED_NEAR_GLOBAL <- as.integer(get_arg("--expected_near_global", "195"))
EXPECTED_REALM_CORE <- as.integer(get_arg("--expected_realm_core", "33"))
EXPECTED_MULTI_REALM <- as.integer(get_arg("--expected_multi_realm", "0"))
EXPECTED_REALM_SPECIFIC <- as.integer(get_arg("--expected_realm_specific", "0"))

SCRIPT_NUMBER <- "173"
SCRIPT_BASENAME <-
  "ANNOTATE_FUNCTIONAL_NETWORK_MODULES_FINAL"
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

dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOGS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)
writeLines(OUTDIR, LATEST_FILE)
writeLines("RUNNING", STATUS_FILE)

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

FUNCTIONAL_CATEGORY_ORDER <- c(
  "Anaerobic redox / energy",
  "Transport / membrane systems",
  "Information processing",
  "Amino acid / cofactor biosynthesis",
  "Carbon / carbohydrate metabolism",
  "Signal / stress / metals",
  "Uncharacterized / other"
)
FUNCTIONAL_CATEGORY_COLORS <- c(
  "Anaerobic redox / energy" = "#8C510A",
  "Transport / membrane systems" = "#01665E",
  "Information processing" = "#762A83",
  "Amino acid / cofactor biosynthesis" = "#5AAE61",
  "Carbon / carbohydrate metabolism" = "#1F78B4",
  "Signal / stress / metals" = "#E66101",
  "Uncharacterized / other" = "#BDBDBD"
)

clean_chr <- function(x) {
  y <- trimws(as.character(x))
  y[
    is.na(y) | y == "" |
      tolower(y) %in% c("na", "n/a", "nan", "null", "none", "unknown")
  ] <- NA_character_
  y
}

safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

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

shorten_label <- function(x, max_chars = 70L) {
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

wrap_label <- function(x, width = 32L) {
  vapply(
    as.character(x),
    function(value_i) paste(strwrap(value_i, width = width), collapse = "\n"),
    character(1)
  )
}

save_gg <- function(plot, file_base, width, height, dpi) {
  if (!MAKE_FIGURES) return(invisible(NULL))
  ggplot2::ggsave(
    paste0(file_base, ".pdf"),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf,
    bg = "white"
  )
  ggplot2::ggsave(
    paste0(file_base, ".png"),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white"
  )
  ggplot2::ggsave(
    paste0(file_base, ".tiff"),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    compression = "lzw",
    bg = "white"
  )
}

theme_clean <- function(base_size = 12) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(color = "#111111"),
      axis.title = ggplot2::element_text(color = "#111111", face = "bold"),
      strip.background = ggplot2::element_rect(fill = "#F1F3F4", color = NA),
      strip.text = ggplot2::element_text(face = "bold", color = "#111111"),
      legend.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", color = "#111111"),
      plot.subtitle = ggplot2::element_text(color = "#555555"),
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
    )
}

make_complete_composition <- function(nodes, grouping_column, levels_vector) {
  grouping_values <- sort(unique(as.character(nodes[[grouping_column]])))
  grid <- data.table::CJ(
    module_label = grouping_values,
    category = levels_vector,
    unique = TRUE
  )
  observed <- nodes[
    ,
    .(n_KOs = .N),
    by = c(grouping_column, "category")
  ]
  data.table::setnames(observed, grouping_column, "module_label")
  output <- merge(
    grid,
    observed,
    by = c("module_label", "category"),
    all.x = TRUE,
    sort = FALSE
  )
  output[is.na(n_KOs), n_KOs := 0L]
  output[, module_n_KOs := sum(n_KOs), by = module_label]
  output[, proportion := n_KOs / module_n_KOs]
  output
}

main <- function() {
  validate_run(RUN172, "RUN172")
  validate_run(RUN144_REFERENCE, "RUN144_REFERENCE")
  if (
    TOP_HUBS_PER_MODULE < 1L || TOP_HUBS_FIGURE < 1L ||
      TOP_KEGG_PER_MODULE < 1L || DPI < 72L
  ) {
    stop("Parametros invalidos.")
  }

  log_path <- file.path(LOGS_DIR, "173_log.txt")
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
  cat("173_ANNOTATE_FUNCTIONAL_NETWORK_MODULES_FINAL\n")
  cat("Started: ", as.character(Sys.time()), "\n", sep = "")
  cat("RUN172: ", RUN172, "\n", sep = "")
  cat("RUN144_REFERENCE: ", RUN144_REFERENCE, "\n", sep = "")
  cat("TOP_HUBS_PER_MODULE: ", TOP_HUBS_PER_MODULE, "\n", sep = "")
  cat("TOP_KEGG_PER_MODULE: ", TOP_KEGG_PER_MODULE, "\n", sep = "")
  cat("STRICT_EXPECTED: ", STRICT_EXPECTED, "\n", sep = "")
  cat("OUTDIR: ", OUTDIR, "\n", sep = "")
  cat("============================================================\n\n")

  nodes_path <- file.path(
    RUN172,
    "tables",
    "172_network_nodes_for_plot.tsv"
  )
  edges_path <- file.path(
    RUN172,
    "tables",
    "172_network_edges_for_plot.tsv.gz"
  )
  modules_path <- file.path(RUN172, "tables", "172_curated_modules.tsv")
  reference_nodes_path <- file.path(
    RUN144_REFERENCE,
    "tables",
    "144_network_nodes_annotated.tsv"
  )
  reference_kegg_path <- file.path(
    RUN144_REFERENCE,
    "tables",
    "144_KO_to_KEGG_MODULE_long.tsv"
  )

  cat("[1] Leyendo red curada y snapshot funcional de referencia...\n")
  nodes <- read_tsv(nodes_path)
  edges <- read_tsv(edges_path)
  modules <- read_tsv(modules_path)
  reference_nodes <- read_tsv(reference_nodes_path)
  reference_kegg <- read_tsv(reference_kegg_path)

  require_columns(
    nodes,
    c(
      "KO", "x", "y", "module_label", "module_display",
      "functional_core_class", "source_degree",
      "source_strength_abs", "KO_function"
    ),
    "164 nodes"
  )
  require_columns(
    edges,
    c(
      "KO1", "KO2", "consensus_weight", "mean_abs_weight",
      "edge_stability", "sign_consistency"
    ),
    "164 edges"
  )
  require_columns(
    modules,
    c("module", "module_label", "module_display", "n_KOs"),
    "164 modules"
  )
  require_columns(
    reference_nodes,
    c("KO", "functional_category"),
    "144 reference nodes"
  )
  require_columns(
    reference_kegg,
    c("KO", "KEGG_module_id", "KEGG_module_name"),
    "144 KO-to-KEGG snapshot"
  )

  cat("    current nodes: ", nrow(nodes), "\n", sep = "")
  cat("    current edges: ", nrow(edges), "\n", sep = "")
  cat("    current modules: ", nrow(modules), "\n", sep = "")
  cat("    reference annotated nodes: ", nrow(reference_nodes), "\n", sep = "")
  cat("    reference KO-KEGG links: ", nrow(reference_kegg), "\n", sep = "")

  cat("[2] Alineando anotaciones exclusivamente por KO...\n")
  nodes[, KO := clean_chr(KO)]
  nodes[, x := safe_num(x)]
  nodes[, y := safe_num(y)]
  nodes[, module_label := clean_chr(module_label)]
  nodes[, module_display := clean_chr(module_display)]
  nodes[, functional_core_class := clean_chr(
    functional_core_class
  )]
  nodes[, source_degree := safe_num(source_degree)]
  nodes[, source_strength_abs := safe_num(source_strength_abs)]
  nodes[, KO_function := clean_chr(KO_function)]
  nodes[is.na(KO_function), KO_function := KO]
  nodes <- nodes[!is.na(KO) & !duplicated(KO)]

  reference_nodes[, KO := clean_chr(KO)]
  reference_nodes[, functional_category_reference := clean_chr(
    functional_category
  )]
  if ("KO_function" %in% names(reference_nodes)) {
    reference_nodes[, KO_function_reference := clean_chr(KO_function)]
  } else {
    reference_nodes[, KO_function_reference := NA_character_]
  }
  optional_reference_columns <- intersect(
    c("n_KEGG_modules", "KEGG_module_ids", "KEGG_module_names"),
    names(reference_nodes)
  )
  reference_node_subset <- reference_nodes[
    !is.na(KO) & !duplicated(KO),
    c(
      "KO",
      "functional_category_reference",
      "KO_function_reference",
      optional_reference_columns
    ),
    with = FALSE
  ]
  nodes <- merge(
    nodes,
    reference_node_subset,
    by = "KO",
    all.x = TRUE,
    sort = FALSE
  )
  nodes[
    is.na(KO_function) & !is.na(KO_function_reference),
    KO_function := KO_function_reference
  ]
  nodes[, combined_category := functional_category_reference]
  category_scores <- lapply(nodes$KO_function,function(z) score_functional_category(z,""))
  nodes[, functional_category := vapply(category_scores,function(z) z$selected,character(1))]
  nodes[, functional_category_reference := NULL]
  nodes[, KO_function_reference := NULL]

  invalid_categories <- !nodes$functional_category %in%
    FUNCTIONAL_CATEGORY_ORDER
  nodes[invalid_categories, functional_category := "Uncharacterized / other"]
  nodes[, functional_compartment := CLASS_LABELS[
    functional_core_class
  ]]
  # Module IDs can change after cohort updates; historical biological labels
  # must not be transferred to newly inferred communities.
  nodes[, module_program := paste("Co-abundance module", module_label)]

  reference_match_n <- sum(
    nodes$KO %in% reference_nodes$KO & !is.na(nodes$functional_category)
  )
  cat(
    "    nodes matched to functional snapshot: ",
    reference_match_n, "/", nrow(nodes), "\n",
    sep = ""
  )

  cat("[3] Filtrando el snapshot KO-KEGG a los nodos actuales...\n")
  reference_kegg[, KO := clean_chr(KO)]
  reference_kegg[, KEGG_module_id := clean_chr(KEGG_module_id)]
  reference_kegg[, KEGG_module_name := clean_chr(KEGG_module_name)]
  kegg_map <- reference_kegg[
    KO %in% nodes$KO & !is.na(KO) & !is.na(KEGG_module_id)
  ]
  kegg_map <- kegg_map[
    !duplicated(paste(KO, KEGG_module_id, sep = "||"))
  ]

  node_module_lookup <- nodes[, .(KO, module_label, module_display)]
  kegg_map <- merge(
    kegg_map[, .(KO, KEGG_module_id, KEGG_module_name)],
    node_module_lookup,
    by = "KO",
    all.x = TRUE,
    sort = FALSE
  )

  if (STRICT_EXPECTED && nrow(kegg_map) == 0L) {
    stop("El snapshot 144 no contiene enlaces KO-KEGG para los nodos actuales.")
  }
  if (
    STRICT_EXPECTED &&
      data.table::uniqueN(kegg_map$module_label) != EXPECTED_MODULES
  ) {
    stop("El snapshot KO-KEGG no cubre los seis modulos actuales.")
  }

  kegg_by_node <- if (nrow(kegg_map) > 0L) {
    kegg_map[
      ,
      .(
        n_KEGG_modules_snapshot = data.table::uniqueN(KEGG_module_id),
        KEGG_module_ids_snapshot = paste(
          sort(unique(KEGG_module_id)),
          collapse = ";"
        ),
        KEGG_module_names_snapshot = paste(
          sort(unique(KEGG_module_name[!is.na(KEGG_module_name)])),
          collapse = " | "
        )
      ),
      by = KO
    ]
  } else {
    data.table::data.table(
      KO = character(),
      n_KEGG_modules_snapshot = integer(),
      KEGG_module_ids_snapshot = character(),
      KEGG_module_names_snapshot = character()
    )
  }
  nodes <- merge(nodes, kegg_by_node, by = "KO", all.x = TRUE, sort = FALSE)
  nodes[is.na(n_KEGG_modules_snapshot), n_KEGG_modules_snapshot := 0L]

  cat("    retained KO-KEGG links: ", nrow(kegg_map), "\n", sep = "")
  cat(
    "    nodes with at least one KEGG MODULE: ",
    sum(nodes$n_KEGG_modules_snapshot > 0L), "\n",
    sep = ""
  )

  cat("[4] Anotando aristas y calculando hubs...\n")
  edges[, KO1 := clean_chr(KO1)]
  edges[, KO2 := clean_chr(KO2)]
  edges[, consensus_weight := safe_num(consensus_weight)]
  edges[, mean_abs_weight := safe_num(mean_abs_weight)]
  edges[, edge_stability := safe_num(edge_stability)]
  edges[, sign_consistency := safe_num(sign_consistency)]
  existing_edge_annotations <- intersect(
    c(
      "module_KO1", "module_KO2", "module_label_KO1",
      "module_label_KO2", "is_internal_module_edge"
    ),
    names(edges)
  )
  if (length(existing_edge_annotations) > 0L) {
    edges[, (existing_edge_annotations) := NULL]
  }
  edges <- edges[
    KO1 %in% nodes$KO & KO2 %in% nodes$KO & KO1 != KO2
  ]

  annotation_lookup <- nodes[, .(
    KO,
    module_label,
    module_program,
    functional_core_class,
    functional_category,
    KO_function
  )]
  edges <- merge(
    edges,
    annotation_lookup,
    by.x = "KO1",
    by.y = "KO",
    all.x = TRUE,
    sort = FALSE
  )
  data.table::setnames(
    edges,
    c(
      "module_label", "module_program", "functional_core_class",
      "functional_category", "KO_function"
    ),
    c(
      "module_label_KO1", "module_program_KO1", "core_class_KO1",
      "functional_category_KO1", "KO_function_KO1"
    )
  )
  edges <- merge(
    edges,
    annotation_lookup,
    by.x = "KO2",
    by.y = "KO",
    all.x = TRUE,
    sort = FALSE
  )
  data.table::setnames(
    edges,
    c(
      "module_label", "module_program", "functional_core_class",
      "functional_category", "KO_function"
    ),
    c(
      "module_label_KO2", "module_program_KO2", "core_class_KO2",
      "functional_category_KO2", "KO_function_KO2"
    )
  )
  edges[, is_internal_module_edge := module_label_KO1 == module_label_KO2]

  data.table::setorder(
    nodes,
    module_label,
    -source_strength_abs,
    -source_degree,
    KO
  )
  top_hubs <- nodes[
    ,
    head(.SD, TOP_HUBS_PER_MODULE),
    by = module_label
  ]
  top_hubs[, hub_rank := seq_len(.N), by = module_label]
  top_hubs[, hub_label := paste0(
    KO,
    " - ",
    shorten_label(KO_function, 80L)
  )]

  cat("[5] Resumiendo clases de persistencia y funcion...\n")
  core_nodes <- data.table::copy(nodes)
  core_nodes[, category := functional_core_class]
  core_composition <- make_complete_composition(
    core_nodes,
    "module_label",
    CLASS_LEVELS
  )
  data.table::setnames(
    core_composition,
    "category",
    "functional_core_class"
  )
  core_composition[, functional_compartment := CLASS_LABELS[
    functional_core_class
  ]]

  functional_nodes <- data.table::copy(nodes)
  functional_nodes[, category := functional_category]
  functional_composition <- make_complete_composition(
    functional_nodes,
    "module_label",
    FUNCTIONAL_CATEGORY_ORDER
  )
  data.table::setnames(
    functional_composition,
    "category",
    "functional_category"
  )

  module_totals <- nodes[
    ,
    .(
      n_KOs = .N,
      n_absolute_global_core = sum(functional_core_class == "absolute_global_core"),
      n_near_global = sum(functional_core_class == "near_global"),
      n_realm_core_only = sum(functional_core_class == "realm_core_only"),
      n_multi_realm_non_core = sum(functional_core_class == "multi_realm_non_core"),
      n_realm_specific = sum(functional_core_class == "realm_specific"),
      mean_degree = mean(source_degree),
      mean_strength_abs = mean(source_strength_abs),
      top_hub_KO = KO[which.max(source_strength_abs)],
      top_hub_function = KO_function[which.max(source_strength_abs)]
    ),
    by = .(module_label, module_display, module_program)
  ]
  module_totals[, percent_absolute_global_core := 100 * n_absolute_global_core / n_KOs]
  module_totals[, percent_near_global := 100 * n_near_global / n_KOs]
  module_totals[, percent_realm_core_only := 100 * n_realm_core_only / n_KOs]
  module_totals[, percent_multi_realm_non_core := 100 * n_multi_realm_non_core / n_KOs]
  module_totals[, percent_realm_specific := 100 * n_realm_specific / n_KOs]

  dominant_function <- functional_composition[
    order(module_label, -n_KOs, functional_category),
    .SD[1L],
    by = module_label
  ][
    ,
    .(
      module_label,
      dominant_functional_category = functional_category,
      dominant_functional_n_KOs = n_KOs,
      dominant_functional_proportion = proportion
    )
  ]
  module_summary <- merge(
    module_totals,
    dominant_function,
    by = "module_label",
    all.x = TRUE,
    sort = FALSE
  )
  data.table::setorder(module_summary, module_label)

  cat("[6] Resumiendo KEGG MODULE y solapamiento entre programas...\n")
  kegg_module_summary <- if (nrow(kegg_map) > 0L) {
    kegg_map[
      ,
      .(
        n_KOs = data.table::uniqueN(KO),
        KO_examples = paste(head(sort(unique(KO)), 8L), collapse = ";")
      ),
      by = .(
        module_label,
        module_display,
        KEGG_module_id,
        KEGG_module_name
      )
    ]
  } else {
    data.table::data.table()
  }
  if (nrow(kegg_module_summary) > 0L) {
    data.table::setorder(
      kegg_module_summary,
      module_label,
      -n_KOs,
      KEGG_module_id
    )
    kegg_module_summary[
      ,
      rank_within_network_module := seq_len(.N),
      by = module_label
    ]
  }
  top_kegg_module_summary <- if (nrow(kegg_module_summary) > 0L) {
    kegg_module_summary[
      rank_within_network_module <= TOP_KEGG_PER_MODULE
    ]
  } else {
    data.table::data.table()
  }

  process_sets <- if (nrow(kegg_map) > 0L) {
    kegg_map[
      ,
      .(
        KEGG_module_set = list(sort(unique(KEGG_module_id))),
        KEGG_module_name_set = list(sort(unique(
          KEGG_module_name[!is.na(KEGG_module_name)]
        )))
      ),
      by = module_label
    ]
  } else {
    data.table::data.table()
  }

  process_overlap_rows <- list()
  overlap_index <- 1L
  if (nrow(process_sets) > 0L) {
    for (i in seq_len(nrow(process_sets))) {
      for (j in seq_len(nrow(process_sets))) {
        set_i <- process_sets$KEGG_module_set[[i]]
        set_j <- process_sets$KEGG_module_set[[j]]
        shared_ids <- intersect(set_i, set_j)
        union_ids <- union(set_i, set_j)
        shared_rows <- kegg_map[KEGG_module_id %in% shared_ids]
        shared_names <- sort(unique(
          shared_rows$KEGG_module_name[!is.na(shared_rows$KEGG_module_name)]
        ))
        process_overlap_rows[[overlap_index]] <- data.table::data.table(
          module_i = process_sets$module_label[i],
          module_j = process_sets$module_label[j],
          n_processes_i = length(set_i),
          n_processes_j = length(set_j),
          shared_n = length(shared_ids),
          union_n = length(union_ids),
          jaccard = ifelse(
            length(union_ids) == 0L,
            NA_real_,
            length(shared_ids) / length(union_ids)
          ),
          shared_KEGG_module_ids = ifelse(
            length(shared_ids) == 0L,
            NA_character_,
            paste(shared_ids, collapse = ";")
          ),
          shared_KEGG_module_names = ifelse(
            length(shared_names) == 0L,
            NA_character_,
            paste(shared_names, collapse = " | ")
          )
        )
        overlap_index <- overlap_index + 1L
      }
    }
  }
  process_overlap <- if (length(process_overlap_rows) > 0L) {
    data.table::rbindlist(process_overlap_rows)
  } else {
    data.table::data.table()
  }

  cat("[7] Escribiendo tablas anotadas y source data...\n")
  write_tsv(nodes, file.path(TABLES_DIR, "173_network_nodes_annotated.tsv"))
  data.table::fwrite(
    edges,
    file.path(TABLES_DIR, "173_network_edges_annotated.tsv.gz"),
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )
  write_tsv(
    module_summary,
    file.path(TABLES_DIR, "173_network_module_functional_summary.tsv")
  )
  write_tsv(
    top_hubs,
    file.path(TABLES_DIR, "173_network_module_top_hubs_annotated.tsv")
  )
  write_tsv(
    core_composition,
    file.path(TABLES_DIR, "173_module_core_composition.tsv")
  )
  write_tsv(
    functional_composition,
    file.path(TABLES_DIR, "173_module_functional_category_composition.tsv")
  )
  write_tsv(
    kegg_map,
    file.path(TABLES_DIR, "173_KO_to_KEGG_MODULE_snapshot_long.tsv")
  )
  write_tsv(
    kegg_module_summary,
    file.path(TABLES_DIR, "173_network_module_kegg_module_summary.tsv")
  )
  write_tsv(
    top_kegg_module_summary,
    file.path(TABLES_DIR, "173_network_module_top_kegg_modules.tsv")
  )
  write_tsv(
    process_overlap,
    file.path(TABLES_DIR, "173_network_module_process_overlap.tsv")
  )

  cat("[8] Generando figuras editoriales...\n")
  prevalence_plot_data <- data.table::copy(core_composition)
  prevalence_plot_data[, module_label := factor(
    module_label,
    levels = sort(unique(module_label))
  )]
  prevalence_plot_data[, functional_core_class := factor(
    functional_core_class,
    levels = CLASS_LEVELS
  )]
  prevalence_plot_data[, tile_label := ifelse(
    n_KOs == 0L,
    "",
    paste0(n_KOs, "\n", round(100 * proportion), "%")
  )]
  p_prevalence <- ggplot2::ggplot(
    prevalence_plot_data,
    ggplot2::aes(
      x = functional_core_class,
      y = module_label,
      fill = proportion
    )
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.6) +
    ggplot2::geom_text(
      ggplot2::aes(label = tile_label),
      size = 4,
      lineheight = 0.95
    ) +
    ggplot2::scale_fill_viridis_c(
      option = "cividis",
      limits = c(0, 1),
      name = "Proportion"
    ) +
    ggplot2::labs(
      title = "Functional persistence architecture of network modules",
      subtitle = "Descriptive composition; conditional enrichment is tested downstream",
      x = NULL,
      y = NULL
    ) +
    theme_clean(12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))

  functional_plot_data <- data.table::copy(functional_composition)
  functional_plot_data[, module_label := factor(
    module_label,
    levels = sort(unique(module_label))
  )]
  functional_plot_data[, functional_category := factor(
    functional_category,
    levels = FUNCTIONAL_CATEGORY_ORDER
  )]
  functional_plot_data[, tile_label := ifelse(
    n_KOs == 0L,
    "",
    paste0(n_KOs, "\n", round(100 * proportion), "%")
  )]
  p_functional <- ggplot2::ggplot(
    functional_plot_data,
    ggplot2::aes(
      x = functional_category,
      y = module_label,
      fill = proportion
    )
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.55) +
    ggplot2::geom_text(
      ggplot2::aes(label = tile_label),
      size = 3.3,
      lineheight = 0.9
    ) +
    ggplot2::scale_fill_viridis_c(option = "cividis", name = "Proportion") +
    ggplot2::labs(
      title = "Functional programs represented within network modules",
      x = NULL,
      y = NULL
    ) +
    theme_clean(11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 38, hjust = 1)
    )

  overlap_plot_data <- data.table::copy(process_overlap)
  overlap_plot_data[, module_i := factor(
    module_i,
    levels = sort(unique(module_i))
  )]
  overlap_plot_data[, module_j := factor(
    module_j,
    levels = sort(unique(module_j))
  )]
  overlap_plot_data[, tile_label := paste0(
    shared_n,
    "\nJ=",
    sprintf("%.2f", jaccard)
  )]
  p_overlap <- ggplot2::ggplot(
    overlap_plot_data,
    ggplot2::aes(x = module_j, y = module_i, fill = jaccard)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.6) +
    ggplot2::geom_text(
      ggplot2::aes(label = tile_label),
      size = 3.8,
      lineheight = 0.95
    ) +
    ggplot2::scale_fill_viridis_c(
      option = "cividis",
      limits = c(0, 1),
      name = "Jaccard"
    ) +
    ggplot2::labs(
      title = "Overlap of KEGG MODULE repertoires",
      subtitle = "Cell values show shared modules and Jaccard similarity",
      x = NULL,
      y = NULL
    ) +
    theme_clean(12)

  nodes_plot <- data.table::copy(nodes)
  nodes_plot[, functional_category := factor(
    functional_category,
    levels = FUNCTIONAL_CATEGORY_ORDER
  )]
  edges_plot <- merge(
    edges[, .(KO1, KO2, mean_abs_weight)],
    nodes[, .(KO, x, y)],
    by.x = "KO1",
    by.y = "KO",
    all.x = TRUE,
    sort = FALSE
  )
  data.table::setnames(edges_plot, c("x", "y"), c("x1", "y1"))
  edges_plot <- merge(
    edges_plot,
    nodes[, .(KO, x, y)],
    by.x = "KO2",
    by.y = "KO",
    all.x = TRUE,
    sort = FALSE
  )
  data.table::setnames(edges_plot, c("x", "y"), c("x2", "y2"))
  module_centroids <- nodes[
    ,
    .(x = mean(x), y = mean(y)),
    by = .(module_label, module_program)
  ]
  p_network_category <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edges_plot,
      ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2),
      color = "#777777",
      alpha = 0.14,
      linewidth = 0.25
    ) +
    ggplot2::geom_point(
      data = nodes_plot,
      ggplot2::aes(
        x = x,
        y = y,
        color = functional_category,
        size = source_strength_abs
      ),
      alpha = 0.92
    ) +
    ggplot2::geom_label(
      data = module_centroids,
      ggplot2::aes(x = x, y = y, label = module_label),
      size = 3.6,
      fontface = "bold",
      fill = grDevices::adjustcolor("white", alpha.f = 0.82),
      label.size = 0.2
    ) +
    ggplot2::scale_color_manual(
      values = FUNCTIONAL_CATEGORY_COLORS,
      drop = FALSE,
      name = "Functional category"
    ) +
    ggplot2::scale_size_continuous(range = c(1.5, 5.5), guide = "none") +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = "Functional annotation of the conserved network backbone",
      subtitle = "Edges denote robust CLR co-abundance associations"
    ) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(color = "#555555")
    )

  hubs_figure <- top_hubs[hub_rank <= TOP_HUBS_FIGURE]
  hubs_figure[, function_label := paste0(
    KO,
    " - ",
    wrap_label(shorten_label(KO_function, 46L), 24L)
  )]
  hubs_figure[, function_label := stats::reorder(
    function_label,
    source_strength_abs
  )]
  p_hubs <- ggplot2::ggplot(
    hubs_figure,
    ggplot2::aes(
      x = source_strength_abs,
      y = function_label,
      color = functional_core_class
    )
  ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = source_strength_abs, yend = function_label),
      linewidth = 0.55,
      color = "#BDBDBD"
    ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::facet_wrap(~module_label, scales = "free_y", ncol = 2) +
    ggplot2::scale_color_manual(
      values = CLASS_COLORS,
      breaks = CLASS_LEVELS,
      drop = FALSE,
      name = "Prevalence class"
    ) +
    ggplot2::labs(
      title = "Highest-strength KOs within each network module",
      x = "Absolute node strength",
      y = NULL
    ) +
    theme_clean(10)

  save_gg(
    p_prevalence,
    file.path(FIGURES_DIR, "173_module_core_composition"),
    7.2,
    5.2,
    DPI
  )
  save_gg(
    p_functional,
    file.path(FIGURES_DIR, "173_module_functional_category_composition"),
    11.5,
    5.5,
    DPI
  )
  save_gg(
    p_overlap,
    file.path(FIGURES_DIR, "173_module_KEGG_process_overlap"),
    7,
    6,
    DPI
  )
  save_gg(
    p_network_category,
    file.path(FIGURES_DIR, "173_network_by_functional_category"),
    11,
    8,
    DPI
  )
  save_gg(
    p_hubs,
    file.path(FIGURES_DIR, "173_top_hubs_by_module"),
    12,
    14,
    DPI
  )

  cat("[9] Ejecutando validaciones estrictas...\n")
  class_counts <- table(factor(
    nodes$functional_core_class,
    levels = CLASS_LEVELS
  ))
  duplicate_nodes <- anyDuplicated(nodes$KO)
  edge_keys <- paste(
    pmin(edges$KO1, edges$KO2),
    pmax(edges$KO1, edges$KO2),
    sep = "||"
  )
  duplicate_edges <- anyDuplicated(edge_keys)
  missing_endpoints <- length(setdiff(
    unique(c(edges$KO1, edges$KO2)),
    nodes$KO
  ))
  missing_program_labels <- sum(is.na(nodes$module_program))
  missing_functional_categories <- sum(is.na(nodes$functional_category))
  unexpected_classes <- sum(
    !nodes$functional_core_class %in% CLASS_LEVELS
  )

  validation <- data.frame(
    metric = c(
      "n_network_nodes",
      "n_network_edges",
      "n_network_modules",
      "n_absolute_global_core_nodes",
      "n_near_global_nodes",
      "n_realm_core_only_nodes",
      "n_multi_realm_non_core_nodes",
      "n_realm_specific_nodes",
      "n_nodes_matched_to_reference_snapshot",
      "n_missing_functional_categories",
      "n_missing_module_program_labels",
      "n_unexpected_core_classes",
      "n_duplicate_nodes",
      "n_duplicate_edge_keys",
      "n_missing_edge_endpoints",
      "has_retained_KO_KEGG_links",
      "n_process_overlap_rows"
    ),
    observed = c(
      nrow(nodes),
      nrow(edges),
      data.table::uniqueN(nodes$module_label),
      unname(class_counts["absolute_global_core"]),
      unname(class_counts["near_global"]),
      unname(class_counts["realm_core_only"]),
      unname(class_counts["multi_realm_non_core"]),
      unname(class_counts["realm_specific"]),
      reference_match_n,
      missing_functional_categories,
      missing_program_labels,
      unexpected_classes,
      as.integer(duplicate_nodes > 0L),
      as.integer(duplicate_edges > 0L),
      missing_endpoints,
      as.integer(nrow(kegg_map) > 0L),
      nrow(process_overlap)
    ),
    expected = c(
      EXPECTED_NODES,
      EXPECTED_EDGES,
      EXPECTED_MODULES,
      EXPECTED_ABSOLUTE,
      EXPECTED_NEAR_GLOBAL,
      EXPECTED_REALM_CORE,
      EXPECTED_MULTI_REALM,
      EXPECTED_REALM_SPECIFIC,
      EXPECTED_NODES,
      0, 0, 0, 0, 0, 0,
      1,
      EXPECTED_MODULES^2
    ),
    stringsAsFactors = FALSE
  )
  validation$pass <- validation$observed == validation$expected
  write_tsv(
    validation,
    file.path(TABLES_DIR, "173_expected_validation.tsv")
  )

  if (STRICT_EXPECTED && any(!validation$pass)) {
    stop(
      "Fallo la validacion estricta: ",
      paste(validation$metric[!validation$pass], collapse = ", ")
    )
  }

  figure_bases <- c(
    "173_module_core_composition",
    "173_module_functional_category_composition",
    "173_module_KEGG_process_overlap",
    "173_network_by_functional_category",
    "173_top_hubs_by_module"
  )
  if (!MAKE_FIGURES) figure_bases <- character()
  figure_paths <- as.vector(outer(
    file.path(FIGURES_DIR, figure_bases),
    c(".pdf", ".png", ".tiff"),
    paste0
  ))
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
    file.path(TABLES_DIR, "173_figure_validation.tsv")
  )
  if (any(!figure_validation$pass)) stop("Fallo la escritura de figuras.")

  parameters <- data.frame(
    parameter = c(
      "run172",
      "run144_reference",
      "annotation_source",
      "external_KEGG_query",
      "old_core_classes_reused",
      "top_hubs_per_module",
      "top_kegg_per_module",
      "dpi"
    ),
    value = c(
      RUN172,
      RUN144_REFERENCE,
      "Frozen KO annotations and KO-to-KEGG MODULE links from script 144",
      "No",
      "No",
      TOP_HUBS_PER_MODULE,
      TOP_KEGG_PER_MODULE,
      DPI
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(parameters, file.path(TABLES_DIR, "173_parameters.tsv"))

  methods_text <- c(
    "Functional annotation of consensus network modules",
    "",
    paste0(
      "The major consensus-network modules were retained from scripts 171-172 ",
      "without reinference of edges or community membership."
    ),
    paste0(
      "KO functional descriptions, functional categories and KO-to-KEGG ",
      "MODULE links were recovered by KO identifier from the frozen annotation ",
      "snapshot previously generated by script 144. No new KEGGREST query was ",
      "performed."
    ),
    paste0(
      "Each KO was classified with the final hierarchical functional-persistence ",
      "definition established by script 153."
    ),
    paste0(
      "Module-level persistence-class and functional-category compositions are ",
      "descriptive; conditional enrichment and network null-model support must ",
      "be evaluated separately."
    ),
    paste0(
      "Network modules describe co-abundance programs of encoded functional ",
      "potential and not metabolic activity or direct biochemical interactions."
    )
  )
  writeLines(methods_text, file.path(TABLES_DIR, "173_methods_draft.txt"))

  summary_table <- data.frame(
    item = c(
      "status",
      "run172",
      "run144_reference",
      "n_nodes",
      "n_edges",
      "n_modules",
      "n_absolute_global_core",
      "n_near_global",
      "n_realm_core_only",
      "n_multi_realm_non_core",
      "n_realm_specific",
      "n_reference_matched_nodes",
      "n_KO_KEGG_links",
      "n_nodes_with_KEGG_module",
      "annotation_changed_by_external_query",
      "output_dir"
    ),
    value = c(
      "PASS",
      RUN172,
      RUN144_REFERENCE,
      nrow(nodes),
      nrow(edges),
      data.table::uniqueN(nodes$module_label),
      unname(class_counts["absolute_global_core"]),
      unname(class_counts["near_global"]),
      unname(class_counts["realm_core_only"]),
      unname(class_counts["multi_realm_non_core"]),
      unname(class_counts["realm_specific"]),
      reference_match_n,
      nrow(kegg_map),
      sum(nodes$n_KEGG_modules_snapshot > 0L),
      "No",
      OUTDIR
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(summary_table, file.path(TABLES_DIR, "173_summary.tsv"))

  saveRDS(
    list(
      nodes = as.data.frame(nodes),
      edges = as.data.frame(edges),
      module_summary = as.data.frame(module_summary),
      core_composition = as.data.frame(core_composition),
      functional_composition = as.data.frame(functional_composition),
      KO_to_KEGG_MODULE = as.data.frame(kegg_map),
      process_overlap = as.data.frame(process_overlap),
      parameters = parameters
    ),
    file.path(RDS_DIR, "173_annotated_network_module_objects.rds")
  )

  writeLines("PASS", STATUS_FILE)

  cat("\n============================================================\n")
  cat("STATUS: PASS\n")
  cat("Output dir: ", OUTDIR, "\n", sep = "")
  cat("Annotated nodes: ", nrow(nodes), "\n", sep = "")
  cat("Annotated edges: ", nrow(edges), "\n", sep = "")
  cat("Network modules: ", data.table::uniqueN(nodes$module_label), "\n", sep = "")
  cat(
    "Absolute/near-global/realm-core/multi-realm/realm-specific: ",
    paste(as.integer(class_counts), collapse = "/"), "\n",
    sep = ""
  )
  cat("KO-KEGG MODULE links retained: ", nrow(kegg_map), "\n", sep = "")
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
      file.path(LOGS_DIR, "173_error.txt")
    )
    1L
  }
)

quit(save = "no", status = status)
