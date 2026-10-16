#!/usr/bin/env Rscript

# ============================================================
# 177_PLOT_FUNCTIONAL_NETWORK_PANEL_A_EDITORIAL.R
#
# Figura principal de la sección:
#   Reproducible functional co-abundance structure
#
# Paneles:
#   a) Subred de los módulos mayores de la cohorte vigente
#   b) Especialización funcional validada
#   c) Estabilidad modular bajo perturbación de aristas
#
# Inputs directos: outputs finales de 171, 172, 173 y 174.
# No depende de RUN176 ni de una corrida anterior de 177.
# ============================================================

suppressPackageStartupMessages({
  required_packages <- c("data.table", "ggplot2", "patchwork")
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, quietly = TRUE,
            FUN.VALUE = logical(1))
  ]
  if (length(missing_packages) > 0L) {
    stop("Faltan paquetes requeridos: ", paste(missing_packages, collapse = ", "))
  }
})

# ------------------------------------------------------------
# Argumentos y rutas
# ------------------------------------------------------------

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
RESULTS_DIR <- file.path(BASE, "resultados")

read_latest_by_number <- function(number) {
  latest_files <- list.files(
    RESULTS_DIR,
    pattern = paste0("^LATEST_", number, "_.*\\.txt$"),
    full.names = TRUE
  )
  if (length(latest_files) == 0L) return(NA_character_)
  latest_files <- latest_files[
    order(file.info(latest_files)$mtime, decreasing = TRUE)
  ]
  for (latest_file in latest_files) {
    value <- readLines(latest_file, warn = FALSE, n = 1L)
    if (length(value) == 1L && dir.exists(trimws(value))) {
      return(trimws(value))
    }
  }
  NA_character_
}

RUN171 <- get_arg("--run171", read_latest_by_number("171"))
RUN172 <- get_arg("--run172", read_latest_by_number("172"))
RUN173 <- get_arg("--run173", read_latest_by_number("173"))
RUN174 <- get_arg("--run174", read_latest_by_number("174"))
STUDY_SCORES <- get_arg("--study_scores",file.path(BASE,"study_sensitivity","study_exclusion_module_scores.tsv"))
DPI <- suppressWarnings(as.integer(get_arg("--dpi", "600")))
STRICT_EXPECTED <- as_flag(get_arg("--strict_expected", "1"))

if (!is.finite(DPI) || DPI < 300L) stop("--dpi debe ser >= 300.")

SCRIPT_ID <- "177_PLOT_FUNCTIONAL_NETWORK_PANEL_A_EDITORIAL"
FIGURE_LAYOUT_VERSION <- "network_specialization_stability_current_cohort"
STAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUTDIR <- get_arg(
  "--outdir",
  file.path(RESULTS_DIR, paste0(SCRIPT_ID, "_", STAMP))
)
TABLES_DIR <- file.path(OUTDIR, "tables")
FIGURES_DIR <- file.path(OUTDIR, "figures")
LOGS_DIR <- file.path(OUTDIR, "logs")
RDS_DIR <- file.path(OUTDIR, "rds")

for (directory in c(TABLES_DIR, FIGURES_DIR, LOGS_DIR, RDS_DIR)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

STATUS_FILE <- file.path(OUTDIR, "STATUS.txt")
LOG_FILE <- file.path(LOGS_DIR, "177_log.txt")
ERROR_FILE <- file.path(LOGS_DIR, "177_error.txt")
LATEST_FILE <- file.path(RESULTS_DIR, paste0("LATEST_", SCRIPT_ID, ".txt"))

writeLines("RUNNING", STATUS_FILE)
writeLines(OUTDIR, LATEST_FILE)
writeLines(character(), LOG_FILE)

log_message <- function(...) {
  message_text <- paste0(...)
  cat(message_text, "\n", sep = "")
  cat(message_text, "\n", sep = "", file = LOG_FILE, append = TRUE)
}

assert_directory <- function(path, label) {
  if (is.na(path) || !nzchar(path) || !dir.exists(path)) {
    stop(label, " no existe: ", path)
  }
}

assert_file <- function(path) {
  if (!file.exists(path)) stop("No existe el archivo requerido: ", path)
}

read_tsv <- function(path) {
  assert_file(path)
  data.table::fread(path, check.names = FALSE, showProgress = FALSE,
                    na.strings = c("NA", "NaN", ""))
}

write_tsv <- function(x, path) {
  data.table::fwrite(x, path, sep = "\t", quote = FALSE, na = "NA")
}

normalize_names <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("[^a-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

pick_col <- function(x, patterns, label, required = TRUE) {
  original <- names(x)
  normalized <- normalize_names(original)
  for (pattern in patterns) {
    hit <- grep(pattern, normalized, perl = TRUE)
    if (length(hit) > 0L) return(original[hit[1L]])
  }
  if (required) {
    stop("No se identificó ", label, ". Columnas: ",
         paste(original, collapse = ", "))
  }
  NA_character_
}

safe_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

normalize_module <- function(x) {
  raw <- toupper(trimws(as.character(x)))
  out <- rep(NA_character_, length(raw))
  hit <- grepl("M[0-9]+", raw)
  out[hit] <- sprintf(
    "M%02d",
    as.integer(sub(".*M([0-9]+).*", "\\1", raw[hit]))
  )
  out
}

normalize_tier <- function(x) {
  raw <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(raw))
  out[raw %in% c("absolute_global_core", "absolute global core",
                 "absolute-global-core", "absolute_core")] <-
    "Absolute global core"
  out[raw %in% c("near_global", "near-global", "near global")] <-
    "Near-global"
  out[raw %in% c("realm_core_only", "realm-core", "realm core",
                 "realm_core")] <- "Realm-core"
  out[raw %in% c("multi_realm_non_core", "realm_specific",
                 "peripheral", "non-core", "non_core")] <- "Peripheral"
  out
}

save_plot_multi <- function(plot_object, file_base, width_mm, height_mm, dpi) {
  ggplot2::ggsave(
    paste0(file_base, ".pdf"), plot_object,
    width = width_mm, height = height_mm, units = "mm",
    device = grDevices::cairo_pdf, bg = "white"
  )
  ggplot2::ggsave(
    paste0(file_base, ".png"), plot_object,
    width = width_mm, height = height_mm, units = "mm",
    dpi = dpi, bg = "white"
  )
  ggplot2::ggsave(
    paste0(file_base, ".tiff"), plot_object,
    width = width_mm, height = height_mm, units = "mm",
    dpi = dpi, compression = "lzw", bg = "white"
  )
}

theme_editorial <- function(base_size = 9) {
  ggplot2::theme_classic(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 8, color = "black"),
      axis.title = ggplot2::element_text(size = 9, color = "black"),
      legend.text = ggplot2::element_text(size = 8, color = "black"),
      legend.title = ggplot2::element_text(size = 8, face = "bold"),
      plot.tag = ggplot2::element_text(size = 11, face = "bold"),
      plot.tag.position = "topleft",
      plot.title = ggplot2::element_text(size = 9, face = "bold", hjust = 0),
      plot.caption = ggplot2::element_text(
        size = 7.5, color = "#444444", hjust = 0, lineheight = 0.95,
        margin = ggplot2::margin(t = 3)
      ),
      axis.line = ggplot2::element_line(linewidth = 0.4),
      axis.ticks = ggplot2::element_line(linewidth = 0.35),
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(5, 8, 5, 6)
    )
}

theme_network <- function(base_size = 9) {
  ggplot2::theme_void(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      legend.text = ggplot2::element_text(size = 8, color = "black"),
      legend.title = ggplot2::element_text(size = 8, face = "bold"),
      plot.tag = ggplot2::element_text(size = 11, face = "bold"),
      plot.tag.position = "topleft",
      plot.title = ggplot2::element_text(size = 9, face = "bold", hjust = 0),
      plot.caption = ggplot2::element_text(
        size = 7.5, color = "#444444", hjust = 0, lineheight = 0.95,
        margin = ggplot2::margin(t = 3)
      ),
      plot.margin = ggplot2::margin(5, 6, 5, 6)
    )
}

MODULE_PALETTE <- c(
  "M01" = "#0072B2", "M02" = "#D55E00", "M03" = "#009E73",
  "M04" = "#CC79A7", "M05" = "#E69F00", "M06" = "#56B4E9"
)
TIER_ORDER <- c("Absolute global core", "Near-global", "Realm-core", "Peripheral")
TIER_SHAPES <- c(
  "Absolute global core" = 22,
  "Near-global" = 21,
  "Realm-core" = 24,
  "Peripheral" = 23
)

main <- function() {
  log_message("============================================================")
  log_message(SCRIPT_ID)
  log_message("Started: ", as.character(Sys.time()))
  log_message("RUN171: ", RUN171)
  log_message("RUN172: ", RUN172)
  log_message("RUN173: ", RUN173)
  log_message("RUN174: ", RUN174)
  log_message("FIGURE_LAYOUT_VERSION: ", FIGURE_LAYOUT_VERSION)
  log_message("DPI: ", DPI)
  log_message("STRICT_EXPECTED: ", STRICT_EXPECTED)
  log_message("OUTDIR: ", OUTDIR)
  log_message("============================================================")

  assert_directory(RUN171, "RUN171")
  assert_directory(RUN172, "RUN172")
  assert_directory(RUN173, "RUN173")
  assert_directory(RUN174, "RUN174")
  for (run in c(RUN171, RUN172, RUN173, RUN174)) {
    if (trimws(readLines(file.path(run, "STATUS.txt"), n = 1)) != "PASS") {
      stop("Input stage did not complete successfully: ", run)
    }
  }

  paths <- list(
    metrics = file.path(RUN171, "tables", "171_network_global_metrics.tsv"),
    nodes = file.path(RUN173, "tables", "173_network_nodes_annotated.tsv"),
    edges = file.path(RUN172, "tables", "172_network_edges_for_plot.tsv.gz"),
    enrichment = file.path(
      RUN174, "tables", "174_supported_module_functional_specialization.tsv"
    ),
    stability = file.path(
      RUN174, "tables", "174_module_membership_stability.tsv"
    ),
    perturbation = file.path(
      RUN174, "tables", "174_edge_perturbation_summary.tsv"
    ),
    validation_summary = file.path(
      RUN174, "tables", "174_summary.tsv"
    ),
    annotation_support = file.path(
      RUN174, "tables", "174_annotation_source_support.tsv"
    )
  )
  invisible(lapply(paths, assert_file))

  log_message("[1] Leyendo outputs finales de 171–174...")
  metrics_raw <- read_tsv(paths$metrics)
  nodes_raw <- read_tsv(paths$nodes)
  edges_raw <- read_tsv(paths$edges)
  enrichment_raw <- read_tsv(paths$enrichment)
  stability_raw <- read_tsv(paths$stability)
  perturbation_raw <- read_tsv(paths$perturbation)
  validation_summary_raw <- read_tsv(paths$validation_summary)
  annotation_support <- read_tsv(paths$annotation_support)

  # 174_edge_perturbation_summary.tsv contiene una fila por perturbación
  # (columna `perturbation`), no una columna repetida `n_perturbations`.
  if ("perturbation" %in% names(perturbation_raw)) {
    perturbation_ids <- safe_numeric(perturbation_raw$perturbation)
    perturbation_ids <- perturbation_ids[is.finite(perturbation_ids)]
    n_edge_perturbations <- data.table::uniqueN(perturbation_ids)
  } else {
    n_edge_perturbations <- nrow(perturbation_raw)
  }

  if (!is.finite(n_edge_perturbations) || n_edge_perturbations < 1L) {
    stop("No se pudo contar las perturbaciones de aristas de RUN174.")
  }

  n_edge_boot_summary <- NA_integer_
  if (all(c("item", "value") %in% names(validation_summary_raw))) {
    summary_hit <- validation_summary_raw[
      item %in% c("n_edge_boot", "n_edge_perturbations"), value
    ]
    summary_hit <- as.integer(safe_numeric(summary_hit))
    summary_hit <- summary_hit[is.finite(summary_hit)]
    if (length(summary_hit) > 0L) n_edge_boot_summary <- summary_hit[1L]
  }

  if (is.finite(n_edge_boot_summary) &&
      n_edge_boot_summary != n_edge_perturbations) {
    stop(
      "Inconsistencia en RUN174: la tabla contiene ",
      n_edge_perturbations, " perturbaciones, pero 174_summary.tsv reporta ",
      n_edge_boot_summary, "."
    )
  }

  # ----------------------------------------------------------
  # Red: estandarización y validación
  # ----------------------------------------------------------

  node_id_col <- pick_col(
    nodes_raw, c("^ko$", "^node$", "^node_id$", "^name$", "^id$"),
    "identificador KO"
  )
  x_col <- pick_col(nodes_raw, c("^x$", "^layout_x$", "^x_coord$"),
                    "coordenada x")
  y_col <- pick_col(nodes_raw, c("^y$", "^layout_y$", "^y_coord$"),
                    "coordenada y")
  module_col <- pick_col(
    nodes_raw,
    c("^module_label$", "^module_display$", "^network_module$", "^module$"),
    "módulo de red"
  )
  tier_col <- pick_col(
    nodes_raw,
    c("^functional_core_class$", "^functional_prevalence_class$",
      "^functional_core_class_recommended_minloc3$", "^core_class$",
      "^conservation_tier$"),
    "clase funcional de persistencia"
  )

  nodes <- unique(nodes_raw[, .(
    KO = trimws(as.character(get(node_id_col))),
    x = safe_numeric(get(x_col)),
    y = safe_numeric(get(y_col)),
    module = normalize_module(get(module_col)),
    tier = normalize_tier(get(tier_col))
  )], by = "KO")
  if (anyNA(nodes) || any(!is.finite(nodes$x)) || any(!is.finite(nodes$y))) {
    stop("Hay nodos sin coordenadas, módulo o clase de persistencia.")
  }
  MAJOR_MODULES <- sort(unique(nodes$module))
  if (length(MAJOR_MODULES) > length(MODULE_PALETTE)) {
    MODULE_PALETTE <- setNames(grDevices::hcl.colors(length(MAJOR_MODULES), "Dark 3"),
                               MAJOR_MODULES)
  } else {
    MODULE_PALETTE <- setNames(unname(MODULE_PALETTE[seq_along(MAJOR_MODULES)]),
                               MAJOR_MODULES)
  }
  nodes <- nodes[
    !is.na(KO) & nzchar(KO) & is.finite(x) & is.finite(y) &
      module %in% MAJOR_MODULES & tier %in% TIER_ORDER
  ]

  edge_from_col <- pick_col(
    edges_raw,
    c("^ko1$", "^from$", "^source$", "^node1$", "^v1$"),
    "origen de arista"
  )
  edge_to_col <- pick_col(
    edges_raw,
    c("^ko2$", "^to$", "^target$", "^node2$", "^v2$"),
    "destino de arista"
  )
  edges <- unique(edges_raw[, .(
    from = trimws(as.character(get(edge_from_col))),
    to = trimws(as.character(get(edge_to_col)))
  )])
  edges <- edges[
    from %in% nodes$KO & to %in% nodes$KO & from != to
  ]
  edges[, edge_key := paste(pmin(from, to), pmax(from, to), sep = "||")]
  edges <- edges[!duplicated(edge_key)]

  nodes[, x := as.numeric(scale(x))]
  nodes[, y := as.numeric(scale(y))]
  from_xy <- nodes[, .(from = KO, x, y)]
  to_xy <- nodes[, .(to = KO, xend = x, yend = y)]
  edges_plot <- merge(edges, from_xy, by = "from", all = FALSE)
  edges_plot <- merge(edges_plot, to_xy, by = "to", all = FALSE)

  nodes[, module := factor(module, levels = MAJOR_MODULES)]
  nodes[, tier := factor(tier, levels = TIER_ORDER)]

  module_counts <- nodes[, .(n_KOs = .N), by = module]
  data.table::setorder(module_counts, module)
  tier_counts <- nodes[, .(n_KOs = .N), by = tier]

  observed_modules <- setNames(integer(length(MAJOR_MODULES)), MAJOR_MODULES)
  observed_modules[as.character(module_counts$module)] <- module_counts$n_KOs
  observed_tiers <- setNames(integer(length(TIER_ORDER)), TIER_ORDER)
  observed_tiers[as.character(tier_counts$tier)] <- tier_counts$n_KOs

  metric_value <- function(metric_name) {
    metric_col <- pick_col(metrics_raw, c("^metric$", "^item$"),
                           "nombre de métrica")
    value_col <- pick_col(metrics_raw, c("^value$", "^observed$"),
                          "valor de métrica")
    value <- metrics_raw[
      as.character(get(metric_col)) == metric_name,
      safe_numeric(get(value_col))
    ]
    if (length(value) == 0L) return(NA_real_)
    value[1L]
  }
  full_nodes <- metric_value("n_connected_KOs")
  full_edges <- metric_value("n_consensus_edges")
  # Validate against the current inference, never historical manuscript counts.
  primary_edges <- read_tsv(file.path(RUN171, "tables", "171_consensus_edges.tsv.gz"))
  expected_edges <- primary_edges[KO1 %in% nodes$KO & KO2 %in% nodes$KO]
  validation <- data.table::data.table(
    metric = c("major_module_KOs", "major_modules", "major_subgraph_edges",
               "edge_endpoints_match", "module_sizes_match_stability", "unique_KOs"),
    pass = c(nrow(nodes) == metric_value("n_KOs_in_major_modules"),
             length(MAJOR_MODULES) == metric_value("n_major_modules_ge_min_size"),
             nrow(edges_plot) == nrow(expected_edges),
             setequal(edges$edge_key, paste(pmin(expected_edges$KO1, expected_edges$KO2),
                                           pmax(expected_edges$KO1, expected_edges$KO2), sep = "||")),
             all(observed_modules[stability_raw$module_label] == stability_raw$n_KOs),
             nrow(nodes) == data.table::uniqueN(nodes$KO))
  )
  write_tsv(validation, file.path(TABLES_DIR, "177_expected_validation.tsv"))
  if (anyNA(validation$pass) || !all(validation$pass)) stop("Figure inputs do not match current inference.")

  # ----------------------------------------------------------
  # Panel a: subred funcional
  # ----------------------------------------------------------

  log_message("[2] Construyendo panel a: subred funcional...")
  x_range <- diff(range(nodes$x))
  y_range <- diff(range(nodes$y))
  offsets <- data.table::data.table(
    module = factor(MAJOR_MODULES, levels = MAJOR_MODULES),
    dx = rep_len(c(0.05, 0.10, -0.30, 0.38, 0.25, -0.08), length(MAJOR_MODULES)),
    dy = rep_len(c(0.30, -0.35, 0.17, 0.38, -0.22, -0.23), length(MAJOR_MODULES))
  )
  module_centroids <- nodes[, .(
    x = median(x), y = median(y), n_KOs = .N
  ), by = module]
  module_centroids <- merge(module_centroids, offsets, by = "module",
                            all.x = TRUE, sort = FALSE)
  # A functional name is assigned only to a supported category. Module IDs are
  # size ranks, so names are rebuilt from the current enrichment output.
  module_names <- data.table::data.table(module = MAJOR_MODULES, functional_name = "Diverse functions")
  for (m in MAJOR_MODULES) {
    hit <- annotation_support[module_label == m & primary_supported == TRUE]
    if (nrow(hit)) {
      category <- sub(" systems$", "", hit$category[1L])
      category <- sub("Carbon / carbohydrate metabolism", "Carbon / carbohydrate\nmetabolism", category, fixed=TRUE)
      category <- sub("Amino acid / cofactor biosynthesis", "Amino acid / cofactor\nbiosynthesis", category, fixed = TRUE)

      if (!isTRUE(hit$combined_supported[1L])) category <- paste0(category, "*")
      module_names[module == m, functional_name := category]
    }
  }
  module_centroids <- merge(module_centroids, module_names, by = "module", all.x = TRUE, sort = FALSE)
  module_centroids[, `:=`(
    x_label = x + dx * x_range,
    y_label = y + dy * y_range,
    label = paste0(module, " · ", n_KOs, " KOs\n", functional_name)
  )]

  p_a <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edges_plot,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      color = "#7A7A7A", linewidth = 0.18, alpha = 0.28
    ) +
    ggplot2::geom_point(
      data = nodes,
      ggplot2::aes(x = x, y = y, fill = module, shape = tier),
      color = "#222222", size = 1.55, stroke = 0.22, alpha = 0.96
    ) +
    ggplot2::geom_segment(
      data = module_centroids,
      ggplot2::aes(x = x, y = y, xend = x_label, yend = y_label, colour = module),
      linewidth = 0.35, show.legend = FALSE
    ) +
    ggplot2::geom_label(
      data = module_centroids,
      ggplot2::aes(x = x_label, y = y_label, label = label),
      colour = "black", fill = "white", size = 2.9, fontface = "bold", label.size = 0,
      label.padding = grid::unit(0.8, "mm"), show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(values = MODULE_PALETTE, guide = "none") +
    ggplot2::scale_colour_manual(values = MODULE_PALETTE, guide = "none") +
    ggplot2::scale_shape_manual(values = TIER_SHAPES,
                                name = "Persistence class") +
    ggplot2::coord_equal(clip = "off") +
    ggplot2::labs(tag = "a", title = "Co-abundance network") +
    theme_network() +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.justification = "left"
    ) +
    ggplot2::guides(
      shape = ggplot2::guide_legend(
        order = 2, nrow = 1, byrow = TRUE,
        title.position = "left", title.hjust = 0,
        override.aes = list(fill = "white", size = 3.2, alpha = 1)
      )
    )

  # ----------------------------------------------------------
  # Panel b: especialización funcional validada
  # ----------------------------------------------------------

  log_message("[3] Construyendo panel b: especialización funcional...")
  required_enrichment <- c(
    "module_label", "category", "observed_n", "fold_enrichment",
    "q_BY", "p_maxT_FWER"
  )
  missing_enrichment <- setdiff(required_enrichment, names(enrichment_raw))
  if (length(missing_enrichment) > 0L) {
    stop("Faltan columnas de enriquecimiento: ",
         paste(missing_enrichment, collapse = ", "))
  }

  enrichment <- enrichment_raw[, .(
    module = normalize_module(module_label),
    category = as.character(category),
    observed_n = as.integer(safe_numeric(observed_n)),
    fold_enrichment = safe_numeric(fold_enrichment),
    q_BY = safe_numeric(q_BY),
    p_maxT_FWER = safe_numeric(p_maxT_FWER)
  )]
  enrichment <- enrichment[
    module %in% MAJOR_MODULES &
      q_BY <= 0.05 & p_maxT_FWER <= 0.05
  ]
  enrichment[, category_short := vapply(category, function(z) {
    paste(strwrap(z, width = 24), collapse = "\n")
  }, character(1))]
  enrichment <- merge(enrichment,
    annotation_support[, .(module = module_label, category, combined_supported)],
    by = c("module", "category"), all.x = TRUE, sort = FALSE)
  enrichment[is.na(combined_supported), combined_supported := FALSE]
  enrichment[, row_label := module]
  data.table::setorder(enrichment, module, category)
  enrichment[, row_label := factor(row_label, levels = rev(row_label))]
  enrichment[, module := factor(module, levels = MAJOR_MODULES)]
  enrichment[, point_label := sprintf("%.2f×; %d KOs",
                                       fold_enrichment, observed_n)]


  enrichment[combined_supported == FALSE, point_label := paste0(point_label, "*")]

  if (anyNA(enrichment$module)) stop("Especialización sin módulo válido.")
  write_tsv(enrichment,
            file.path(TABLES_DIR, "177_validated_functional_enrichment_used.tsv"))

  p_b <- ggplot2::ggplot(
    enrichment,
    ggplot2::aes(x = fold_enrichment, y = row_label, color = module)
  ) +
    ggplot2::geom_vline(xintercept = 1, color = "#777777", linewidth = 0.45,
                        linetype = 2) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 1, xend = fold_enrichment,
                   y = row_label, yend = row_label),
      color = "#C7C7C7", linewidth = 0.80
    ) +
    ggplot2::geom_point(
      data = enrichment[combined_supported == TRUE],
      ggplot2::aes(size = observed_n, fill = module), shape = 21, stroke = 0.55,
      alpha = 0.97, show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data=enrichment[combined_supported==FALSE],
      ggplot2::aes(size=observed_n),shape=21,stroke=0.75,fill="white",show.legend=FALSE
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = point_label), hjust = 0, nudge_x = 0.24,
      color = "black", size = 2.75
    ) +
    ggplot2::scale_color_manual(values = MODULE_PALETTE, guide = "none") +
    ggplot2::scale_fill_manual(values = MODULE_PALETTE, guide = "none") +
    ggplot2::scale_size_continuous(
      range = c(3.2, 4.3), guide = "none"
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0.95, max(c(4.55, enrichment$fold_enrichment * 1.45))),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
      x = "Fold enrichment", y = NULL, tag = "b", title = "Functional enrichment",
      caption = "* Support varies with annotation source"
    ) +
    theme_editorial() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 8),
      plot.margin = ggplot2::margin(5, 14, 5, 6)
    )
  if (nrow(enrichment) == 0L) {
    p_b <- ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0,
      label = "No category passed both adjusted thresholds", size = 3) +
      ggplot2::labs(tag = "b") + ggplot2::theme_void()
  }

  # ----------------------------------------------------------
  # Panel c: estabilidad de módulos
  # ----------------------------------------------------------

  log_message("[4] Construyendo panel c: estabilidad modular...")
  required_stability <- c("module_label", "n_KOs", "mean_stability")
  missing_stability <- setdiff(required_stability, names(stability_raw))
  if (length(missing_stability) > 0L) {
    stop("Faltan columnas de estabilidad: ",
         paste(missing_stability, collapse = ", "))
  }
  stability <- stability_raw[, .(
    module = normalize_module(module_label),
    n_KOs = as.integer(safe_numeric(n_KOs)),
    mean_stability = safe_numeric(mean_stability)
  )]
  stability <- stability[module %in% MAJOR_MODULES]
  stability[, module := factor(
    module,
    levels = stability[order(mean_stability), as.character(module)]
  )]
  stability[, value_label := sprintf("%.3f", mean_stability)]
  if (!setequal(as.character(stability$module), MAJOR_MODULES)) {
    stop("La tabla de estabilidad no cubre todos los módulos mayores.")
  }
  write_tsv(stability,
            file.path(TABLES_DIR, "177_module_membership_stability_used.tsv"))

  study_scores <- read_tsv(STUDY_SCORES)
  study_scores <- study_scores[grepl("^PRJ",study_id)]
  stopifnot(data.table::uniqueN(study_scores$study_id)==43L)
  study_stats <- study_scores[,.(study_median=median(within_module_coassignment),
    study_min=min(within_module_coassignment),study_max=max(within_module_coassignment),
    n_study_units=data.table::uniqueN(study_id)),by=.(module=module_label)]
  stability[,module:=as.character(module)]
  stability <- merge(stability,study_stats,by="module",all.x=TRUE,sort=FALSE)
  stability[,y:=match(module,rev(MAJOR_MODULES))]
  stopifnot(!anyNA(stability))
  write_tsv(stability,file.path(TABLES_DIR,"177_module_membership_stability_used.tsv"))
  p_c <- ggplot2::ggplot(stability,ggplot2::aes(colour=module)) +
    ggplot2::geom_segment(ggplot2::aes(x=study_min,xend=study_max,y=y-0.12,yend=y-0.12),
      linewidth=0.55,alpha=0.60) +
    ggplot2::geom_point(ggplot2::aes(x=mean_stability,y=y+0.12),shape=16,size=2.7) +
    ggplot2::geom_point(ggplot2::aes(x=study_median,y=y-0.12),shape=23,fill="white",size=2.6,stroke=0.65) +
    ggplot2::scale_color_manual(values=MODULE_PALETTE,guide="none") +
    ggplot2::scale_x_continuous(limits=c(0,1.035),breaks=c(0,.25,.5,.75,1),
      expand=ggplot2::expansion(mult=c(0,0))) +
    ggplot2::scale_y_continuous(breaks=seq_along(MAJOR_MODULES),labels=rev(MAJOR_MODULES),
      limits=c(.55,length(MAJOR_MODULES)+.45),expand=ggplot2::expansion(mult=c(0,0))) +
    ggplot2::labs(x="Within-module co-assignment",y=NULL,tag="c",title="Membership sensitivity",
      subtitle="Circles: edge perturbation\nDiamonds: BioProject omission (range)") +
    theme_editorial() +
    ggplot2::theme(plot.margin=ggplot2::margin(5,14,5,6),
                   plot.subtitle=ggplot2::element_text(size=7.6,colour="black",lineheight=1.05))

  # ----------------------------------------------------------
  # Ensamblaje y exportación
  # ----------------------------------------------------------

  log_message("[5] Ensamblando figura principal y exportando...")
  design <- "AB\nAC"
  composite <- p_a + p_b + p_c +
    patchwork::plot_layout(
      design = design,
      widths = c(1.55, 1),
      heights = c(1, 1),
      guides = "collect"
    ) &
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.justification = "center"
    )

  figure_specs <- list(
    list(composite, "177_reproducible_functional_coabundance_structure", 190, 125),
    list(p_a, "177_panel_a_functional_coabundance_network", 135, 140),
    list(p_b, "177_panel_b_validated_functional_specialization", 110, 65),
    list(p_c, "177_panel_c_module_stability", 110, 70)
  )
  for (spec in figure_specs) {
    save_plot_multi(
      spec[[1L]], file.path(FIGURES_DIR, spec[[2L]]),
      width_mm = spec[[3L]], height_mm = spec[[4L]], dpi = DPI
    )
  }

  write_tsv(nodes, file.path(TABLES_DIR, "177_network_nodes_used.tsv"))
  data.table::fwrite(
    edges_plot,
    file.path(TABLES_DIR, "177_network_edges_used.tsv.gz"),
    sep = "\t", quote = FALSE, na = "NA"
  )
  write_tsv(module_centroids,
            file.path(TABLES_DIR, "177_network_module_label_positions.tsv"))

  panel_inventory <- data.table::data.table(
    panel = c("a", "b", "c"),
    content = c(
      "Major-module functional co-abundance network",
      "Validated functional specialization",
      "Membership sensitivity to edge perturbation and BioProject omission"
    ),
    source_run = c("RUN172/RUN173", "RUN174", "RUN174"),
    included_in_main_figure = TRUE
  )
  write_tsv(
    panel_inventory,
    file.path(TABLES_DIR, "177_main_figure_panel_inventory.tsv")
  )

  summary_table <- data.table::data.table(
    item = c(
      "status", "figure_layout_version",
      "run171", "run172", "run173", "run174",
      "n_plotted_major_module_KOs",
      "n_plotted_major_module_associations",
      "n_full_consensus_KOs", "n_full_consensus_associations",
      "n_major_modules", "n_supported_specializations",
      "n_edge_perturbations", "n_figure_panels",
      "n_figure_files", "output_directory"
    ),
    value = as.character(c(
      "PASS", FIGURE_LAYOUT_VERSION,
      RUN171, RUN172, RUN173, RUN174,
      nrow(nodes), nrow(edges_plot), as.integer(full_nodes),
      as.integer(full_edges), length(unique(nodes$module)), nrow(enrichment),
      n_edge_perturbations, 3L, length(figure_specs) * 3L, OUTDIR
    ))
  )
  write_tsv(summary_table, file.path(TABLES_DIR, "177_summary.tsv"))

  saveRDS(
    list(
      paths = paths,
      nodes = nodes,
      edges = edges_plot,
      enrichment = enrichment,
      stability = stability,
      module_centroids = module_centroids,
      plots = list(a = p_a, b = p_b, c = p_c,
                   composite = composite)
    ),
    file.path(RDS_DIR, "177_figure_objects.rds")
  )

  writeLines("PASS", STATUS_FILE)
  log_message("============================================================")
  log_message("STATUS: PASS")
  log_message("Output dir: ", OUTDIR)
  log_message("Plotted KOs: ", nrow(nodes))
  log_message("Plotted associations: ", nrow(edges_plot))
  log_message("Supported specializations: ", nrow(enrichment))
  log_message("Edge perturbations: ", n_edge_perturbations)
  log_message("Figure files: ", length(figure_specs) * 3L)
  log_message("Finished: ", as.character(Sys.time()))
  log_message("============================================================")
}

status <- tryCatch(
  {
    main()
    "PASS"
  },
  error = function(e) {
    writeLines("FAIL", STATUS_FILE)
    writeLines(conditionMessage(e), ERROR_FILE)
    cat("ERROR: ", conditionMessage(e), "\n", sep = "")
    cat("ERROR: ", conditionMessage(e), "\n", sep = "",
        file = LOG_FILE, append = TRUE)
    "FAIL"
  }
)

if (status != "PASS") quit(save = "no", status = 1L)
