#!/usr/bin/env Rscript

# ============================================================
# 154_FIG4B_FUNCTIONAL_PCA_LOADINGS_RECONCILED.R
#
# Objetivo
#   Generar el panel b de la Figura 4 mediante:
#   1) PCA global de perfiles KO transformados por CLR.
#   2) Inset con la contribucion ponderada por clase funcional
#      entre los 100 loadings absolutos mas altos de PC1.
#
# Contratos metodologicos
#   - Llave primaria: run_accession.
#   - Matriz funcional canónica: MATRIX_KEGG_ko.tsv.
#   - Cohorte funcional definida por el run 153: 194 muestras.
#   - No se aplica balanceo geografico ni submuestreo de Guadalupe.
#   - Clasificacion KO reconciliada del run 153.
#   - CLR con pseudoconteo 1; PCA centrada sin reescalado.
#   - Sin elipses, centroides, titulo ni letra de panel.
#   - Salida de 87 x 100 mm y tipografia de 8-11 pt.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# -----------------------------
# Argumentos
# -----------------------------
get_arg <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- which(args == flag)
  if (length(hit) == 0L || hit[length(hit)] >= length(args)) return(default)
  args[hit[length(hit)] + 1L]
}

BASE <- "/home/fjbalvino/Biogeografia_2026"
RESULTS_DIR <- file.path(BASE, "resultados")

read_latest <- function(prefix) {
  path <- file.path(RESULTS_DIR, paste0("LATEST_", prefix, ".txt"))
  if (!file.exists(path)) return(NA_character_)
  value <- readLines(path, warn = FALSE)
  if (length(value) == 0L) return(NA_character_)
  trimws(value[1L])
}

SCRIPT_ID <- "154_FIG4B_FUNCTIONAL_PCA_LOADINGS_RECONCILED"

RUN153 <- get_arg(
  "--run153",
  read_latest("153_FIG4A_FUNCTIONAL_KO_CLASSES_RECONCILED")
)

KO_MATRIX <- get_arg(
  "--ko_matrix",
  "/data/Magenta/results/035.metagenome_assembly_annotation/final_results/MATRIX_KEGG_ko.tsv"
)

TOP_N <- as.integer(get_arg("--top_n", "100"))
PSEUDOCOUNT <- as.numeric(get_arg("--pseudocount", "1"))
EXPECTED_SAMPLES <- as.integer(get_arg("--expected_samples", "194"))
EXPECTED_KOS <- as.integer(get_arg("--expected_kos", "17508"))
STRICT_EXPECTED <- as.integer(get_arg("--strict_expected", "1")) == 1L
ORIENT_AXES <- as.integer(get_arg("--orient_axes", "1")) == 1L
PANEL_WIDTH_MM <- as.numeric(get_arg("--panel_width_mm", "87"))
PANEL_HEIGHT_MM <- as.numeric(get_arg("--panel_height_mm", "100"))
BASE_FONT_SIZE <- as.numeric(get_arg("--base_font_size", "8"))
TEXT_SIZE_MM <- BASE_FONT_SIZE / 2.845276
DPI <- as.integer(get_arg("--dpi", "600"))

if (!is.finite(TOP_N) || TOP_N < 1L) stop("top_n debe ser >= 1.")
if (!is.finite(PSEUDOCOUNT) || PSEUDOCOUNT <= 0) {
  stop("pseudocount debe ser > 0.")
}
if (!is.finite(BASE_FONT_SIZE) || BASE_FONT_SIZE < 8 || BASE_FONT_SIZE > 11) {
  stop("base_font_size debe estar entre 8 y 11 pt.")
}
if (!is.finite(PANEL_WIDTH_MM) || PANEL_WIDTH_MM <= 0 ||
    !is.finite(PANEL_HEIGHT_MM) || PANEL_HEIGHT_MM <= 0) {
  stop("Las dimensiones del panel deben ser positivas.")
}

# -----------------------------
# Directorios y log
# -----------------------------
STAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUTDIR <- file.path(RESULTS_DIR, paste0(SCRIPT_ID, "_", STAMP))
FIG_DIR <- file.path(OUTDIR, "figures")
TAB_DIR <- file.path(OUTDIR, "tables")
RDS_DIR <- file.path(OUTDIR, "rds")
LOG_DIR <- file.path(OUTDIR, "logs")

for (d in c(OUTDIR, FIG_DIR, TAB_DIR, RDS_DIR, LOG_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

LATEST_FILE <- file.path(RESULTS_DIR, paste0("LATEST_", SCRIPT_ID, ".txt"))
STATUS_FILE <- file.path(OUTDIR, "STATUS.txt")
LOG_FILE <- file.path(LOG_DIR, "154_run.log")

writeLines(OUTDIR, LATEST_FILE)
writeLines("RUNNING", STATUS_FILE)

log_msg <- function(...) {
  msg <- paste0(..., collapse = "")
  line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg)
  cat(line, "\n", sep = "")
  cat(line, "\n", file = LOG_FILE, append = TRUE, sep = "")
}

write_tsv <- function(x, path) {
  fwrite(as.data.table(x), path, sep = "\t", quote = FALSE, na = "NA")
}

write_status <- function(status, extra = character()) {
  writeLines(
    c(
      status,
      paste0("timestamp=", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      paste0("output_directory=", OUTDIR),
      extra
    ),
    STATUS_FILE
  )
}

# -----------------------------
# Utilidades de IDs y matriz KO
# -----------------------------
clean_chr <- function(x) {
  y <- trimws(as.character(x))
  y[y == "" | tolower(y) %in% c("na", "n/a", "nan", "null", "none")] <- NA_character_
  y
}

norm_id <- function(x) {
  y <- tolower(clean_chr(x))
  gsub("[^a-z0-9]", "", y)
}

match_ids_strict <- function(input_ids, reference_ids) {
  input_ids <- clean_chr(input_ids)
  reference_ids <- clean_chr(reference_ids)

  if (anyDuplicated(reference_ids[!is.na(reference_ids)])) {
    stop("La referencia contiene run_accession duplicados.")
  }
  if (anyDuplicated(norm_id(reference_ids)[!is.na(reference_ids)])) {
    stop("Los run_accession colisionan despues de normalizar IDs.")
  }

  out <- data.table(
    input_id = input_ids,
    matched_id = NA_character_,
    match_mode = NA_character_
  )

  idx <- match(input_ids, reference_ids)
  hit <- !is.na(idx)
  out$matched_id[hit] <- reference_ids[idx[hit]]
  out$match_mode[hit] <- "exact"

  remaining <- which(is.na(out$matched_id) & !is.na(input_ids))
  if (length(remaining) > 0L) {
    idx <- match(tolower(input_ids[remaining]), tolower(reference_ids))
    hit <- !is.na(idx)
    out$matched_id[remaining[hit]] <- reference_ids[idx[hit]]
    out$match_mode[remaining[hit]] <- "case_insensitive"
  }

  remaining <- which(is.na(out$matched_id) & !is.na(input_ids))
  if (length(remaining) > 0L) {
    idx <- match(norm_id(input_ids[remaining]), norm_id(reference_ids))
    hit <- !is.na(idx)
    out$matched_id[remaining[hit]] <- reference_ids[idx[hit]]
    out$match_mode[remaining[hit]] <- "normalized"
  }

  out
}

numeric_matrix_from_df <- function(df) {
  out <- as.matrix(data.frame(
    lapply(df, function(x) suppressWarnings(as.numeric(as.character(x)))),
    check.names = FALSE
  ))
  out[is.na(out)] <- 0
  storage.mode(out) <- "numeric"
  out
}

load_ko_matrix <- function(path, reference_samples) {
  if (!file.exists(path)) stop("No existe la matriz KO: ", path)

  raw <- as.data.frame(
    fread(path, data.table = FALSE, check.names = FALSE),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  ko_hits <- vapply(
    raw,
    function(x) sum(grepl("^K[0-9]{5}$", clean_chr(x))),
    integer(1)
  )

  if (max(ko_hits) < 100L) {
    stop("No fue posible identificar una columna de KOs en la matriz.")
  }

  ko_col <- names(ko_hits)[which.max(ko_hits)]
  ko_ids <- clean_chr(raw[[ko_col]])
  keep_ko <- grepl("^K[0-9]{5}$", ko_ids)
  raw <- raw[keep_ko, , drop = FALSE]
  ko_ids <- ko_ids[keep_ko]

  candidate_cols <- setdiff(names(raw), ko_col)
  matching <- match_ids_strict(candidate_cols, reference_samples)
  matching <- matching[!is.na(matched_id)]

  if (nrow(matching) == 0L) {
    stop("No se encontraron columnas asociables a run_accession.")
  }
  if (anyDuplicated(matching$matched_id)) {
    stop("Varias columnas KO se asociaron al mismo run_accession.")
  }

  mat <- numeric_matrix_from_df(raw[, matching$input_id, drop = FALSE])
  colnames(mat) <- matching$matched_id

  if (anyDuplicated(ko_ids)) {
    mat <- rowsum(mat, group = ko_ids, reorder = FALSE)
  } else {
    rownames(mat) <- ko_ids
  }

  list(matrix = mat, matching = matching, ko_column = ko_col)
}

# -----------------------------
# Paletas y tema
# -----------------------------
REALM_ORDER <- c(
  "Central Indo-Pacific",
  "Indo-Pacific",
  "Temperate Northern Atlantic",
  "Tropical Atlantic",
  "Tropical Eastern Pacific",
  "Western Indo-Pacific"
)

REALM_PALETTE <- c(
  "Central Indo-Pacific" = "#66C2A5",
  "Indo-Pacific" = "#FC8D62",
  "Temperate Northern Atlantic" = "#8DA0CB",
  "Tropical Atlantic" = "#E78AC3",
  "Tropical Eastern Pacific" = "#A6D854",
  "Western Indo-Pacific" = "#FFD92F"
)

REALM_LABELS <- c(
  "Central Indo-Pacific" = "Central Indo-Pacific",
  "Indo-Pacific" = "Indo-Pacific",
  "Temperate Northern Atlantic" = "Temperate N. Atlantic",
  "Tropical Atlantic" = "Tropical Atlantic",
  "Tropical Eastern Pacific" = "Tropical E. Pacific",
  "Western Indo-Pacific" = "Western Indo-Pacific"
)

CLASS_ORDER <- c(
  "absolute_global_core",
  "near_global",
  "realm_core_only",
  "functional_periphery"
)

CLASS_PALETTE <- c(
  "absolute_global_core" = "#7F0000",
  "near_global" = "#C62828",
  "realm_core_only" = "#EF5350",
  "functional_periphery" = "#7F7F7F"
)

CLASS_LABELS <- c(
  "absolute_global_core" = "Absolute",
  "near_global" = "Near-global",
  "realm_core_only" = "Realm-core",
  "functional_periphery" = "Periphery"
)

theme_panel <- function(base_size = BASE_FONT_SIZE) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      text = element_text(color = "black"),
      axis.text = element_text(size = base_size, color = "black"),
      axis.title = element_text(size = base_size, color = "black"),
      axis.line = element_line(linewidth = 0.42, color = "black"),
      axis.ticks = element_line(linewidth = 0.35, color = "black"),
      legend.text = element_text(size = base_size, color = "black"),
      legend.key.height = grid::unit(3.0, "mm"),
      legend.key.width = grid::unit(3.2, "mm"),
      legend.spacing.x = grid::unit(0.8, "mm"),
      legend.spacing.y = grid::unit(0.2, "mm"),
      plot.margin = margin(3, 4, 2, 3)
    )
}

# -----------------------------
# Analisis
# -----------------------------
main <- function() {
  log_msg("Starting ", SCRIPT_ID)
  log_msg("Run 153: ", RUN153)
  log_msg("KO matrix: ", KO_MATRIX)
  log_msg("Output: ", OUTDIR)

  if (is.na(RUN153) || !nzchar(RUN153) || !dir.exists(RUN153)) {
    stop("run153 no es valido: ", RUN153)
  }

  object_file <- file.path(
    RUN153,
    "rds",
    "153_Figure_4a_analysis_objects.rds"
  )
  if (!file.exists(object_file)) stop("No existe: ", object_file)

  obj153 <- readRDS(object_file)
  required_objects <- c("primary_annotation", "metadata_aligned", "panel_summary")
  missing_objects <- setdiff(required_objects, names(obj153))
  if (length(missing_objects) > 0L) {
    stop("El objeto 153 no contiene: ", paste(missing_objects, collapse = ", "))
  }

  metadata <- as.data.table(obj153$metadata_aligned)
  annotation <- as.data.table(obj153$primary_annotation)

  required_meta <- c("run_accession", "realm", "locality")
  required_ann <- c("KO", "functional_core_class")
  if (length(setdiff(required_meta, names(metadata))) > 0L) {
    stop("metadata_aligned no contiene las columnas requeridas.")
  }
  if (length(setdiff(required_ann, names(annotation))) > 0L) {
    stop("primary_annotation no contiene KO y functional_core_class.")
  }

  metadata[, run_accession := clean_chr(run_accession)]
  metadata[, realm := clean_chr(realm)]
  annotation[, KO := clean_chr(KO)]

  if (anyNA(metadata$run_accession) || anyDuplicated(metadata$run_accession)) {
    stop("run_accession faltantes o duplicados en metadata 153.")
  }
  if (anyNA(metadata$realm)) stop("Hay muestras funcionales sin marine realm.")
  if (anyNA(annotation$KO) || anyDuplicated(annotation$KO)) {
    stop("KO faltantes o duplicados en la clasificacion 153.")
  }

  unknown_realms <- setdiff(unique(metadata$realm), REALM_ORDER)
  if (length(unknown_realms) > 0L) {
    stop("Marine realms sin color: ", paste(unknown_realms, collapse = ", "))
  }

  loaded <- load_ko_matrix(KO_MATRIX, metadata$run_accession)
  mat <- loaded$matrix
  write_tsv(loaded$matching, file.path(TAB_DIR, "154_sample_matching.tsv"))

  sample_idx <- match(metadata$run_accession, colnames(mat))
  if (anyNA(sample_idx)) {
    stop("Faltan muestras del run153 en la matriz KO.")
  }
  mat <- mat[, sample_idx, drop = FALSE]
  if (!identical(colnames(mat), metadata$run_accession)) {
    stop("Fallo de alineacion por run_accession.")
  }

  ko_idx <- match(annotation$KO, rownames(mat))
  if (anyNA(ko_idx)) {
    stop("Faltan KOs clasificados por 153 en la matriz funcional.")
  }
  mat <- mat[ko_idx, , drop = FALSE]
  if (!identical(rownames(mat), annotation$KO)) {
    stop("Fallo de alineacion entre matriz KO y clasificacion 153.")
  }

  keep_nonzero <- rowSums(mat > 0) > 0L
  if (!all(keep_nonzero)) {
    stop("La clasificacion 153 contiene KOs completamente nulos.")
  }

  n_samples <- ncol(mat)
  n_kos <- nrow(mat)
  log_msg("Functional samples: ", n_samples)
  log_msg("Non-zero KOs: ", n_kos)

  if (STRICT_EXPECTED && n_samples != EXPECTED_SAMPLES) {
    stop("Se esperaban ", EXPECTED_SAMPLES, " muestras; se encontraron ", n_samples, ".")
  }
  if (STRICT_EXPECTED && n_kos != EXPECTED_KOS) {
    stop("Se esperaban ", EXPECTED_KOS, " KOs; se encontraron ", n_kos, ".")
  }

  n_guadeloupe <- if ("guadeloupe_flag" %in% names(metadata)) {
    metadata[guadeloupe_flag == TRUE, .N]
  } else {
    metadata[grepl("guadeloupe", tolower(locality)), .N]
  }
  log_msg("Guadeloupe samples retained: ", n_guadeloupe)

  realm_counts <- metadata[, .(n_metagenomes = .N), by = realm]
  realm_counts[, realm := factor(realm, levels = REALM_ORDER)]
  setorder(realm_counts, realm)
  realm_counts[, realm := as.character(realm)]
  write_tsv(realm_counts, file.path(TAB_DIR, "154_realm_sample_counts.tsv"))

  # CLR por muestra y PCA centrada por KO.
  counts_samples <- t(mat)
  clr_matrix <- log(counts_samples + PSEUDOCOUNT)
  clr_matrix <- sweep(clr_matrix, 1L, rowMeans(clr_matrix), "-")

  if (any(!is.finite(clr_matrix))) stop("CLR produjo valores no finitos.")

  log_msg("Running PCA on CLR-transformed KO profiles")
  pca <- prcomp(clr_matrix, center = TRUE, scale. = FALSE)
  variance_fraction <- pca$sdev^2 / sum(pca$sdev^2)

  scores <- data.table(
    run_accession = rownames(pca$x),
    PC1 = pca$x[, 1L],
    PC2 = pca$x[, 2L]
  )
  scores <- merge(
    scores,
    metadata[, .(run_accession, realm, locality)],
    by = "run_accession",
    all.x = TRUE,
    sort = FALSE
  )
  scores <- scores[match(metadata$run_accession, run_accession)]
  if (!identical(scores$run_accession, metadata$run_accession)) {
    stop("Fallo de alineacion final entre scores PCA y metadata.")
  }

  loadings <- data.table(
    KO = rownames(pca$rotation),
    PC1_loading = pca$rotation[, 1L],
    PC2_loading = pca$rotation[, 2L]
  )

  orientation <- data.table(
    axis = c("PC1", "PC2"),
    multiplier = c(1, 1),
    rule = c(
      "Temperate Northern Atlantic centroid > Central Indo-Pacific centroid",
      "Central Indo-Pacific centroid > Tropical Atlantic centroid"
    )
  )

  realm_mean <- function(axis_name, realm_name) {
    scores[realm == realm_name, mean(get(axis_name), na.rm = TRUE)]
  }

  if (ORIENT_AXES) {
    if (realm_mean("PC1", "Temperate Northern Atlantic") <
        realm_mean("PC1", "Central Indo-Pacific")) {
      scores[, PC1 := -PC1]
      loadings[, PC1_loading := -PC1_loading]
      orientation[axis == "PC1", multiplier := -1]
    }

    if (realm_mean("PC2", "Central Indo-Pacific") <
        realm_mean("PC2", "Tropical Atlantic")) {
      scores[, PC2 := -PC2]
      loadings[, PC2_loading := -PC2_loading]
      orientation[axis == "PC2", multiplier := -1]
    }
  }

  scores[, realm := factor(realm, levels = REALM_ORDER)]
  loadings[, abs_PC1_loading := abs(PC1_loading)]
  loadings[, rank_abs_PC1 := frank(-abs_PC1_loading, ties.method = "first")]

  loadings <- merge(
    loadings,
    annotation[, .(KO, functional_core_class)],
    by = "KO",
    all.x = TRUE,
    sort = FALSE
  )
  if (anyNA(loadings$functional_core_class)) {
    stop("Hay loadings sin clasificacion funcional reconciliada.")
  }

  loadings[, plot_class := fifelse(
    functional_core_class %in% c("multi_realm_non_core", "realm_specific"),
    "functional_periphery",
    functional_core_class
  )]

  class_counts_current <- loadings[, .(n_KOs = .N), by = plot_class]
  expected_summary <- as.data.table(obj153$panel_summary)[, .(plot_class, n_KOs)]
  check_counts <- merge(
    expected_summary,
    class_counts_current,
    by = "plot_class",
    all = TRUE,
    suffixes = c("_run153", "_run154")
  )
  check_counts[, matches :=
    !is.na(n_KOs_run153) &
    !is.na(n_KOs_run154) &
    n_KOs_run153 == n_KOs_run154]
  write_tsv(check_counts, file.path(TAB_DIR, "154_classification_count_audit.tsv"))
  if (!all(check_counts$matches)) {
    stop("Los conteos funcionales no coinciden con el run 153.")
  }

  top100 <- head(loadings[order(rank_abs_PC1)], TOP_N)
  top_contribution <- top100[, .(
    n_KOs = .N,
    weighted_loading = sum(abs_PC1_loading)
  ), by = plot_class]
  top_contribution[, contribution_fraction :=
    weighted_loading / sum(weighted_loading)]

  all_classes <- data.table(plot_class = CLASS_ORDER)
  top_contribution <- merge(
    all_classes,
    top_contribution,
    by = "plot_class",
    all.x = TRUE,
    sort = FALSE
  )
  top_contribution[is.na(n_KOs), n_KOs := 0L]
  top_contribution[is.na(weighted_loading), weighted_loading := 0]
  top_contribution[is.na(contribution_fraction), contribution_fraction := 0]
  top_contribution[, percent_label := sprintf("%.1f%%", 100 * contribution_fraction)]
  top_contribution[, plot_class := factor(plot_class, levels = CLASS_ORDER)]
  setorder(top_contribution, plot_class)
  top_contribution[, plot_class := as.character(plot_class)]

  variance_table <- data.table(
    PC = paste0("PC", seq_along(variance_fraction)),
    variance_fraction = variance_fraction,
    variance_percent = 100 * variance_fraction
  )

  write_tsv(scores, file.path(TAB_DIR, "154_PCA_sample_scores.tsv"))
  write_tsv(variance_table, file.path(TAB_DIR, "154_PCA_variance_explained.tsv"))
  write_tsv(loadings, file.path(TAB_DIR, "154_PC1_loadings_with_classes.tsv"))
  write_tsv(top100, file.path(TAB_DIR, "154_top100_PC1_loadings.tsv"))
  write_tsv(top_contribution, file.path(TAB_DIR, "154_top100_class_contribution.tsv"))
  write_tsv(orientation, file.path(TAB_DIR, "154_axis_orientation.tsv"))

  pc1_percent <- 100 * variance_fraction[1L]
  pc2_percent <- 100 * variance_fraction[2L]
  log_msg("PC1 variance: ", sprintf("%.3f%%", pc1_percent))
  log_msg("PC2 variance: ", sprintf("%.3f%%", pc2_percent))

  if (abs(pc1_percent - 56.5) > 1.0 || abs(pc2_percent - 10.6) > 1.0) {
    log_msg(
      "WARNING: PCA variance differs from the previous 56.5%/10.6% by >1 percentage point."
    )
  }

  realm_breaks <- REALM_ORDER[REALM_ORDER %in% as.character(scores$realm)]

  p_main <- ggplot(scores, aes(PC1, PC2, color = realm)) +
    geom_hline(yintercept = 0, linewidth = 0.25, linetype = 2, color = "#D0D0D0") +
    geom_vline(xintercept = 0, linewidth = 0.25, linetype = 2, color = "#D0D0D0") +
    geom_point(size = 1.55, alpha = 0.78) +
    scale_color_manual(
      values = REALM_PALETTE,
      breaks = realm_breaks,
      labels = unname(REALM_LABELS[realm_breaks]),
      name = NULL,
      drop = FALSE
    ) +
    labs(
      x = sprintf("PC1 (%.1f%%)", pc1_percent),
      y = sprintf("PC2 (%.1f%%)", pc2_percent)
    ) +
    coord_cartesian(clip = "off") +
    theme_panel() +
    theme(
      legend.position = "bottom",
      legend.justification = "center",
      legend.direction = "horizontal",
      legend.margin = margin(1, 0, 0, 0),
      plot.margin = margin(3, 3, 1, 3)
    ) +
    guides(
      color = guide_legend(
        ncol = 2,
        byrow = TRUE,
        override.aes = list(size = 2.2, alpha = 1)
      )
    )

  contribution_plot <- copy(top_contribution[contribution_fraction > 0])
  contribution_plot[, plot_class := factor(
    plot_class,
    levels = c(
      "functional_periphery",
      "absolute_global_core",
      "near_global",
      "realm_core_only"
    )
  )]

  p_inset <- ggplot(
    contribution_plot,
    aes(contribution_fraction, plot_class, fill = plot_class)
  ) +
    geom_col(width = 0.58, color = "white", linewidth = 0.3) +
    geom_text(
      aes(x = contribution_fraction + 0.025, label = percent_label),
      hjust = 0,
      size = TEXT_SIZE_MM,
      color = "black"
    ) +
    scale_x_continuous(
      # Espacio adicional a la derecha para que los porcentajes no se recorten.
      limits = c(0, 1.30),
      breaks = c(0, 0.5, 1),
      labels = percent_format(accuracy = 1),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_discrete(labels = function(x) unname(CLASS_LABELS[x])) +
    scale_fill_manual(values = CLASS_PALETTE, guide = "none") +
    labs(x = "Top-100 loading\ncontribution", y = NULL) +
    theme_classic(base_size = BASE_FONT_SIZE, base_family = "sans") +
    theme(
      text = element_text(color = "black"),
      axis.text = element_text(size = BASE_FONT_SIZE, color = "black"),
      axis.title.x = element_text(
        size = BASE_FONT_SIZE,
        lineheight = 0.90,
        margin = margin(t = 2)
      ),
      axis.line = element_line(linewidth = 0.35, color = "black"),
      axis.ticks = element_line(linewidth = 0.3, color = "black"),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(2, 5, 2, 2)
    )

  combined <- p_main + inset_element(
    p_inset,
    # El inset se ensancha hacia la izquierda para conservar texto de 8 pt.
    left = 0.40,
    bottom = 0.60,
    right = 0.94,
    top = 0.995,
    align_to = "panel",
    on_top = TRUE
  )

  figure_base <- file.path(FIG_DIR, "Figure_4b_functional_PCA_PC1_loadings")
  ggsave(
    paste0(figure_base, ".png"),
    combined,
    width = PANEL_WIDTH_MM,
    height = PANEL_HEIGHT_MM,
    units = "mm",
    dpi = DPI,
    bg = "white"
  )
  ggsave(
    paste0(figure_base, ".pdf"),
    combined,
    width = PANEL_WIDTH_MM,
    height = PANEL_HEIGHT_MM,
    units = "mm",
    device = grDevices::cairo_pdf,
    bg = "white"
  )
  ggsave(
    paste0(figure_base, ".tiff"),
    combined,
    width = PANEL_WIDTH_MM,
    height = PANEL_HEIGHT_MM,
    units = "mm",
    dpi = DPI,
    compression = "lzw",
    bg = "white"
  )

  parameters <- data.table(
    parameter = c(
      "run153",
      "ko_matrix",
      "sample_key",
      "n_samples",
      "n_KOs",
      "pseudocount",
      "PCA_center",
      "PCA_scale",
      "top_n_PC1_loadings",
      "geographic_balancing_applied",
      "guadeloupe_subsampling_applied",
      "guadeloupe_samples_retained",
      "panel_width_mm",
      "panel_height_mm",
      "base_font_size",
      "dpi"
    ),
    value = c(
      RUN153,
      KO_MATRIX,
      "run_accession",
      n_samples,
      n_kos,
      PSEUDOCOUNT,
      "TRUE",
      "FALSE",
      TOP_N,
      "FALSE",
      "FALSE",
      n_guadeloupe,
      PANEL_WIDTH_MM,
      PANEL_HEIGHT_MM,
      BASE_FONT_SIZE,
      DPI
    )
  )
  write_tsv(parameters, file.path(TAB_DIR, "154_parameters.tsv"))

  run_summary <- data.table(
    status = "PASS",
    n_samples = n_samples,
    n_KOs = n_kos,
    n_realms = uniqueN(metadata$realm),
    n_guadeloupe = n_guadeloupe,
    PC1_variance_percent = pc1_percent,
    PC2_variance_percent = pc2_percent,
    top_n = TOP_N,
    figure_png = paste0(figure_base, ".png")
  )
  write_tsv(run_summary, file.path(TAB_DIR, "154_run_summary.tsv"))

  saveRDS(
    list(
      PCA_scores = scores,
      PCA_variance = variance_table,
      PC1_loadings = loadings,
      top100_loadings = top100,
      top100_class_contribution = top_contribution,
      metadata_aligned = metadata,
      axis_orientation = orientation,
      parameters = parameters
    ),
    file.path(RDS_DIR, "154_Figure_4b_analysis_objects.rds")
  )

  capture.output(sessionInfo(), file = file.path(LOG_DIR, "154_sessionInfo.txt"))
  write_status(
    "PASS",
    c(
      paste0("n_samples=", n_samples),
      paste0("n_KOs=", n_kos),
      paste0("PC1_percent=", sprintf("%.4f", pc1_percent)),
      paste0("PC2_percent=", sprintf("%.4f", pc2_percent)),
      paste0("figure=", paste0(figure_base, ".png"))
    )
  )

  log_msg("Figure: ", paste0(figure_base, ".png"))
  log_msg("=== SCRIPT COMPLETED SUCCESSFULLY ===")
}

tryCatch(
  main(),
  error = function(e) {
    write_status("FAILED", paste0("error=", conditionMessage(e)))
    log_msg("ERROR: ", conditionMessage(e))
    stop(e)
  }
)
