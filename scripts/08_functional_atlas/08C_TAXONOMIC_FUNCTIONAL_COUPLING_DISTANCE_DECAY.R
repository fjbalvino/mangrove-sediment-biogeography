#!/usr/bin/env Rscript

# ============================================================================
# 155_FIG4CD_TAXONOMIC_FUNCTIONAL_COUPLING_DECAY_RECONCILED.R
#
# Objetivo
#   Generar los paneles c y d de la Figura 4 usando una cohorte y definiciones
#   reconciliadas:
#
#   Panel c
#     Acoplamiento entre las disimilitudes taxonomica y funcional mediante
#     distancias de Aitchison calculadas sobre perfiles CLR (pseudoconteo = 1).
#     La asociacion se evalua con Mantel de Pearson y 9,999 permutaciones de
#     etiquetas de metagenomas.
#
#   Panel d
#     Comparacion del decaimiento geografico de la similitud taxonomica y
#     funcional. Las similitudes son 1 - Bray-Curtis calculado sobre perfiles
#     de abundancia relativa. Las pendientes se estiman contra log10(km), y su
#     diferencia se prueba mediante 9,999 permutaciones pareadas de etiquetas.
#
# Contratos metodologicos
#   - Llave primaria: run_accession.
#   - Cohorte funcional del run 153: 194 metagenomas de calidad de ensamblaje.
#   - No se aplica submuestreo de Guadalupe ni balanceo geografico.
#   - La misma cohorte y el mismo orden se usan para taxonomia y funcion.
#   - Se comparan perfiles completos; las clases core/periphery no intervienen.
#   - Panel c: Aitchison sobre CLR con pseudoconteo 1.
#   - Panel d: Bray-Curtis sobre abundancias relativas, para obtener similitud
#     acotada entre 0 y 1 y comparable entre taxonomia y funcion.
#   - Pares separados por <=0.1 km se excluyen del panel d.
#   - La prueba de pendientes es una permutacion de etiquetas de muestras, no
#     una prueba "site-block".
#   - Sensibilidad: una muestra por sitio definido por coordenadas.
#   - Sin titulo ni letra de panel; salidas individuales de 87 x 100 mm.
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(phyloseq)
  library(vegan)
  library(geosphere)
  library(ggplot2)
  library(scales)
  library(patchwork)
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
  x <- readLines(path, warn = FALSE)
  if (length(x) == 0L) return(NA_character_)
  trimws(x[1L])
}

SCRIPT_ID <- "155_FIG4CD_TAXONOMIC_FUNCTIONAL_COUPLING_DECAY_RECONCILED"

RUN153 <- get_arg(
  "--run153",
  read_latest("153_FIG4A_FUNCTIONAL_KO_CLASSES_RECONCILED")
)

KO_MATRIX <- get_arg(
  "--ko_matrix",
  "/data/Magenta/results/035.metagenome_assembly_annotation/final_results/MATRIX_KEGG_ko.tsv"
)

TAXONOMIC_PHYLOSEQ <- get_arg(
  "--taxonomic_phyloseq",
  "/data/Magenta/phyloseq_MAGENTA.rds"
)

METADATA_FILE <- get_arg(
  "--metadata",
  "/data/Magenta/data/metadata/Metadata_2026_copia.csv"
)

PSEUDOCOUNT <- as.numeric(get_arg("--pseudocount", "1"))
MIN_GEO_KM <- as.numeric(get_arg("--min_geo_km", "0.1"))
N_PERM <- as.integer(get_arg("--permutations", "9999"))
SENSITIVITY_REPS <- as.integer(get_arg("--site_sensitivity_reps", "1000"))
COORD_DIGITS <- as.integer(get_arg("--coordinate_digits", "6"))
EXPECTED_SAMPLES <- as.integer(get_arg("--expected_samples", "194"))
EXPECTED_KOS <- as.integer(get_arg("--expected_kos", "17508"))
STRICT_EXPECTED <- as.integer(get_arg("--strict_expected", "1")) == 1L
SEED <- as.integer(get_arg("--seed", "123"))
PANEL_WIDTH_MM <- as.numeric(get_arg("--panel_width_mm", "87"))
PANEL_HEIGHT_MM <- as.numeric(get_arg("--panel_height_mm", "100"))
BASE_FONT_SIZE <- as.numeric(get_arg("--base_font_size", "8"))
DPI <- as.integer(get_arg("--dpi", "600"))
TEXT_SIZE_MM <- BASE_FONT_SIZE / 2.845276

if (!is.finite(PSEUDOCOUNT) || PSEUDOCOUNT <= 0) {
  stop("pseudocount debe ser > 0.")
}
if (!is.finite(MIN_GEO_KM) || MIN_GEO_KM < 0) {
  stop("min_geo_km debe ser >= 0.")
}
if (!is.finite(N_PERM) || N_PERM < 99L) {
  stop("permutations debe ser >= 99.")
}
if (!is.finite(SENSITIVITY_REPS) || SENSITIVITY_REPS < 0L) {
  stop("site_sensitivity_reps debe ser >= 0.")
}
if (!is.finite(BASE_FONT_SIZE) || BASE_FONT_SIZE < 8 || BASE_FONT_SIZE > 11) {
  stop("base_font_size debe estar entre 8 y 11 pt.")
}
if (!is.finite(PANEL_WIDTH_MM) || PANEL_WIDTH_MM <= 0 ||
    !is.finite(PANEL_HEIGHT_MM) || PANEL_HEIGHT_MM <= 0) {
  stop("Las dimensiones del panel deben ser positivas.")
}

set.seed(SEED)

# -----------------------------
# Directorios, log y estado
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
writeLines(OUTDIR, LATEST_FILE)

LOG_FILE <- file.path(LOG_DIR, "155_run.log")
log_con <- file(LOG_FILE, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")

on.exit({
  while (sink.number(type = "message") > 0L) sink(type = "message")
  while (sink.number(type = "output") > 0L) sink(type = "output")
  close(log_con)
}, add = TRUE)

log_msg <- function(...) {
  cat(
    sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0(..., collapse = ""),
    "
",
    sep = ""
  )
}

write_tsv <- function(x, path) {
  data.table::fwrite(
    as.data.table(x),
    file = path,
    sep = "	",
    quote = FALSE,
    na = "NA"
  )
}

write_status <- function(status, extra = character()) {
  writeLines(
    c(
      status,
      paste0("timestamp=", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      paste0("output_directory=", OUTDIR),
      extra
    ),
    file.path(OUTDIR, "STATUS.txt")
  )
}

# -----------------------------
# Utilidades de IDs, columnas y tablas
# -----------------------------
is_missing_like <- function(x) {
  y <- trimws(as.character(x))
  is.na(x) | y == "" | tolower(y) %in% c("na", "n/a", "nan", "null", "none")
}

clean_chr <- function(x) {
  y <- trimws(as.character(x))
  y[is_missing_like(y)] <- NA_character_
  y
}

norm_name <- function(x) {
  y <- tolower(trimws(as.character(x)))
  y <- gsub("[^a-z0-9]+", "_", y)
  gsub("^_+|_+$", "", y)
}

norm_id <- function(x) {
  y <- clean_chr(x)
  y <- toupper(y)
  gsub("[^A-Z0-9]+", "", y)
}

first_nonmissing <- function(x) {
  y <- clean_chr(x)
  y <- y[!is.na(y)]
  if (length(y) == 0L) return(NA_character_)
  y[1L]
}

numeric_safe <- function(x) {
  suppressWarnings(as.numeric(gsub(",", ".", clean_chr(x), fixed = TRUE)))
}

resolve_col <- function(dt, candidates, label, required = TRUE) {
  nm <- names(dt)
  nn <- norm_name(nm)
  cc <- norm_name(candidates)
  idx <- match(cc, nn)
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0L) {
    if (required) {
      stop("No se encontro la columna ", label, ". Candidatas: ",
           paste(candidates, collapse = ", "))
    }
    return(NA_character_)
  }
  nm[idx[1L]]
}

read_table_auto <- function(path) {
  if (!file.exists(path)) stop("No existe: ", path)
  ext <- tolower(tools::file_ext(path))
  sep <- if (ext %in% c("tsv", "txt")) "	" else ","
  data.table::fread(path, sep = sep, data.table = TRUE, check.names = FALSE)
}

match_ids_strict <- function(query_ids, reference_ids, context) {
  if (anyNA(query_ids) || anyNA(reference_ids)) {
    stop(context, ": hay IDs faltantes.")
  }
  if (anyDuplicated(reference_ids)) {
    stop(context, ": la referencia contiene run_accession duplicados.")
  }
  ref_norm <- norm_id(reference_ids)
  if (anyDuplicated(ref_norm)) {
    stop(context, ": run_accession colisionan despues de normalizar IDs.")
  }
  idx <- match(norm_id(query_ids), ref_norm)
  if (anyNA(idx)) {
    missing <- query_ids[is.na(idx)]
    stop(
      context, ": faltan ", length(missing), " IDs. Primeros: ",
      paste(head(missing, 20L), collapse = ", ")
    )
  }
  idx
}

# -----------------------------
# Lectura robusta de la matriz KO
# -----------------------------
load_ko_matrix <- function(path, reference_ids) {
  dt <- read_table_auto(path)
  if (nrow(dt) == 0L || ncol(dt) < 2L) stop("La matriz KO esta vacia.")

  column_hits <- which(norm_id(names(dt)) %in% norm_id(reference_ids))
  row_candidate <- resolve_col(
    dt,
    c("run_accession", "sample_id", "sample", "id"),
    "ID de muestra en matriz KO",
    required = FALSE
  )

  if (length(column_hits) > 0L) {
    ko_col <- setdiff(seq_along(dt), column_hits)[1L]
    if (is.na(ko_col)) stop("No se encontro la columna de KO.")
    ko_ids <- clean_chr(dt[[ko_col]])
    sample_names_raw <- names(dt)[column_hits]
    sample_idx <- match_ids_strict(
      sample_names_raw,
      reference_ids,
      "Columnas de la matriz KO"
    )
    canonical_samples <- reference_ids[sample_idx]
    if (anyDuplicated(canonical_samples)) {
      stop("Varias columnas KO se asociaron al mismo run_accession.")
    }
    mat <- as.matrix(dt[, ..column_hits])
    storage.mode(mat) <- "numeric"
    rownames(mat) <- ko_ids
    colnames(mat) <- canonical_samples
    orientation <- "features_by_samples"
  } else if (!is.na(row_candidate)) {
    sample_ids_raw <- clean_chr(dt[[row_candidate]])
    row_hits <- which(norm_id(sample_ids_raw) %in% norm_id(reference_ids))
    if (length(row_hits) == 0L) {
      stop("No se asociaron filas KO por run_accession.")
    }
    feature_cols <- setdiff(seq_along(dt), match(row_candidate, names(dt)))
    sample_idx <- match_ids_strict(
      sample_ids_raw[row_hits],
      reference_ids,
      "Filas de la matriz KO"
    )
    canonical_samples <- reference_ids[sample_idx]
    if (anyDuplicated(canonical_samples)) {
      stop("Varias filas KO se asociaron al mismo run_accession.")
    }
    tmp <- as.matrix(dt[row_hits, ..feature_cols])
    storage.mode(tmp) <- "numeric"
    mat <- t(tmp)
    rownames(mat) <- names(dt)[feature_cols]
    colnames(mat) <- canonical_samples
    orientation <- "samples_by_features"
  } else {
    stop("No se pudo resolver la orientacion de MATRIX_KEGG_ko.tsv.")
  }

  if (anyNA(rownames(mat)) || any(!nzchar(rownames(mat)))) {
    stop("La matriz KO contiene IDs de KO faltantes.")
  }
  if (anyDuplicated(rownames(mat))) {
    log_msg("Collapsing duplicated KO identifiers by sum")
    mat <- t(rowsum(t(mat), group = rownames(mat), reorder = FALSE))
  }
  if (any(!is.finite(mat)) || any(mat < 0)) {
    stop("La matriz KO contiene valores negativos o no finitos.")
  }

  matching <- data.table(
    matrix_sample = if (orientation == "features_by_samples") {
      sample_names_raw
    } else {
      sample_ids_raw[row_hits]
    },
    run_accession = canonical_samples,
    orientation = orientation
  )

  list(matrix = mat, matching = matching)
}

# -----------------------------
# Transformaciones y estadistica
# -----------------------------
clr_transform <- function(sample_by_feature, pseudocount) {
  x <- as.matrix(sample_by_feature)
  if (any(!is.finite(x)) || any(x < 0)) {
    stop("La matriz contiene valores negativos o no finitos antes de CLR.")
  }
  lx <- log(x + pseudocount)
  sweep(lx, 1L, rowMeans(lx), "-")
}

relative_abundance <- function(sample_by_feature) {
  x <- as.matrix(sample_by_feature)
  totals <- rowSums(x)
  if (any(!is.finite(totals)) || any(totals <= 0)) {
    stop("Hay muestras con abundancia total no positiva.")
  }
  sweep(x, 1L, totals, "/")
}

slope_fast <- function(y, x) {
  keep <- is.finite(y) & is.finite(x)
  y <- y[keep]
  x <- x[keep]
  denom <- sum((x - mean(x))^2)
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  sum((x - mean(x)) * (y - mean(y))) / denom
}

format_p <- function(p, digits = 3L) {
  threshold <- 10^(-digits)
  if (!is.finite(p)) return("P = NA")
  if (p < threshold) return(paste0("P < ", format(threshold, scientific = FALSE)))
  paste0("P = ", formatC(p, format = "f", digits = digits))
}

permutation_slope_difference <- function(
    sim_tax_mat,
    sim_func_mat,
    geo_mat,
    min_geo_km,
    nperm,
    seed) {
  n <- nrow(geo_mat)
  if (!identical(dim(sim_tax_mat), dim(sim_func_mat)) ||
      !identical(dim(sim_tax_mat), dim(geo_mat))) {
    stop("Matrices incompatibles en la prueba de pendientes.")
  }

  upper <- upper.tri(geo_mat)
  eligible <- upper & is.finite(geo_mat) & geo_mat > min_geo_km
  log_geo <- log10(geo_mat[eligible])
  sim_tax <- sim_tax_mat[eligible]
  sim_func <- sim_func_mat[eligible]

  slope_tax <- slope_fast(sim_tax, log_geo)
  slope_func <- slope_fast(sim_func, log_geo)
  observed <- slope_tax - slope_func

  delta_mat <- sim_tax_mat - sim_func_mat
  permuted <- numeric(nperm)
  set.seed(seed)

  for (i in seq_len(nperm)) {
    perm <- sample.int(n)
    delta_perm <- delta_mat[perm, perm, drop = FALSE]
    permuted[i] <- slope_fast(delta_perm[eligible], log_geo)
    if (i %% 1000L == 0L) log_msg("Slope permutations: ", i, "/", nperm)
  }

  p_two_sided <- (sum(abs(permuted) >= abs(observed)) + 1) / (nperm + 1)
  p_tax_more_negative <- (sum(permuted <= observed) + 1) / (nperm + 1)

  list(
    slope_taxonomic = slope_tax,
    slope_functional = slope_func,
    difference_tax_minus_function = observed,
    p_two_sided = p_two_sided,
    p_taxonomic_more_negative = p_tax_more_negative,
    n_pairs = sum(eligible),
    eligible = eligible,
    log_geo = log_geo,
    permuted = permuted
  )
}

# -----------------------------
# Tema grafico
# -----------------------------
theme_panel <- function() {
  theme_classic(base_size = BASE_FONT_SIZE, base_family = "sans") +
    theme(
      text = element_text(color = "black"),
      axis.title = element_text(size = BASE_FONT_SIZE, color = "black"),
      axis.text = element_text(size = BASE_FONT_SIZE, color = "black"),
      axis.line = element_line(linewidth = 0.45, color = "black"),
      axis.ticks = element_line(linewidth = 0.35, color = "black"),
      legend.text = element_text(size = BASE_FONT_SIZE, color = "black"),
      legend.key.height = grid::unit(3.5, "mm"),
      legend.key.width = grid::unit(7, "mm"),
      plot.margin = margin(3, 3, 3, 3)
    )
}

PLOT_COLORS <- c(
  Taxonomic = "#0072B2",
  Functional = "#009E73"
)

# -----------------------------
# Analisis principal
# -----------------------------
write_status("RUNNING")

result <- tryCatch({
  log_msg("Starting ", SCRIPT_ID)
  log_msg("Run153: ", RUN153)
  log_msg("KO matrix: ", KO_MATRIX)
  log_msg("Taxonomic phyloseq: ", TAXONOMIC_PHYLOSEQ)
  log_msg("Metadata: ", METADATA_FILE)
  log_msg("Output: ", OUTDIR)

  if (is.na(RUN153) || !nzchar(RUN153) || !dir.exists(RUN153)) {
    stop("run153 no es valido: ", RUN153)
  }
  if (!file.exists(KO_MATRIX)) stop("No existe KO matrix: ", KO_MATRIX)
  if (!file.exists(TAXONOMIC_PHYLOSEQ)) {
    stop("No existe phyloseq taxonomico: ", TAXONOMIC_PHYLOSEQ)
  }
  if (!file.exists(METADATA_FILE)) stop("No existe metadata: ", METADATA_FILE)

  object153_file <- file.path(
    RUN153,
    "rds",
    "153_Figure_4a_analysis_objects.rds"
  )
  if (!file.exists(object153_file)) stop("No existe: ", object153_file)

  obj153 <- readRDS(object153_file)
  if (!"metadata_aligned" %in% names(obj153)) {
    stop("El objeto 153 no contiene metadata_aligned.")
  }

  cohort <- as.data.table(obj153$metadata_aligned)
  required_cohort <- c("run_accession", "realm", "locality")
  if (length(setdiff(required_cohort, names(cohort))) > 0L) {
    stop("metadata_aligned del run153 no contiene las columnas requeridas.")
  }
  cohort[, run_accession := clean_chr(run_accession)]
  cohort[, realm := clean_chr(realm)]
  cohort[, locality := clean_chr(locality)]
  if (anyNA(cohort$run_accession) || anyDuplicated(cohort$run_accession)) {
    stop("La cohorte 153 contiene run_accession faltantes o duplicados.")
  }
  if (anyNA(cohort$realm) || anyNA(cohort$locality)) {
    stop("La cohorte 153 contiene realm/locality faltantes.")
  }
  cohort_ids <- cohort$run_accession

  # 1. Matriz funcional, exactamente en el orden del run153.
  loaded_ko <- load_ko_matrix(KO_MATRIX, cohort_ids)
  write_tsv(loaded_ko$matching, file.path(TAB_DIR, "155_KO_sample_matching.tsv"))
  ko_mat <- loaded_ko$matrix
  ko_idx <- match_ids_strict(cohort_ids, colnames(ko_mat), "Cohorte vs KO matrix")
  ko_mat <- ko_mat[, ko_idx, drop = FALSE]
  colnames(ko_mat) <- cohort_ids
  ko_mat <- ko_mat[rowSums(ko_mat) > 0, , drop = FALSE]

  if (STRICT_EXPECTED && length(cohort_ids) != EXPECTED_SAMPLES) {
    stop("Se esperaban ", EXPECTED_SAMPLES, " muestras; se encontraron ",
         length(cohort_ids), ".")
  }
  if (STRICT_EXPECTED && nrow(ko_mat) != EXPECTED_KOS) {
    stop("Se esperaban ", EXPECTED_KOS, " KOs; se encontraron ",
         nrow(ko_mat), ".")
  }

  # 2. Matriz taxonomica alineada por run_accession.
  ps <- readRDS(TAXONOMIC_PHYLOSEQ)
  ps_sample_names <- sample_names(ps)
  sd <- as.data.table(as(sample_data(ps), "data.frame"), keep.rownames = "ps_sample_name")
  sd <- sd[match(ps_sample_names, ps_sample_name)]
  if (!identical(sd$ps_sample_name, ps_sample_names)) {
    stop("No fue posible alinear sample_data con sample_names del phyloseq.")
  }
  run_col_ps <- resolve_col(
    sd,
    c("run_accession", "Run_accession", "RUN_ACCESSION"),
    "run_accession en phyloseq",
    required = FALSE
  )
  ps_run_ids <- if (is.na(run_col_ps)) {
    clean_chr(ps_sample_names)
  } else {
    clean_chr(sd[[run_col_ps]])
  }
  if (length(ps_run_ids) != length(ps_sample_names)) {
    stop("Longitudes incompatibles entre sample_data y sample_names.")
  }
  if (anyNA(ps_run_ids) || anyDuplicated(ps_run_ids)) {
    stop("run_accession faltantes o duplicados en el phyloseq.")
  }

  otu <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) otu <- t(otu)
  otu <- otu[ps_sample_names, , drop = FALSE]
  rownames(otu) <- ps_run_ids

  tax_idx <- match_ids_strict(cohort_ids, rownames(otu), "Cohorte vs phyloseq")
  tax_mat <- otu[tax_idx, , drop = FALSE]
  rownames(tax_mat) <- cohort_ids
  tax_mat <- tax_mat[, colSums(tax_mat) > 0, drop = FALSE]
  storage.mode(tax_mat) <- "numeric"

  if (!identical(rownames(tax_mat), colnames(ko_mat))) {
    stop("Fallo de alineacion final taxonomia-funcion por run_accession.")
  }
  if (any(!is.finite(tax_mat)) || any(tax_mat < 0)) {
    stop("La matriz taxonomica contiene valores negativos o no finitos.")
  }

  # 3. Metadata geografica deduplicada por run_accession.
  metadata_raw <- read_table_auto(METADATA_FILE)
  run_col <- resolve_col(metadata_raw, c("run_accession"), "run_accession")
  lon_col <- resolve_col(
    metadata_raw,
    c(
      "longitude", "lon", "long", "Longitude",
      "decimalLongitude", "decimal_longitude"
    ),
    "longitude"
  )
  lat_col <- resolve_col(
    metadata_raw,
    c("latitude", "lat", "Latitude", "decimalLatitude", "decimal_latitude"),
    "latitude"
  )

  geo_raw <- metadata_raw[, .(
    run_accession = clean_chr(get(run_col)),
    longitude = numeric_safe(get(lon_col)),
    latitude = numeric_safe(get(lat_col))
  )]
  geo_raw <- geo_raw[!is.na(run_accession)]

  geo_audit <- geo_raw[, .(
    n_rows = .N,
    n_complete_coordinates = sum(is.finite(longitude) & is.finite(latitude)),
    n_unique_coordinate_pairs = uniqueN(
      paste(longitude[is.finite(longitude) & is.finite(latitude)],
            latitude[is.finite(longitude) & is.finite(latitude)], sep = "|")
    ),
    longitude_min = if (any(is.finite(longitude))) min(longitude, na.rm = TRUE) else NA_real_,
    longitude_max = if (any(is.finite(longitude))) max(longitude, na.rm = TRUE) else NA_real_,
    latitude_min = if (any(is.finite(latitude))) min(latitude, na.rm = TRUE) else NA_real_,
    latitude_max = if (any(is.finite(latitude))) max(latitude, na.rm = TRUE) else NA_real_
  ), by = run_accession]
  write_tsv(geo_audit, file.path(TAB_DIR, "155_metadata_coordinate_audit.tsv"))

  geo_unique <- geo_raw[, .(
    longitude = if (any(is.finite(longitude))) median(longitude, na.rm = TRUE) else NA_real_,
    latitude = if (any(is.finite(latitude))) median(latitude, na.rm = TRUE) else NA_real_
  ), by = run_accession]

  geo_idx <- match(norm_id(cohort_ids), norm_id(geo_unique$run_accession))
  cohort[, longitude := geo_unique$longitude[geo_idx]]
  cohort[, latitude := geo_unique$latitude[geo_idx]]
  cohort[, has_coordinates :=
    is.finite(longitude) & is.finite(latitude) &
    longitude >= -180 & longitude <= 180 &
    latitude >= -90 & latitude <= 90]
  cohort[, site_id := fifelse(
    has_coordinates == TRUE,
    paste0(
      format(round(longitude, COORD_DIGITS), nsmall = COORD_DIGITS), "|",
      format(round(latitude, COORD_DIGITS), nsmall = COORD_DIGITS)
    ),
    NA_character_
  )]
  cohort[, geographic_balancing_applied := FALSE]
  cohort[, guadeloupe_subsampling_applied := FALSE]
  write_tsv(cohort, file.path(TAB_DIR, "155_cohort_audit.tsv"))

  n_samples <- nrow(cohort)
  n_taxa <- ncol(tax_mat)
  n_kos <- nrow(ko_mat)
  n_geo_samples <- sum(cohort$has_coordinates)
  n_sites <- uniqueN(cohort$site_id[cohort$has_coordinates == TRUE])

  log_msg("Common samples: ", n_samples)
  log_msg("Taxonomic features retained in common cohort: ", n_taxa)
  log_msg("KOs retained in common cohort: ", n_kos)
  log_msg("Samples with coordinates: ", n_geo_samples)
  log_msg("Coordinate-defined sites: ", n_sites)

  if (n_geo_samples < 3L) stop("No hay suficientes muestras con coordenadas.")

  # 4. Panel c: Aitchison taxonomia vs funcion.
  ko_samples <- t(ko_mat)
  tax_clr <- clr_transform(tax_mat, PSEUDOCOUNT)
  func_clr <- clr_transform(ko_samples, PSEUDOCOUNT)

  dist_tax_ait <- stats::dist(tax_clr, method = "euclidean")
  dist_func_ait <- stats::dist(func_clr, method = "euclidean")

  set.seed(SEED)
  mantel_c <- vegan::mantel(
    dist_tax_ait,
    dist_func_ait,
    method = "pearson",
    permutations = N_PERM,
    na.rm = FALSE
  )

  panel_c_pairs <- data.table(
    taxonomic_aitchison = as.vector(dist_tax_ait),
    functional_aitchison = as.vector(dist_func_ait)
  )
  model_c <- lm(functional_aitchison ~ taxonomic_aitchison, data = panel_c_pairs)
  write_tsv(panel_c_pairs, file.path(TAB_DIR, "155_panel_c_Aitchison_pairs.tsv"))

  # 5. Panel d: Bray-Curtis de abundancia relativa y distancia geografica.
  tax_rel <- relative_abundance(tax_mat)
  func_rel <- relative_abundance(ko_samples)
  dist_tax_bray_all <- as.matrix(vegdist(tax_rel, method = "bray"))
  dist_func_bray_all <- as.matrix(vegdist(func_rel, method = "bray"))
  rownames(dist_tax_bray_all) <- colnames(dist_tax_bray_all) <- cohort_ids
  rownames(dist_func_bray_all) <- colnames(dist_func_bray_all) <- cohort_ids

  geo_ids <- cohort[has_coordinates == TRUE, run_accession]
  geo_order <- match(geo_ids, cohort_ids)
  coords <- as.matrix(cohort[has_coordinates == TRUE, .(longitude, latitude)])
  rownames(coords) <- geo_ids
  geo_mat <- geosphere::distm(coords, fun = geosphere::distHaversine) / 1000
  rownames(geo_mat) <- colnames(geo_mat) <- geo_ids

  dist_tax_bray <- dist_tax_bray_all[geo_order, geo_order, drop = FALSE]
  dist_func_bray <- dist_func_bray_all[geo_order, geo_order, drop = FALSE]
  sim_tax_mat <- 1 - dist_tax_bray
  sim_func_mat <- 1 - dist_func_bray

  slope_test <- permutation_slope_difference(
    sim_tax_mat = sim_tax_mat,
    sim_func_mat = sim_func_mat,
    geo_mat = geo_mat,
    min_geo_km = MIN_GEO_KM,
    nperm = N_PERM,
    seed = SEED + 1L
  )

  eligible <- slope_test$eligible
  panel_d_pairs <- rbind(
    data.table(
      geographic_distance_km = geo_mat[eligible],
      log10_geographic_distance = log10(geo_mat[eligible]),
      similarity = sim_tax_mat[eligible],
      profile = "Taxonomic"
    ),
    data.table(
      geographic_distance_km = geo_mat[eligible],
      log10_geographic_distance = log10(geo_mat[eligible]),
      similarity = sim_func_mat[eligible],
      profile = "Functional"
    )
  )
  panel_d_pairs[, profile := factor(profile, levels = c("Taxonomic", "Functional"))]
  write_tsv(panel_d_pairs, file.path(TAB_DIR, "155_panel_d_distance_decay_pairs.tsv"))
  write_tsv(
    data.table(permuted_slope_difference = slope_test$permuted),
    file.path(TAB_DIR, "155_panel_d_permuted_slope_differences.tsv")
  )

  # 6. Sensibilidad: una muestra por sitio definido por coordenadas.
  sensitivity <- data.table()
  if (SENSITIVITY_REPS > 0L) {
    geo_cohort <- cohort[has_coordinates == TRUE]
    geo_global_idx <- match(geo_cohort$run_accession, cohort_ids)
    tax_ait_mat_all <- as.matrix(dist_tax_ait)
    func_ait_mat_all <- as.matrix(dist_func_ait)
    set.seed(SEED + 2L)

    sensitivity_list <- vector("list", SENSITIVITY_REPS)
    for (b in seq_len(SENSITIVITY_REPS)) {
      selected_rows <- geo_cohort[, sample(.I, 1L), by = site_id]$V1
      selected_global <- geo_global_idx[selected_rows]
      selected_geo <- selected_rows

      upper_c <- upper.tri(tax_ait_mat_all[selected_global, selected_global, drop = FALSE])
      r_c <- suppressWarnings(cor(
        tax_ait_mat_all[selected_global, selected_global, drop = FALSE][upper_c],
        func_ait_mat_all[selected_global, selected_global, drop = FALSE][upper_c],
        method = "pearson"
      ))

      geo_b <- geo_mat[selected_geo, selected_geo, drop = FALSE]
      tax_sim_b <- sim_tax_mat[selected_geo, selected_geo, drop = FALSE]
      func_sim_b <- sim_func_mat[selected_geo, selected_geo, drop = FALSE]
      eligible_b <- upper.tri(geo_b) & is.finite(geo_b) & geo_b > MIN_GEO_KM
      x_b <- log10(geo_b[eligible_b])
      slope_tax_b <- slope_fast(tax_sim_b[eligible_b], x_b)
      slope_func_b <- slope_fast(func_sim_b[eligible_b], x_b)

      sensitivity_list[[b]] <- data.table(
        replicate = b,
        n_sites = length(selected_rows),
        coupling_r = r_c,
        taxonomic_slope = slope_tax_b,
        functional_slope = slope_func_b,
        slope_difference_tax_minus_function = slope_tax_b - slope_func_b
      )
    }
    sensitivity <- rbindlist(sensitivity_list)
    write_tsv(sensitivity, file.path(TAB_DIR, "155_one_sample_per_site_sensitivity.tsv"))

    sensitivity_summary <- rbindlist(lapply(
      c(
        "coupling_r",
        "taxonomic_slope",
        "functional_slope",
        "slope_difference_tax_minus_function"
      ),
      function(metric) {
        x <- sensitivity[[metric]]
        data.table(
          metric = metric,
          n_replicates = sum(is.finite(x)),
          median = median(x, na.rm = TRUE),
          q025 = quantile(x, 0.025, na.rm = TRUE, names = FALSE),
          q975 = quantile(x, 0.975, na.rm = TRUE, names = FALSE),
          minimum = min(x, na.rm = TRUE),
          maximum = max(x, na.rm = TRUE)
        )
      }
    ))
    write_tsv(
      sensitivity_summary,
      file.path(TAB_DIR, "155_one_sample_per_site_sensitivity_summary.tsv")
    )
  }

  # 7. Tablas de estadisticas.
  stats_table <- data.table(
    panel = c("c", "d", "d", "d"),
    analysis = c(
      "Taxonomic-functional Aitchison coupling",
      "Taxonomic geographic distance-decay",
      "Functional geographic distance-decay",
      "Slope difference: taxonomic minus functional"
    ),
    n_samples = c(n_samples, n_geo_samples, n_geo_samples, n_geo_samples),
    n_pairs = c(
      nrow(panel_c_pairs),
      slope_test$n_pairs,
      slope_test$n_pairs,
      slope_test$n_pairs
    ),
    statistic = c(
      unname(mantel_c$statistic),
      slope_test$slope_taxonomic,
      slope_test$slope_functional,
      slope_test$difference_tax_minus_function
    ),
    statistic_name = c("Mantel r", "slope", "slope", "delta slope"),
    p_value = c(
      mantel_c$signif,
      NA_real_,
      NA_real_,
      slope_test$p_two_sided
    ),
    permutations = c(N_PERM, NA_integer_, NA_integer_, N_PERM)
  )
  write_tsv(stats_table, file.path(TAB_DIR, "155_Figure_4cd_statistics.tsv"))

  # 8. Figuras.
  c_label <- paste0(
    "Mantel r = ", sprintf("%.3f", unname(mantel_c$statistic)),
    ", ", format_p(mantel_c$signif)
  )

  c_x <- min(panel_c_pairs$taxonomic_aitchison) +
    0.03 * diff(range(panel_c_pairs$taxonomic_aitchison))
  c_y <- max(panel_c_pairs$functional_aitchison) -
    0.04 * diff(range(panel_c_pairs$functional_aitchison))

  p_c <- ggplot(
    panel_c_pairs,
    aes(taxonomic_aitchison, functional_aitchison)
  ) +
    geom_point(
      size = 0.60,
      alpha = 0.20,
      color = "#4C78A8",
      stroke = 0
    ) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      linewidth = 0.75,
      color = "#1F3A5F"
    ) +
    annotate(
      "text",
      x = c_x,
      y = c_y,
      label = c_label,
      hjust = 0,
      vjust = 1,
      size = TEXT_SIZE_MM,
      color = "black"
    ) +
    labs(
      x = "Taxonomic Aitchison dissimilarity",
      y = "Functional Aitchison dissimilarity"
    ) +
    theme_panel()

  tax_slope_label <- paste0(
    "Taxonomic slope = ", sprintf("%.3f", slope_test$slope_taxonomic)
  )
  func_slope_label <- paste0(
    "Functional slope = ", sprintf("%.3f", slope_test$slope_functional)
  )
  slope_difference_label <- paste0(
    "Slope difference = ",
    sprintf("%.3f", slope_test$difference_tax_minus_function),
    ", ", format_p(slope_test$p_two_sided)
  )

  x_min_d <- min(panel_d_pairs$geographic_distance_km)
  x_max_d <- max(panel_d_pairs$geographic_distance_km)
  x_grid <- exp(seq(log(x_min_d), log(x_max_d), length.out = 250L))
  line_data <- rbind(
    data.table(
      geographic_distance_km = x_grid,
      similarity = coef(lm(
        similarity ~ log10_geographic_distance,
        data = panel_d_pairs[profile == "Taxonomic"]
      ))[1L] + slope_test$slope_taxonomic * log10(x_grid),
      profile = factor("Taxonomic", levels = c("Taxonomic", "Functional"))
    ),
    data.table(
      geographic_distance_km = x_grid,
      similarity = coef(lm(
        similarity ~ log10_geographic_distance,
        data = panel_d_pairs[profile == "Functional"]
      ))[1L] + slope_test$slope_functional * log10(x_grid),
      profile = factor("Functional", levels = c("Taxonomic", "Functional"))
    )
  )

  p_d <- ggplot(
    panel_d_pairs,
    aes(geographic_distance_km, similarity, color = profile)
  ) +
    geom_point(size = 0.55, alpha = 0.18, stroke = 0, show.legend = FALSE) +
    geom_line(
      data = line_data,
      aes(geographic_distance_km, similarity, color = profile),
      linewidth = 0.80,
      inherit.aes = FALSE
    ) +
    annotate(
      "text",
      x = x_min_d * 1.18,
      y = 0.22,
      label = tax_slope_label,
      hjust = 0,
      vjust = 0,
      size = TEXT_SIZE_MM,
      color = PLOT_COLORS[["Taxonomic"]]
    ) +
    annotate(
      "text",
      x = x_min_d * 1.18,
      y = 0.17,
      label = func_slope_label,
      hjust = 0,
      vjust = 0,
      size = TEXT_SIZE_MM,
      color = PLOT_COLORS[["Functional"]]
    ) +
    annotate(
      "text",
      x = x_min_d * 1.18,
      y = 0.12,
      label = slope_difference_label,
      hjust = 0,
      vjust = 0,
      size = TEXT_SIZE_MM,
      color = "black"
    ) +
    scale_x_log10(
      breaks = scales::log_breaks(base = 10),
      labels = scales::label_number(big.mark = ",")
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.25),
      expand = expansion(mult = c(0, 0.02))
    ) +
    scale_color_manual(
      values = PLOT_COLORS,
      breaks = c("Taxonomic", "Functional"),
      name = NULL
    ) +
    labs(
      x = expression("Geographic distance (km, "*log[10]*" scale)"),
      y = "Community similarity (1 - Bray-Curtis)"
    ) +
    theme_panel() +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.margin = margin(1, 0, 0, 0)
    ) +
    guides(
      color = guide_legend(
        override.aes = list(linewidth = 1.0, alpha = 1)
      )
    )

  figure_c_base <- file.path(
    FIG_DIR,
    "Figure_4c_taxonomic_functional_Aitchison_coupling"
  )
  figure_d_base <- file.path(
    FIG_DIR,
    "Figure_4d_taxonomic_functional_distance_decay"
  )

  save_panel <- function(plot, path_base, width_mm, height_mm) {
    ggsave(
      paste0(path_base, ".png"), plot,
      width = width_mm, height = height_mm, units = "mm", dpi = DPI,
      bg = "white"
    )
    ggsave(
      paste0(path_base, ".pdf"), plot,
      width = width_mm, height = height_mm, units = "mm",
      device = grDevices::cairo_pdf, bg = "white"
    )
    ggsave(
      paste0(path_base, ".tiff"), plot,
      width = width_mm, height = height_mm, units = "mm", dpi = DPI,
      compression = "lzw", bg = "white"
    )
  }

  save_panel(p_c, figure_c_base, PANEL_WIDTH_MM, PANEL_HEIGHT_MM)
  save_panel(p_d, figure_d_base, PANEL_WIDTH_MM, PANEL_HEIGHT_MM)

  combined <- p_c + p_d + plot_layout(ncol = 2, guides = "keep")
  combined_base <- file.path(FIG_DIR, "Figure_4cd_combined")
  save_panel(combined, combined_base, 2 * PANEL_WIDTH_MM, PANEL_HEIGHT_MM)

  # 9. Resumen, parametros y objeto reproducible.
  parameters <- data.table(
    parameter = c(
      "run153",
      "KO_matrix",
      "taxonomic_phyloseq",
      "metadata",
      "primary_key",
      "functional_cohort_rule",
      "Guadeloupe_subsampling",
      "panel_c_distance",
      "panel_c_transformation",
      "panel_c_test",
      "panel_d_distance",
      "panel_d_similarity",
      "panel_d_predictor",
      "panel_d_slope_test",
      "minimum_geographic_distance_km",
      "permutations",
      "site_sensitivity_replicates",
      "seed"
    ),
    value = c(
      RUN153,
      KO_MATRIX,
      TAXONOMIC_PHYLOSEQ,
      METADATA_FILE,
      "run_accession",
      "All 194 assembly-QC-passing metagenomes in MATRIX_KEGG_ko.tsv",
      "none",
      "Aitchison (Euclidean distance in CLR space)",
      paste0("CLR; pseudocount=", PSEUDOCOUNT),
      paste0("Pearson Mantel; ", N_PERM, " sample-label permutations"),
      "Bray-Curtis on relative abundance",
      "1 - Bray-Curtis",
      "log10 geographic distance in km",
      paste0("paired two-sided sample-label permutation; ", N_PERM, " permutations"),
      as.character(MIN_GEO_KM),
      as.character(N_PERM),
      as.character(SENSITIVITY_REPS),
      as.character(SEED)
    )
  )
  write_tsv(parameters, file.path(TAB_DIR, "155_parameters.tsv"))

  run_summary <- data.table(
    metric = c(
      "functional_samples",
      "samples_with_coordinates",
      "coordinate_defined_sites",
      "taxonomic_features",
      "KOs",
      "panel_c_pairs",
      "panel_c_Mantel_r",
      "panel_c_Mantel_P",
      "panel_d_pairs_per_profile",
      "panel_d_taxonomic_slope",
      "panel_d_functional_slope",
      "panel_d_slope_difference_tax_minus_function",
      "panel_d_slope_difference_P_two_sided"
    ),
    value = c(
      n_samples,
      n_geo_samples,
      n_sites,
      n_taxa,
      n_kos,
      nrow(panel_c_pairs),
      unname(mantel_c$statistic),
      mantel_c$signif,
      slope_test$n_pairs,
      slope_test$slope_taxonomic,
      slope_test$slope_functional,
      slope_test$difference_tax_minus_function,
      slope_test$p_two_sided
    )
  )
  write_tsv(run_summary, file.path(TAB_DIR, "155_run_summary.tsv"))

  saveRDS(
    list(
      cohort = cohort,
      taxonomic_counts = tax_mat,
      functional_KO_counts = ko_mat,
      taxonomic_CLR = tax_clr,
      functional_CLR = func_clr,
      taxonomic_Aitchison = dist_tax_ait,
      functional_Aitchison = dist_func_ait,
      taxonomic_Bray_Curtis = dist_tax_bray_all,
      functional_Bray_Curtis = dist_func_bray_all,
      geographic_distance_km = geo_mat,
      panel_c_Mantel = mantel_c,
      panel_c_model = model_c,
      panel_d_slope_test = slope_test,
      site_sensitivity = sensitivity,
      parameters = parameters,
      run_summary = run_summary
    ),
    file.path(RDS_DIR, "155_Figure_4cd_analysis_objects.rds")
  )

  capture.output(sessionInfo(), file = file.path(OUTDIR, "sessionInfo.txt"))

  log_msg("Panel c Mantel r: ", sprintf("%.6f", unname(mantel_c$statistic)),
          "; P = ", mantel_c$signif)
  log_msg("Panel d taxonomic slope: ", sprintf("%.6f", slope_test$slope_taxonomic))
  log_msg("Panel d functional slope: ", sprintf("%.6f", slope_test$slope_functional))
  log_msg("Panel d slope difference: ",
          sprintf("%.6f", slope_test$difference_tax_minus_function),
          "; P = ", slope_test$p_two_sided)
  log_msg("Figure c: ", paste0(figure_c_base, ".png"))
  log_msg("Figure d: ", paste0(figure_d_base, ".png"))

  write_status(
    "PASS",
    c(
      paste0("functional_samples=", n_samples),
      paste0("samples_with_coordinates=", n_geo_samples),
      paste0("coordinate_defined_sites=", n_sites),
      paste0("panel_c_mantel_r=", unname(mantel_c$statistic)),
      paste0("panel_c_mantel_p=", mantel_c$signif),
      paste0("panel_d_taxonomic_slope=", slope_test$slope_taxonomic),
      paste0("panel_d_functional_slope=", slope_test$slope_functional),
      paste0("panel_d_slope_difference_p=", slope_test$p_two_sided)
    )
  )

  log_msg("=== SCRIPT COMPLETED SUCCESSFULLY ===")
  TRUE
}, error = function(e) {
  log_msg("ERROR: ", conditionMessage(e))
  write_status("FAIL", paste0("error=", conditionMessage(e)))
  stop(e)
})

invisible(result)
