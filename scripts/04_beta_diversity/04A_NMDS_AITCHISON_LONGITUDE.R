#!/usr/bin/env Rscript
# -*- coding: UTF-8 -*-

# Figure 2d — Global taxonomic composition
# Aitchison NMDS with marine-realm structure and longitude isolines.

options(stringsAsFactors = FALSE, scipen = 999, bitmapType = "cairo")

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(mgcv)
  library(ggplot2)
  library(ggnewscale)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(data.table)
  library(scales)
})

## -----------------------------------------------------------------------------
## 1. Canonical configuration
## -----------------------------------------------------------------------------

SEED <- 12345L
RAREFACTION_DEPTH <- 18145L
CLR_PSEUDOCOUNT <- 1
NMDS_K <- 2L
NMDS_TRY <- 50L
NMDS_TRYMAX <- 500L
LONGITUDE_GAM_K <- 30L
GRID_N <- 220L

EXPECTED_PHYLOSEQ_SAMPLES <- 401L
EXPECTED_COHORT_SAMPLES <- 332L
EXPECTED_ANALYSIS_SAMPLES <- 314L
EXPECTED_COORDINATE_SITES <- 130L
EXPECTED_COHORT_TAXA_NONZERO <- 24030L

PROJECT_DIR <- Sys.getenv("BIOGEOGRAPHY_BASE", unset = "/home/fjbalvino/Biogeografia_2026")
RESULTS_DIR <- file.path(PROJECT_DIR, "resultados")

PHYLOSEQ_PATH <- Sys.getenv("MAGENTA_PHYLOSEQ_401", unset = "/data/Magenta/phyloseq_MAGENTA.rds")
METADATA_PATH <- Sys.getenv(
  "MAGENTA_METADATA_CSV",
  unset = "/data/Magenta/data/metadata/Metadata_2026_environmental_reproduced.csv"
)
COHORT_POINTER <- file.path(
  RESULTS_DIR,
  "LATEST_04_ALPHA_DIVERSITY_HILL_SET18_PANEL_A.txt"
)
COHORT_FILENAME <- "04a_alpha_diversity_filtered_samplewise.csv"

SCRIPT_ID <- "05"
SCRIPT_BASENAME <- "NMDS_AITCHISON_PANEL_D_SET18"
RUN_STAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
RUN_NAME <- paste(SCRIPT_ID, SCRIPT_BASENAME, RUN_STAMP, sep = "_")

OUT_DIR <- file.path(RESULTS_DIR, RUN_NAME)
TABLE_DIR <- file.path(OUT_DIR, "tables")
PLOT_DIR <- file.path(OUT_DIR, "plots")
RDS_DIR <- file.path(OUT_DIR, "rds")
LOG_DIR <- file.path(OUT_DIR, "logs")

invisible(lapply(
  c(OUT_DIR, TABLE_DIR, PLOT_DIR, RDS_DIR, LOG_DIR),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

LATEST_POINTER <- file.path(
  RESULTS_DIR,
  paste0("LATEST_", SCRIPT_ID, "_", SCRIPT_BASENAME, ".txt")
)
writeLines(normalizePath(OUT_DIR, mustWork = FALSE), LATEST_POINTER)

LOG_FILE <- file.path(LOG_DIR, "05_panel_d.log")
sink(LOG_FILE, split = TRUE)
on.exit({
  if (sink.number() > 0L) sink()
}, add = TRUE)

log_msg <- function(...) {
  cat(sprintf(
    "[%s] %s\n",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    paste(..., collapse = "")
  ))
}

stopf <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)

resolve_required_column <- function(df, candidates, label) {
  nms <- names(df)
  hit <- match(tolower(candidates), tolower(nms), nomatch = 0L)
  hit <- hit[hit > 0L]
  if (length(hit) == 0L) {
    stopf(
      "No se encontró la columna %s. Candidatas: %s",
      label,
      paste(candidates, collapse = ", ")
    )
  }
  nms[hit[1L]]
}

resolve_optional_column <- function(df, candidates) {
  nms <- names(df)
  hit <- match(tolower(candidates), tolower(nms), nomatch = 0L)
  hit <- hit[hit > 0L]
  if (length(hit) == 0L) NA_character_ else nms[hit[1L]]
}

clean_character <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | !nzchar(x) | tolower(x) %in% c("na", "nan", "null")] <- NA_character_
  x
}

first_non_missing_character <- function(x) {
  x <- clean_character(x)
  x <- x[!is.na(x)]
  if (length(x) == 0L) NA_character_ else x[1L]
}

first_non_missing_numeric <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else x[1L]
}

n_unique_character <- function(x) {
  x <- clean_character(x)
  dplyr::n_distinct(x[!is.na(x)])
}

n_unique_numeric <- function(x, digits = 8L) {
  x <- x[is.finite(x)]
  dplyr::n_distinct(round(x, digits))
}

normalize_realm <- function(x) {
  raw <- clean_character(x)
  key <- gsub("[^a-z]+", " ", tolower(raw))
  key <- trimws(gsub(" +", " ", key))

  out <- rep(NA_character_, length(key))
  out[key %in% c("cip", "central indo pacific")] <- "Central Indo-Pacific"
  out[key %in% c("ip", "indo pacific")] <- "Indo-Pacific"
  out[key %in% c("tna", "temperate northern atlantic")] <-
    "Temperate Northern Atlantic"
  out[key %in% c("ta", "tropical atlantic")] <- "Tropical Atlantic"
  out[key %in% c("tep", "tropical eastern pacific")] <-
    "Tropical Eastern Pacific"
  out[key %in% c("wip", "western indo pacific")] <- "Western Indo-Pacific"
  out
}

point_in_polygon <- function(px, py, poly_x, poly_y) {
  n <- length(poly_x)
  if (n < 3L) return(rep(FALSE, length(px)))

  inside <- rep(FALSE, length(px))
  j <- n
  for (i in seq_len(n)) {
    crosses_y <- (poly_y[i] > py) != (poly_y[j] > py)
    denom <- poly_y[j] - poly_y[i]
    denom[abs(denom) < .Machine$double.eps] <- .Machine$double.eps
    x_cross <- (poly_x[j] - poly_x[i]) * (py - poly_y[i]) / denom + poly_x[i]
    inside <- xor(inside, crosses_y & px < x_cross)
    j <- i
  }
  inside
}

contour_lines_to_df <- function(contours) {
  if (length(contours) == 0L) {
    return(data.frame(
      NMDS1 = numeric(),
      NMDS2 = numeric(),
      longitude = numeric(),
      contour_group = character()
    ))
  }

  dplyr::bind_rows(lapply(seq_along(contours), function(i) {
    data.frame(
      NMDS1 = contours[[i]]$x,
      NMDS2 = contours[[i]]$y,
      longitude = contours[[i]]$level,
      contour_group = paste0("contour_", i),
      stringsAsFactors = FALSE
    )
  }))
}

log_msg("Starting ", RUN_NAME)
log_msg("Phyloseq: ", PHYLOSEQ_PATH)
log_msg("Metadata: ", METADATA_PATH)

## -----------------------------------------------------------------------------
## 2. Resolve the canonical Guadeloupe Set #18 cohort
## -----------------------------------------------------------------------------

required_files <- c(PHYLOSEQ_PATH, METADATA_PATH, COHORT_POINTER)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stopf("Faltan inputs requeridos: %s", paste(missing_files, collapse = "; "))
}

cohort_dir <- trimws(readLines(COHORT_POINTER, warn = FALSE, n = 1L))
if (!nzchar(cohort_dir) || !dir.exists(cohort_dir)) {
  stopf("Puntero de cohorte inválido: %s", COHORT_POINTER)
}

COHORT_PATH <- file.path(cohort_dir, COHORT_FILENAME)
if (!file.exists(COHORT_PATH)) {
  stopf("No existe la tabla maestra de la cohorte: %s", COHORT_PATH)
}

cohort <- data.table::fread(COHORT_PATH, data.table = FALSE, check.names = FALSE)
if (!"run_accession" %in% names(cohort)) {
  stopf("La tabla maestra no contiene run_accession.")
}

cohort$run_accession <- clean_character(cohort$run_accession)
if (anyNA(cohort$run_accession) || anyDuplicated(cohort$run_accession)) {
  stopf("La tabla maestra contiene run_accession faltantes o duplicados.")
}

cohort_ids <- cohort$run_accession
if (length(cohort_ids) != EXPECTED_COHORT_SAMPLES) {
  stopf(
    "La cohorte debe contener %d muestras; se encontraron %d.",
    EXPECTED_COHORT_SAMPLES,
    length(cohort_ids)
  )
}

## -----------------------------------------------------------------------------
## 3. Load phyloseq and construct the canonical count matrix
## -----------------------------------------------------------------------------

phy0 <- readRDS(PHYLOSEQ_PATH)
if (!inherits(phy0, "phyloseq")) stopf("El RDS no contiene un objeto phyloseq.")
if (nsamples(phy0) != EXPECTED_PHYLOSEQ_SAMPLES) {
  stopf(
    "El phyloseq debe contener %d muestras; contiene %d.",
    EXPECTED_PHYLOSEQ_SAMPLES,
    nsamples(phy0)
  )
}
if (anyDuplicated(sample_names(phy0))) stopf("sample_names contiene duplicados.")

missing_in_phyloseq <- setdiff(cohort_ids, sample_names(phy0))
if (length(missing_in_phyloseq) > 0L) {
  data.table::fwrite(
    data.frame(run_accession = missing_in_phyloseq),
    file.path(TABLE_DIR, "05_cohort_ids_missing_in_phyloseq.csv")
  )
  stopf("Hay %d accesiones ausentes en el phyloseq.", length(missing_in_phyloseq))
}

otu_counts <- as(otu_table(phy0), "matrix")
if (taxa_are_rows(phy0)) otu_counts <- t(otu_counts)
otu_counts <- otu_counts[cohort_ids, , drop = FALSE]
otu_counts <- otu_counts[, colSums(otu_counts) > 0, drop = FALSE]
cohort_taxa_nonzero <- ncol(otu_counts)

if (cohort_taxa_nonzero != EXPECTED_COHORT_TAXA_NONZERO) {
  stopf(
    "La cohorte debe contener %d taxa no nulos; contiene %d.",
    EXPECTED_COHORT_TAXA_NONZERO,
    cohort_taxa_nonzero
  )
}
if (any(!is.finite(otu_counts)) || any(otu_counts < 0)) {
  stopf("La matriz de conteos contiene valores no finitos o negativos.")
}
if (max(abs(otu_counts - round(otu_counts))) > 1e-8) {
  stopf("La matriz no contiene conteos enteros y no puede rarefactarse.")
}
storage.mode(otu_counts) <- "integer"

## -----------------------------------------------------------------------------
## 4. Curate current metadata without arbitrary duplicate removal
## -----------------------------------------------------------------------------

meta_raw <- data.table::fread(
  METADATA_PATH,
  data.table = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA", "NaN", "NULL")
)
if (!"run_accession" %in% names(meta_raw)) {
  stopf("La metadata no contiene run_accession.")
}

lat_col <- resolve_required_column(
  meta_raw,
  c("lat", "latitude", "decimalLatitude", "latitude_dd"),
  "de latitud"
)
lon_col <- resolve_required_column(
  meta_raw,
  c("lon", "longitude", "long", "decimalLongitude", "longitude_dd"),
  "de longitud"
)
realm_col <- resolve_required_column(meta_raw, c("REALM", "realm"), "REALM")

province_col <- resolve_optional_column(meta_raw, c("PROVINCE", "province"))
ecoregion_col <- resolve_optional_column(
  meta_raw,
  c("Ecoregion", "ECOREGION", "ecoregion")
)
country_col <- resolve_optional_column(meta_raw, c("Country", "COUNTRY", "country"))
study_col <- resolve_optional_column(
  meta_raw,
  c("study_accession", "BioProject", "bioproject", "study", "dataset")
)

optional_character <- function(df, column_name) {
  if (is.na(column_name)) rep(NA_character_, nrow(df)) else clean_character(df[[column_name]])
}

meta_work <- data.frame(
  run_accession = clean_character(meta_raw$run_accession),
  latitude = suppressWarnings(as.numeric(meta_raw[[lat_col]])),
  longitude = suppressWarnings(as.numeric(meta_raw[[lon_col]])),
  REALM_raw = clean_character(meta_raw[[realm_col]]),
  PROVINCE = optional_character(meta_raw, province_col),
  Ecoregion = optional_character(meta_raw, ecoregion_col),
  Country = optional_character(meta_raw, country_col),
  study_id = optional_character(meta_raw, study_col),
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(
    !is.na(run_accession),
    run_accession %in% cohort_ids
  )

metadata_audit <- meta_work %>%
  dplyr::group_by(run_accession) %>%
  dplyr::summarise(
    n_metadata_rows = dplyr::n(),
    n_unique_latitude = n_unique_numeric(latitude),
    n_unique_longitude = n_unique_numeric(longitude),
    n_unique_realm = n_unique_character(REALM_raw),
    n_unique_province = n_unique_character(PROVINCE),
    n_unique_ecoregion = n_unique_character(Ecoregion),
    n_unique_country = n_unique_character(Country),
    n_unique_study = n_unique_character(study_id),
    .groups = "drop"
  )

data.table::fwrite(
  metadata_audit,
  file.path(TABLE_DIR, "05_metadata_duplicate_audit.csv")
)

identity_conflicts <- metadata_audit %>%
  dplyr::filter(
    n_unique_latitude > 1L |
      n_unique_longitude > 1L |
      n_unique_realm > 1L |
      n_unique_province > 1L |
      n_unique_ecoregion > 1L |
      n_unique_country > 1L |
      n_unique_study > 1L
  )

if (nrow(identity_conflicts) > 0L) {
  data.table::fwrite(
    identity_conflicts,
    file.path(TABLE_DIR, "05_metadata_identity_conflicts_STOP.csv")
  )
  stopf(
    paste0(
      "Hay %d run_accession con conflictos de coordenadas o identidad regional. ",
      "Se detuvo el análisis."
    ),
    nrow(identity_conflicts)
  )
}

meta_one <- meta_work %>%
  dplyr::group_by(run_accession) %>%
  dplyr::summarise(
    latitude = first_non_missing_numeric(latitude),
    longitude = first_non_missing_numeric(longitude),
    REALM_raw = first_non_missing_character(REALM_raw),
    PROVINCE = first_non_missing_character(PROVINCE),
    Ecoregion = first_non_missing_character(Ecoregion),
    Country = first_non_missing_character(Country),
    study_id = first_non_missing_character(study_id),
    .groups = "drop"
  ) %>%
  dplyr::mutate(REALM = normalize_realm(REALM_raw))

unrecognized_realms <- meta_one %>%
  dplyr::filter(!is.na(REALM_raw), is.na(REALM))
if (nrow(unrecognized_realms) > 0L) {
  data.table::fwrite(
    unrecognized_realms,
    file.path(TABLE_DIR, "05_unrecognized_REALM_values_STOP.csv")
  )
  stopf("Hay %d valores REALM no reconocidos.", nrow(unrecognized_realms))
}

sample_table <- data.frame(run_accession = cohort_ids) %>%
  dplyr::left_join(meta_one, by = "run_accession") %>%
  dplyr::mutate(
    coordinate_valid = is.finite(latitude) & is.finite(longitude) &
      latitude >= -90 & latitude <= 90 &
      longitude >= -180 & longitude <= 180,
    realm_valid = !is.na(REALM)
  )

excluded_samples <- sample_table %>%
  dplyr::filter(!coordinate_valid | !realm_valid)
if (nrow(excluded_samples) > 0L) {
  data.table::fwrite(
    excluded_samples,
    file.path(TABLE_DIR, "05_samples_excluded_coordinates_or_REALM.csv")
  )
}

sample_table <- sample_table %>%
  dplyr::filter(coordinate_valid, realm_valid) %>%
  dplyr::mutate(
    site_id = paste0(
      sprintf("%.5f", latitude),
      "_",
      sprintf("%.5f", longitude)
    )
  )

analysis_ids <- sample_table$run_accession
if (length(analysis_ids) != EXPECTED_ANALYSIS_SAMPLES) {
  stopf(
    "El panel D debe usar %d muestras con coordenadas y REALM; se encontraron %d.",
    EXPECTED_ANALYSIS_SAMPLES,
    length(analysis_ids)
  )
}
if (dplyr::n_distinct(sample_table$site_id) != EXPECTED_COORDINATE_SITES) {
  stopf(
    "Se esperaban %d sitios definidos por coordenadas; se encontraron %d.",
    EXPECTED_COORDINATE_SITES,
    dplyr::n_distinct(sample_table$site_id)
  )
}

realm_levels <- c(
  "Central Indo-Pacific",
  "Indo-Pacific",
  "Temperate Northern Atlantic",
  "Tropical Atlantic",
  "Tropical Eastern Pacific",
  "Western Indo-Pacific"
)
sample_table$REALM <- factor(sample_table$REALM, levels = realm_levels)
sample_table <- sample_table[match(analysis_ids, sample_table$run_accession), , drop = FALSE]

data.table::fwrite(
  sample_table,
  file.path(TABLE_DIR, "05_samples_used_Figure_2d.csv")
)

realm_counts <- sample_table %>%
  dplyr::count(REALM, name = "n_metagenomes", .drop = FALSE) %>%
  dplyr::left_join(
    sample_table %>%
      dplyr::distinct(REALM, site_id) %>%
      dplyr::count(REALM, name = "n_coordinate_sites", .drop = FALSE),
    by = "REALM"
  )
data.table::fwrite(realm_counts, file.path(TABLE_DIR, "05_REALM_counts.csv"))

## -----------------------------------------------------------------------------
## 5. Rarefaction, CLR transformation, and distance matrices
## -----------------------------------------------------------------------------

otu_counts <- otu_counts[analysis_ids, , drop = FALSE]
sample_depths <- rowSums(otu_counts)
below_depth <- names(sample_depths)[sample_depths < RAREFACTION_DEPTH]
if (length(below_depth) > 0L) {
  data.table::fwrite(
    data.frame(run_accession = below_depth, library_size = sample_depths[below_depth]),
    file.path(TABLE_DIR, "05_samples_below_rarefaction_depth_STOP.csv")
  )
  stopf("Hay %d muestras por debajo de %d lecturas.", length(below_depth), RAREFACTION_DEPTH)
}

set.seed(SEED)
otu_rarefied <- vegan::rrarefy(otu_counts, sample = RAREFACTION_DEPTH)
otu_rarefied <- otu_rarefied[, colSums(otu_rarefied) > 0, drop = FALSE]
if (any(rowSums(otu_rarefied) != RAREFACTION_DEPTH)) {
  stopf("La rarefacción no produjo exactamente %d lecturas por muestra.", RAREFACTION_DEPTH)
}

otu_relative <- sweep(otu_rarefied, 1L, rowSums(otu_rarefied), "/")
log_counts <- log(otu_rarefied + CLR_PSEUDOCOUNT)
clr_matrix <- log_counts - rowMeans(log_counts)

dist_bray <- vegan::vegdist(otu_relative, method = "bray")
dist_aitchison <- stats::dist(clr_matrix, method = "euclidean")

## -----------------------------------------------------------------------------
## 6. Aitchison NMDS and deterministic axis orientation
## -----------------------------------------------------------------------------

log_msg("Running Aitchison NMDS with trymax = ", NMDS_TRYMAX)
set.seed(SEED + 10L)
nmds <- vegan::metaMDS(
  dist_aitchison,
  k = NMDS_K,
  try = NMDS_TRY,
  trymax = NMDS_TRYMAX,
  autotransform = FALSE,
  noshare = FALSE,
  wascores = FALSE,
  trace = FALSE
)

site_scores <- as.data.frame(vegan::scores(nmds, display = "sites"))
site_scores$run_accession <- rownames(site_scores)
site_scores <- site_scores[match(analysis_ids, site_scores$run_accession), , drop = FALSE]

cor_axis1_lon <- stats::cor(
  site_scores$NMDS1,
  sample_table$longitude,
  use = "complete.obs"
)
axis1_sign <- if (is.finite(cor_axis1_lon) && cor_axis1_lon > 0) -1 else 1

cor_axis2_lat <- stats::cor(
  site_scores$NMDS2,
  sample_table$latitude,
  use = "complete.obs"
)
axis2_sign <- if (is.finite(cor_axis2_lat) && cor_axis2_lat < 0) -1 else 1

site_scores$NMDS1 <- site_scores$NMDS1 * axis1_sign
site_scores$NMDS2 <- site_scores$NMDS2 * axis2_sign

plot_data <- sample_table %>%
  dplyr::left_join(site_scores, by = "run_accession")

if (anyNA(plot_data$NMDS1) || anyNA(plot_data$NMDS2)) {
  stopf("Falló el alineamiento entre scores NMDS y run_accession.")
}

## -----------------------------------------------------------------------------
## 7. Longitude surface masked to the observed NMDS convex hull
## -----------------------------------------------------------------------------

gam_k_used <- min(
  LONGITUDE_GAM_K,
  max(10L, floor(nrow(plot_data) / 5L))
)

fit_longitude <- mgcv::gam(
  longitude ~ s(NMDS1, NMDS2, bs = "tp", k = gam_k_used),
  data = plot_data,
  method = "REML",
  na.action = na.exclude
)

gx <- seq(min(plot_data$NMDS1), max(plot_data$NMDS1), length.out = GRID_N)
gy <- seq(min(plot_data$NMDS2), max(plot_data$NMDS2), length.out = GRID_N)
surface_grid <- expand.grid(NMDS1 = gx, NMDS2 = gy)
surface_grid$longitude_predicted <- as.numeric(stats::predict(
  fit_longitude,
  newdata = surface_grid,
  type = "response"
))

hull_index <- grDevices::chull(plot_data$NMDS1, plot_data$NMDS2)
hull_x <- plot_data$NMDS1[hull_index]
hull_y <- plot_data$NMDS2[hull_index]
surface_grid$inside_hull <- point_in_polygon(
  surface_grid$NMDS1,
  surface_grid$NMDS2,
  hull_x,
  hull_y
)
surface_grid$longitude_predicted[!surface_grid$inside_hull] <- NA_real_

z_matrix <- matrix(
  surface_grid$longitude_predicted,
  nrow = length(gx),
  ncol = length(gy)
)
contour_levels <- pretty(
  range(surface_grid$longitude_predicted, na.rm = TRUE),
  n = 7L
)
contour_levels <- contour_levels[
  contour_levels >= min(surface_grid$longitude_predicted, na.rm = TRUE) &
    contour_levels <= max(surface_grid$longitude_predicted, na.rm = TRUE)
]

contours <- grDevices::contourLines(
  x = gx,
  y = gy,
  z = z_matrix,
  levels = contour_levels
)
contour_data <- contour_lines_to_df(contours)

contour_labels <- contour_data %>%
  dplyr::group_by(longitude) %>%
  dplyr::slice_max(order_by = NMDS2, n = 1L, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(label = sprintf("%g", longitude))

data.table::fwrite(
  surface_grid,
  file.path(TABLE_DIR, "05_longitude_surface_grid.csv.gz"),
  compress = "gzip"
)
data.table::fwrite(
  contour_data,
  file.path(TABLE_DIR, "05_longitude_contour_lines.csv")
)

longitude_gam_statistics <- data.frame(
  response = "longitude",
  n = stats::nobs(fit_longitude),
  k_used = gam_k_used,
  edf = summary(fit_longitude)$s.table[1L, "edf"],
  ref_df = summary(fit_longitude)$s.table[1L, "Ref.df"],
  F = summary(fit_longitude)$s.table[1L, "F"],
  p_value = summary(fit_longitude)$s.table[1L, "p-value"],
  deviance_explained = summary(fit_longitude)$dev.expl,
  adjusted_r_squared = summary(fit_longitude)$r.sq
)
data.table::fwrite(
  longitude_gam_statistics,
  file.path(TABLE_DIR, "05_longitude_GAM_statistics.csv")
)

writeLines(
  capture.output(mgcv::k.check(fit_longitude)),
  file.path(LOG_DIR, "05_longitude_GAM_k_check.txt")
)

## -----------------------------------------------------------------------------
## 8. Publication panel
## -----------------------------------------------------------------------------

realm_colors <- c(
  "Central Indo-Pacific" = "#1B9E77",
  "Indo-Pacific" = "#E66101",
  "Temperate Northern Atlantic" = "#4C78A8",
  "Tropical Atlantic" = "#D95F9F",
  "Tropical Eastern Pacific" = "#66A61E",
  "Western Indo-Pacific" = "#E6AB02"
)

stress_label <- sprintf("Aitchison NMDS; stress = %.3f", nmds$stress)

p_panel_d <- ggplot() +
  geom_path(
    data = contour_data,
    aes(
      x = NMDS1,
      y = NMDS2,
      group = contour_group,
      color = longitude
    ),
    linewidth = 0.45,
    alpha = 0.80
  ) +
  geom_label(
    data = contour_labels,
    aes(x = NMDS1, y = NMDS2, label = label, color = longitude),
    size = 2.6,
    linewidth = 0,
    label.padding = grid::unit(0.06, "lines"),
    fill = alpha("white", 0.85),
    show.legend = FALSE
  ) +
  scale_color_viridis_c(
    name = "Longitude (°)",
    option = "plasma",
    direction = -1,
    guide = "none"
  ) +
  ggnewscale::new_scale_color() +
  geom_point(
    data = plot_data,
    aes(x = NMDS1, y = NMDS2, fill = REALM, color = REALM),
    shape = 21,
    size = 2.15,
    stroke = 0.38,
    alpha = 0.78
  ) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = stress_label,
    hjust = -0.02,
    vjust = 1.25,
    size = 3.4
  ) +
  scale_color_manual(values = realm_colors, drop = FALSE, guide = "none") +
  scale_fill_manual(
    values = realm_colors,
    drop = FALSE,
    name = "Marine realm"
  ) +
  labs(
    x = "NMDS1",
    y = "NMDS2"
  ) +
  coord_equal(clip = "off") +
  theme_classic(base_size = 11) +
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9.5, color = "black"),
    axis.line = element_line(linewidth = 0.45, color = "black"),
    axis.ticks = element_line(linewidth = 0.4, color = "black"),
    legend.title = element_text(face = "bold", size = 9.5),
    legend.text = element_text(size = 8.5),
    legend.key.height = grid::unit(0.42, "cm"),
    plot.margin = margin(t = 4, r = 7, b = 5, l = 7)
  )

ggsave(
  file.path(PLOT_DIR, "Figure_2d_global_taxonomic_composition.png"),
  p_panel_d,
  width = 7.4,
  height = 3.35,
  units = "in",
  dpi = 600,
  bg = "white"
)
ggsave(
  file.path(PLOT_DIR, "Figure_2d_global_taxonomic_composition.pdf"),
  p_panel_d,
  width = 7.4,
  height = 3.35,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)
ggsave(
  file.path(PLOT_DIR, "Figure_2d_global_taxonomic_composition.tiff"),
  p_panel_d,
  width = 7.4,
  height = 3.35,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

source_data <- plot_data %>%
  dplyr::select(
    run_accession,
    site_id,
    REALM,
    PROVINCE,
    Ecoregion,
    Country,
    study_id,
    latitude,
    longitude,
    NMDS1,
    NMDS2
  )
data.table::fwrite(
  source_data,
  file.path(TABLE_DIR, "05_source_data_Figure_2d.csv")
)

## -----------------------------------------------------------------------------
## 9. Save synchronized objects and provenance
## -----------------------------------------------------------------------------

nmds_statistics <- data.frame(
  distance = "Aitchison (Euclidean distance on CLR)",
  pseudocount = CLR_PSEUDOCOUNT,
  k = NMDS_K,
  stress = nmds$stress,
  tries = if (!is.null(nmds$tries)) nmds$tries else NA_integer_,
  converged = if (!is.null(nmds$converged)) as.character(nmds$converged) else NA_character_,
  n_samples = nrow(plot_data),
  n_taxa_after_rarefaction = ncol(otu_rarefied),
  axis1_sign = axis1_sign,
  axis2_sign = axis2_sign
)
data.table::fwrite(nmds_statistics, file.path(TABLE_DIR, "05_NMDS_statistics.csv"))

run_summary <- data.frame(
  status = "PASS",
  run_name = RUN_NAME,
  phyloseq_path = PHYLOSEQ_PATH,
  metadata_path = METADATA_PATH,
  cohort_path = COHORT_PATH,
  phyloseq_samples = nsamples(phy0),
  balanced_cohort_samples = length(cohort_ids),
  analysis_samples = length(analysis_ids),
  coordinate_defined_sites = dplyr::n_distinct(sample_table$site_id),
  cohort_taxa_nonzero = cohort_taxa_nonzero,
  taxa_nonzero_after_rarefaction = ncol(otu_rarefied),
  rarefaction_depth = RAREFACTION_DEPTH,
  clr_pseudocount = CLR_PSEUDOCOUNT,
  NMDS_stress = nmds$stress,
  seed = SEED,
  latitude_column = lat_col,
  longitude_column = lon_col,
  realm_column = realm_col,
  province_column = ifelse(is.na(province_col), "NOT_FOUND", province_col),
  ecoregion_column = ifelse(is.na(ecoregion_col), "NOT_FOUND", ecoregion_col),
  country_column = ifelse(is.na(country_col), "NOT_FOUND", country_col),
  study_column = ifelse(is.na(study_col), "NOT_FOUND", study_col),
  stringsAsFactors = FALSE
)
data.table::fwrite(run_summary, file.path(TABLE_DIR, "05_run_summary.csv"))

analysis_bundle <- list(
  sample_table = sample_table,
  source_data = source_data,
  dist_aitchison = dist_aitchison,
  dist_bray = dist_bray,
  nmds = nmds,
  oriented_site_scores = site_scores,
  longitude_gam = fit_longitude,
  parameters = list(
    seed = SEED,
    rarefaction_depth = RAREFACTION_DEPTH,
    clr_pseudocount = CLR_PSEUDOCOUNT,
    cohort_path = COHORT_PATH,
    phyloseq_path = PHYLOSEQ_PATH,
    metadata_path = METADATA_PATH
  )
)

saveRDS(
  analysis_bundle,
  file.path(RDS_DIR, "05_panel_d_analysis_objects.rds")
)
saveRDS(p_panel_d, file.path(RDS_DIR, "05_Figure_2d_ggplot_object.rds"))

input_paths <- c(PHYLOSEQ_PATH, METADATA_PATH, COHORT_POINTER, COHORT_PATH)
input_provenance <- data.frame(
  input = c("phyloseq", "metadata", "cohort_pointer", "cohort_table"),
  path = input_paths,
  size_bytes = as.numeric(file.info(input_paths)$size),
  modified_time = as.character(file.info(input_paths)$mtime),
  stringsAsFactors = FALSE
)
data.table::fwrite(input_provenance, file.path(TABLE_DIR, "05_input_provenance.csv"))

writeLines(capture.output(sessionInfo()), file.path(LOG_DIR, "05_sessionInfo.txt"))

log_msg("Analysis samples: ", nrow(plot_data))
log_msg("Coordinate-defined sites: ", dplyr::n_distinct(plot_data$site_id))
log_msg("Aitchison NMDS stress: ", sprintf("%.6f", nmds$stress))
log_msg("Figure: ", file.path(PLOT_DIR, "Figure_2d_global_taxonomic_composition.png"))
log_msg("=== SCRIPT COMPLETED SUCCESSFULLY ===")
