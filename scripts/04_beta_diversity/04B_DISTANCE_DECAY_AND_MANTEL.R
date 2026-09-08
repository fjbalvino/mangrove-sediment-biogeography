#!/usr/bin/env Rscript
# -*- coding: UTF-8 -*-

# Figure 2c — Taxonomic distance–decay
# Abundance-weighted Bray–Curtis versus incidence-based Sørensen.
#
# Canonical inputs (2026-08-17):
#   - phyloseq: /data/Magenta/phyloseq_MAGENTA.rds (401 metagenomes)
#   - metadata: /data/Magenta/data/metadata/Metadata_2026_copia.csv
#   - balanced cohort: inherited from Figure 2a, Guadeloupe Set #18 (332 samples)
#
# Statistical design:
#   - deterministic rarefaction to 18,145 reads, matching Figure 2a;
#   - sample–sample geographic distances (Haversine, km);
#   - similarities regressed on log10 geographic distance;
#   - Mantel-type sample-label permutations for spatial inference;
#   - paired metric contrast assessed by permuting the Bray–Curtis minus
#     Sørensen similarity matrix against geographic distance;
#   - one-sample-per-coordinate sensitivity to reduce within-site replication;
#   - leave-one-study-out sensitivity when a study/BioProject field is present.

options(stringsAsFactors = FALSE, scipen = 999, bitmapType = "cairo")

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(geosphere)
  library(ggplot2)
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
N_PERM <- 9999L
N_SITE_BALANCED <- 1000L
RAREFACTION_DEPTH <- 18145L
MIN_GEO_KM <- 0.1

EXPECTED_PHYLOSEQ_SAMPLES <- 401L
EXPECTED_COHORT_SAMPLES <- 332L
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

SCRIPT_ID <- "06"
SCRIPT_BASENAME <- "MANTEL_DDP_ANALYSIS_PANEL_C_SET18"
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

LOG_FILE <- file.path(LOG_DIR, paste0(SCRIPT_ID, "_panel_c.log"))
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

stopf <- function(fmt, ...) {
  stop(sprintf(fmt, ...), call. = FALSE)
}

format_p <- function(p) {
  if (!is.finite(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  sprintf("= %.3f", p)
}

format_signed <- function(x, digits = 3L) {
  out <- sprintf(paste0("%.", digits, "f"), x)
  sub("-", "\u2212", out, fixed = TRUE)
}

resolve_required_column <- function(df, candidates, label) {
  nms <- names(df)
  hit <- match(tolower(candidates), tolower(nms), nomatch = 0L)
  hit <- hit[hit > 0L]
  if (length(hit) == 0L) {
    stopf(
      "No se encontró la columna %s. Candidatas evaluadas: %s",
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
  if (length(hit) == 0L) return(NA_character_)
  nms[hit[1L]]
}

first_non_missing_numeric <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else x[1L]
}

first_non_missing_character <- function(x) {
  x <- trimws(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0L) NA_character_ else x[1L]
}

two_sided_permutation_p <- function(mantel_result) {
  observed <- unname(mantel_result$statistic)
  permuted <- mantel_result$perm
  if (length(permuted) == 0L || !is.finite(observed)) return(NA_real_)
  (1 + sum(abs(permuted) >= abs(observed) - sqrt(.Machine$double.eps))) /
    (length(permuted) + 1)
}

linear_summary <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 3L || stats::var(x) <= 0 || stats::var(y) <= 0) {
    return(c(slope = NA_real_, r_squared = NA_real_))
  }
  slope <- stats::cov(x, y) / stats::var(x)
  r_squared <- stats::cor(x, y)^2
  c(slope = unname(slope), r_squared = unname(r_squared))
}

extract_distance_decay_stats <- function(indices, bray_matrix, sorensen_matrix,
                                         geo_matrix, min_geo_km) {
  g <- geo_matrix[indices, indices, drop = FALSE]
  b <- bray_matrix[indices, indices, drop = FALSE]
  s <- sorensen_matrix[indices, indices, drop = FALSE]

  keep <- upper.tri(g) & is.finite(g) & g > min_geo_km
  x <- log10(g[keep])
  y_bray <- 1 - b[keep]
  y_sorensen <- 1 - s[keep]

  fit_bray <- linear_summary(x, y_bray)
  fit_sorensen <- linear_summary(x, y_sorensen)

  data.frame(
    n_samples = length(indices),
    n_pairs = length(x),
    slope_bray_curtis = fit_bray[["slope"]],
    r2_bray_curtis = fit_bray[["r_squared"]],
    slope_sorensen = fit_sorensen[["slope"]],
    r2_sorensen = fit_sorensen[["r_squared"]],
    slope_difference_bray_minus_sorensen =
      fit_bray[["slope"]] - fit_sorensen[["slope"]]
  )
}

log_msg("Starting ", RUN_NAME)
log_msg("Canonical phyloseq: ", PHYLOSEQ_PATH)
log_msg("Canonical metadata: ", METADATA_PATH)

## -----------------------------------------------------------------------------
## 2. Resolve and validate the balanced Figure 2a cohort
## -----------------------------------------------------------------------------

required_files <- c(PHYLOSEQ_PATH, METADATA_PATH, COHORT_POINTER)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stopf("Faltan inputs requeridos: %s", paste(missing_files, collapse = "; "))
}

cohort_dir <- trimws(readLines(COHORT_POINTER, warn = FALSE, n = 1L))
if (!nzchar(cohort_dir) || !dir.exists(cohort_dir)) {
  stopf("El puntero de la cohorte no contiene un directorio válido: %s", COHORT_POINTER)
}

COHORT_PATH <- file.path(cohort_dir, COHORT_FILENAME)
if (!file.exists(COHORT_PATH)) {
  stopf("No se encontró la tabla maestra de la cohorte: %s", COHORT_PATH)
}

cohort <- data.table::fread(COHORT_PATH, data.table = FALSE, check.names = FALSE)
if (!"run_accession" %in% names(cohort)) {
  stopf("La tabla maestra no contiene la llave canónica run_accession.")
}

cohort$run_accession <- trimws(as.character(cohort$run_accession))
if (anyNA(cohort$run_accession) || any(!nzchar(cohort$run_accession))) {
  stopf("La tabla maestra contiene run_accession faltantes o vacíos.")
}
if (anyDuplicated(cohort$run_accession)) {
  duplicated_ids <- unique(cohort$run_accession[duplicated(cohort$run_accession)])
  stopf(
    "La tabla maestra contiene run_accession duplicados: %s",
    paste(head(duplicated_ids, 20L), collapse = ", ")
  )
}

cohort_ids <- cohort$run_accession
if (length(cohort_ids) != EXPECTED_COHORT_SAMPLES) {
  stopf(
    "La cohorte canónica debe contener %d muestras, pero contiene %d.",
    EXPECTED_COHORT_SAMPLES,
    length(cohort_ids)
  )
}

log_msg("Balanced cohort inherited from Figure 2a: n = ", length(cohort_ids))
log_msg("Cohort table: ", COHORT_PATH)

## -----------------------------------------------------------------------------
## 3. Load canonical phyloseq and align strictly by run_accession
## -----------------------------------------------------------------------------

phy0 <- readRDS(PHYLOSEQ_PATH)
if (!inherits(phy0, "phyloseq")) {
  stopf("El objeto leído no es de clase phyloseq: %s", PHYLOSEQ_PATH)
}
if (nsamples(phy0) != EXPECTED_PHYLOSEQ_SAMPLES) {
  stopf(
    "El phyloseq canónico debe contener %d muestras, pero contiene %d.",
    EXPECTED_PHYLOSEQ_SAMPLES,
    nsamples(phy0)
  )
}
if (anyDuplicated(sample_names(phy0))) {
  stopf("El phyloseq contiene sample_names duplicados.")
}

missing_in_phyloseq <- setdiff(cohort_ids, sample_names(phy0))
if (length(missing_in_phyloseq) > 0L) {
  data.table::fwrite(
    data.frame(run_accession = missing_in_phyloseq),
    file.path(TABLE_DIR, "06_cohort_ids_missing_in_phyloseq.csv")
  )
  stopf(
    "%d muestras de la cohorte no están presentes en el phyloseq canónico.",
    length(missing_in_phyloseq)
  )
}

otu_counts <- as(otu_table(phy0), "matrix")
if (taxa_are_rows(phy0)) otu_counts <- t(otu_counts)
otu_counts <- otu_counts[cohort_ids, , drop = FALSE]
otu_counts <- otu_counts[, colSums(otu_counts) > 0, drop = FALSE]
cohort_taxa_nonzero <- ncol(otu_counts)

if (nrow(otu_counts) != EXPECTED_COHORT_SAMPLES) {
  stopf("Falló la alineación de filas de la matriz OTU con la cohorte.")
}
if (ncol(otu_counts) != EXPECTED_COHORT_TAXA_NONZERO) {
  stopf(
    paste0(
      "La cohorte debe contener %d taxa no nulos, pero contiene %d. ",
      "Esto sugiere que cambió el objeto o la cohorte canónica."
    ),
    EXPECTED_COHORT_TAXA_NONZERO,
    ncol(otu_counts)
  )
}

if (any(!is.finite(otu_counts)) || any(otu_counts < 0)) {
  stopf("La matriz OTU contiene valores no finitos o negativos.")
}
if (max(abs(otu_counts - round(otu_counts))) > 1e-8) {
  stopf("La matriz OTU no contiene conteos enteros; no se puede rarefactar.")
}
storage.mode(otu_counts) <- "integer"

## -----------------------------------------------------------------------------
## 4. Curate current coordinates and audit duplicated metadata rows
## -----------------------------------------------------------------------------

meta_raw <- data.table::fread(
  METADATA_PATH,
  data.table = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA", "NaN", "NULL")
)
if (!"run_accession" %in% names(meta_raw)) {
  stopf("La metadata actual no contiene la llave canónica run_accession.")
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
study_col <- resolve_optional_column(
  meta_raw,
  c("study_accession", "BioProject", "bioproject", "study", "dataset")
)

meta_geo <- data.frame(
  run_accession = trimws(as.character(meta_raw$run_accession)),
  latitude = suppressWarnings(as.numeric(meta_raw[[lat_col]])),
  longitude = suppressWarnings(as.numeric(meta_raw[[lon_col]])),
  study_id = if (is.na(study_col)) {
    rep(NA_character_, nrow(meta_raw))
  } else {
    trimws(as.character(meta_raw[[study_col]]))
  },
  stringsAsFactors = FALSE
)
meta_geo <- meta_geo %>%
  filter(!is.na(run_accession), nzchar(run_accession), run_accession %in% cohort_ids)

metadata_audit <- meta_geo %>%
  group_by(run_accession) %>%
  summarise(
    n_metadata_rows = n(),
    n_unique_latitude = n_distinct(round(latitude[is.finite(latitude)], 8L)),
    n_unique_longitude = n_distinct(round(longitude[is.finite(longitude)], 8L)),
    n_unique_study = n_distinct(study_id[!is.na(study_id) & nzchar(study_id)]),
    .groups = "drop"
  )

data.table::fwrite(
  metadata_audit,
  file.path(TABLE_DIR, "06_metadata_duplicate_coordinate_study_audit.csv")
)

coordinate_conflicts <- metadata_audit %>%
  filter(n_unique_latitude > 1L | n_unique_longitude > 1L)
study_conflicts <- metadata_audit %>% filter(n_unique_study > 1L)

if (nrow(coordinate_conflicts) > 0L) {
  data.table::fwrite(
    coordinate_conflicts,
    file.path(TABLE_DIR, "06_coordinate_conflicts_STOP.csv")
  )
  stopf(
    paste0(
      "Hay %d run_accession con coordenadas conflictivas entre filas de metadata. ",
      "Se detuvo el análisis para no promediar coordenadas arbitrariamente."
    ),
    nrow(coordinate_conflicts)
  )
}
if (nrow(study_conflicts) > 0L) {
  data.table::fwrite(
    study_conflicts,
    file.path(TABLE_DIR, "06_study_conflicts_STOP.csv")
  )
  stopf(
    "%d run_accession tienen más de un estudio/BioProject asignado.",
    nrow(study_conflicts)
  )
}

meta_one <- meta_geo %>%
  group_by(run_accession) %>%
  summarise(
    latitude = first_non_missing_numeric(latitude),
    longitude = first_non_missing_numeric(longitude),
    study_id = first_non_missing_character(study_id),
    .groups = "drop"
  )

missing_metadata_ids <- setdiff(cohort_ids, meta_one$run_accession)
if (length(missing_metadata_ids) > 0L) {
  data.table::fwrite(
    data.frame(run_accession = missing_metadata_ids),
    file.path(TABLE_DIR, "06_cohort_ids_missing_in_current_metadata.csv")
  )
}

sample_table <- data.frame(run_accession = cohort_ids) %>%
  left_join(meta_one, by = "run_accession") %>%
  mutate(
    coordinate_complete = is.finite(latitude) & is.finite(longitude),
    coordinate_valid = coordinate_complete &
      latitude >= -90 & latitude <= 90 &
      longitude >= -180 & longitude <= 180
  )

invalid_coordinates <- sample_table %>% filter(!coordinate_valid)
if (nrow(invalid_coordinates) > 0L) {
  data.table::fwrite(
    invalid_coordinates,
    file.path(TABLE_DIR, "06_samples_excluded_missing_or_invalid_coordinates.csv")
  )
  log_msg(
    "Samples excluded for missing/invalid coordinates: ",
    nrow(invalid_coordinates)
  )
}

sample_table <- sample_table %>%
  filter(coordinate_valid) %>%
  mutate(
    site_id = paste0(
      sprintf("%.5f", latitude),
      "_",
      sprintf("%.5f", longitude)
    )
  )

analysis_ids <- sample_table$run_accession
otu_counts <- otu_counts[analysis_ids, , drop = FALSE]

if (length(analysis_ids) < 3L) {
  stopf("Quedaron menos de tres muestras con coordenadas válidas.")
}

data.table::fwrite(
  sample_table,
  file.path(TABLE_DIR, "06_samples_used_Figure_2c.csv")
)

log_msg("Samples with valid coordinates: n = ", length(analysis_ids))
log_msg("Coordinate-defined sampling sites: n = ", n_distinct(sample_table$site_id))

## -----------------------------------------------------------------------------
## 5. Rarefy counts and calculate Bray–Curtis and Sørensen dissimilarities
## -----------------------------------------------------------------------------

sample_depths <- rowSums(otu_counts)
insufficient_depth <- names(sample_depths)[sample_depths < RAREFACTION_DEPTH]
if (length(insufficient_depth) > 0L) {
  data.table::fwrite(
    data.frame(
      run_accession = insufficient_depth,
      library_size = sample_depths[insufficient_depth]
    ),
    file.path(TABLE_DIR, "06_samples_below_rarefaction_depth_STOP.csv")
  )
  stopf(
    "%d muestras están por debajo de la profundidad de rarefacción %d.",
    length(insufficient_depth),
    RAREFACTION_DEPTH
  )
}

set.seed(SEED)
otu_rarefied <- vegan::rrarefy(otu_counts, sample = RAREFACTION_DEPTH)
otu_rarefied <- otu_rarefied[, colSums(otu_rarefied) > 0, drop = FALSE]

if (any(rowSums(otu_rarefied) != RAREFACTION_DEPTH)) {
  stopf("La rarefacción no produjo exactamente %d lecturas por muestra.", RAREFACTION_DEPTH)
}

otu_relative <- sweep(otu_rarefied, 1L, rowSums(otu_rarefied), "/")

log_msg("Calculating Bray-Curtis and Sorensen dissimilarities...")
dist_bray <- vegan::vegdist(otu_relative, method = "bray")
dist_sorensen <- vegan::vegdist(otu_rarefied, method = "bray", binary = TRUE)

bray_matrix <- as.matrix(dist_bray)
sorensen_matrix <- as.matrix(dist_sorensen)

## -----------------------------------------------------------------------------
## 6. Geographic distance and pairwise source data
## -----------------------------------------------------------------------------

sample_table <- sample_table[match(analysis_ids, sample_table$run_accession), , drop = FALSE]
coords <- as.matrix(sample_table[, c("longitude", "latitude")])
rownames(coords) <- analysis_ids

geo_matrix <- geosphere::distm(coords, fun = geosphere::distHaversine) / 1000
dimnames(geo_matrix) <- list(analysis_ids, analysis_ids)

pair_index <- which(
  upper.tri(geo_matrix) & is.finite(geo_matrix) & geo_matrix > MIN_GEO_KM,
  arr.ind = TRUE
)
if (nrow(pair_index) < 3L) {
  stopf("Quedaron menos de tres pares con distancia geográfica > %.3f km.", MIN_GEO_KM)
}

pair_base <- data.frame(
  sample1 = rownames(geo_matrix)[pair_index[, 1L]],
  sample2 = colnames(geo_matrix)[pair_index[, 2L]],
  geo_km = geo_matrix[pair_index],
  log_geo = log10(geo_matrix[pair_index]),
  bray_dissimilarity = bray_matrix[pair_index],
  sorensen_dissimilarity = sorensen_matrix[pair_index],
  stringsAsFactors = FALSE
) %>%
  mutate(
    bray_similarity = 1 - bray_dissimilarity,
    sorensen_similarity = 1 - sorensen_dissimilarity,
    delta_similarity_bray_minus_sorensen =
      bray_similarity - sorensen_similarity
  )

pair_long <- bind_rows(
  pair_base %>%
    transmute(
      sample1,
      sample2,
      geo_km,
      log_geo,
      metric = "Bray–Curtis",
      dissimilarity = bray_dissimilarity,
      similarity = bray_similarity
    ),
  pair_base %>%
    transmute(
      sample1,
      sample2,
      geo_km,
      log_geo,
      metric = "Sørensen",
      dissimilarity = sorensen_dissimilarity,
      similarity = sorensen_similarity
    )
) %>%
  mutate(metric = factor(metric, levels = c("Bray–Curtis", "Sørensen")))

data.table::fwrite(
  pair_long,
  file.path(TABLE_DIR, "06_source_data_Figure_2c_pairwise.csv.gz"),
  compress = "gzip"
)

log_msg("Eligible sample pairs (> ", MIN_GEO_KM, " km): n = ", nrow(pair_base))

## -----------------------------------------------------------------------------
## 7. Regression slopes and sample-label permutation tests
## -----------------------------------------------------------------------------

fit_bray <- linear_summary(pair_base$log_geo, pair_base$bray_similarity)
fit_sorensen <- linear_summary(pair_base$log_geo, pair_base$sorensen_similarity)

geo_log_matrix <- matrix(
  NA_real_,
  nrow = nrow(geo_matrix),
  ncol = ncol(geo_matrix),
  dimnames = dimnames(geo_matrix)
)
eligible_geo <- is.finite(geo_matrix) & geo_matrix > MIN_GEO_KM
geo_log_matrix[eligible_geo] <- log10(geo_matrix[eligible_geo])
diag(geo_log_matrix) <- 0
dist_geo_log <- stats::as.dist(geo_log_matrix)

delta_similarity_matrix <- (1 - bray_matrix) - (1 - sorensen_matrix)
diag(delta_similarity_matrix) <- 0
dist_delta_similarity <- stats::as.dist(delta_similarity_matrix)

set.seed(SEED + 100L)
permutation_matrix <- permute::shuffleSet(
  n = length(analysis_ids),
  nset = N_PERM,
  control = permute::how(nperm = N_PERM)
)

log_msg("Running sample-label permutation tests: nperm = ", N_PERM)
mantel_bray <- vegan::mantel(
  dist_bray,
  dist_geo_log,
  method = "pearson",
  permutations = permutation_matrix,
  na.rm = TRUE,
  parallel = 1
)
mantel_sorensen <- vegan::mantel(
  dist_sorensen,
  dist_geo_log,
  method = "pearson",
  permutations = permutation_matrix,
  na.rm = TRUE,
  parallel = 1
)
mantel_metric_contrast <- vegan::mantel(
  dist_delta_similarity,
  dist_geo_log,
  method = "pearson",
  permutations = permutation_matrix,
  na.rm = TRUE,
  parallel = 1
)

p_bray <- two_sided_permutation_p(mantel_bray)
p_sorensen <- two_sided_permutation_p(mantel_sorensen)
p_metric_contrast <- two_sided_permutation_p(mantel_metric_contrast)

regression_statistics <- data.frame(
  metric = c("Bray-Curtis", "Sorensen"),
  community_basis = c("relative abundance", "presence-absence"),
  n_samples = length(analysis_ids),
  n_pairs = nrow(pair_base),
  rarefaction_depth = RAREFACTION_DEPTH,
  min_geographic_distance_km = MIN_GEO_KM,
  slope_similarity_per_log10_km = c(
    fit_bray[["slope"]],
    fit_sorensen[["slope"]]
  ),
  r_squared_descriptive = c(
    fit_bray[["r_squared"]],
    fit_sorensen[["r_squared"]]
  ),
  mantel_r_dissimilarity_vs_log10_distance = c(
    unname(mantel_bray$statistic),
    unname(mantel_sorensen$statistic)
  ),
  permutation_p_two_sided = c(p_bray, p_sorensen),
  n_permutations = N_PERM
)

metric_contrast_statistics <- data.frame(
  contrast = "Bray-Curtis minus Sorensen similarity slope",
  slope_difference = fit_bray[["slope"]] - fit_sorensen[["slope"]],
  mantel_r_delta_similarity_vs_log10_distance =
    unname(mantel_metric_contrast$statistic),
  permutation_p_two_sided = p_metric_contrast,
  n_permutations = N_PERM,
  interpretation = paste(
    "Two-sided sample-label permutation of the paired similarity contrast;",
    "ordinary pairwise LM P-values are not used."
  )
)

data.table::fwrite(
  regression_statistics,
  file.path(TABLE_DIR, "06_distance_decay_regression_and_permutation_statistics.csv")
)
data.table::fwrite(
  metric_contrast_statistics,
  file.path(TABLE_DIR, "06_bray_vs_sorensen_slope_contrast.csv")
)

## -----------------------------------------------------------------------------
## 8. Site-balanced and leave-one-study-out robustness analyses
## -----------------------------------------------------------------------------

site_groups <- split(seq_len(nrow(sample_table)), sample_table$site_id)
set.seed(SEED + 200L)
site_balanced_results <- bind_rows(lapply(seq_len(N_SITE_BALANCED), function(i) {
  selected <- vapply(site_groups, function(z) sample(z, 1L), integer(1L))
  out <- extract_distance_decay_stats(
    indices = selected,
    bray_matrix = bray_matrix,
    sorensen_matrix = sorensen_matrix,
    geo_matrix = geo_matrix,
    min_geo_km = MIN_GEO_KM
  )
  out$replicate <- i
  out
})) %>%
  select(replicate, everything())

site_balanced_summary <- site_balanced_results %>%
  select(
    slope_bray_curtis,
    slope_sorensen,
    slope_difference_bray_minus_sorensen
  ) %>%
  pivot_longer(everything(), names_to = "estimand", values_to = "estimate") %>%
  group_by(estimand) %>%
  summarise(
    n_replicates = sum(is.finite(estimate)),
    median = median(estimate, na.rm = TRUE),
    q025 = quantile(estimate, 0.025, na.rm = TRUE),
    q975 = quantile(estimate, 0.975, na.rm = TRUE),
    proportion_negative = mean(estimate < 0, na.rm = TRUE),
    .groups = "drop"
  )

data.table::fwrite(
  site_balanced_results,
  file.path(TABLE_DIR, "06_site_balanced_one_sample_per_coordinate_replicates.csv")
)
data.table::fwrite(
  site_balanced_summary,
  file.path(TABLE_DIR, "06_site_balanced_one_sample_per_coordinate_summary.csv")
)

study_values <- sort(unique(sample_table$study_id[
  !is.na(sample_table$study_id) & nzchar(sample_table$study_id)
]))

if (length(study_values) >= 2L) {
  leave_one_study_out <- bind_rows(lapply(study_values, function(study) {
    indices <- which(is.na(sample_table$study_id) | sample_table$study_id != study)
    out <- extract_distance_decay_stats(
      indices = indices,
      bray_matrix = bray_matrix,
      sorensen_matrix = sorensen_matrix,
      geo_matrix = geo_matrix,
      min_geo_km = MIN_GEO_KM
    )
    out$study_excluded <- study
    out$n_samples_excluded <- sum(sample_table$study_id == study, na.rm = TRUE)
    out
  })) %>%
    select(study_excluded, n_samples_excluded, everything())
} else {
  leave_one_study_out <- data.frame(
    study_excluded = character(),
    n_samples_excluded = integer(),
    n_samples = integer(),
    n_pairs = integer(),
    slope_bray_curtis = numeric(),
    r2_bray_curtis = numeric(),
    slope_sorensen = numeric(),
    r2_sorensen = numeric(),
    slope_difference_bray_minus_sorensen = numeric()
  )
  log_msg("Leave-one-study-out skipped: fewer than two resolved study levels.")
}

data.table::fwrite(
  leave_one_study_out,
  file.path(TABLE_DIR, "06_leave_one_study_out_sensitivity.csv")
)

## -----------------------------------------------------------------------------
## 9. Figure 2c
## -----------------------------------------------------------------------------

color_values <- c(
  "Bray–Curtis" = "#009E73",
  "Sørensen" = "#E64B00"
)

label_bray <- paste0(
  "Bray–Curtis: r = ",
  sprintf("%.3f", unname(mantel_bray$statistic)),
  ", slope = ",
  format_signed(fit_bray[["slope"]]),
  ", P ",
  format_p(p_bray)
)
label_sorensen <- paste0(
  "Sørensen: r = ",
  sprintf("%.3f", unname(mantel_sorensen$statistic)),
  ", slope = ",
  format_signed(fit_sorensen[["slope"]]),
  ", P ",
  format_p(p_sorensen)
)
label_contrast <- paste0(
  "Slope contrast: Δ = ",
  format_signed(fit_bray[["slope"]] - fit_sorensen[["slope"]]),
  ", P ",
  format_p(p_metric_contrast)
)

x_annotation <- unname(stats::quantile(pair_base$log_geo, 0.01, na.rm = TRUE))
x_break_values <- c(1, 10, 100, 1000, 10000)

p_ddp <- ggplot(
  pair_long,
  aes(x = log_geo, y = similarity, color = metric)
) +
  geom_point(size = 0.65, alpha = 0.17, stroke = 0) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    linewidth = 1.15
  ) +
  annotate(
    "text",
    x = x_annotation,
    y = 0.22,
    label = label_bray,
    hjust = 0,
    vjust = 0.5,
    size = 3.55,
    color = color_values[["Bray–Curtis"]]
  ) +
  annotate(
    "text",
    x = x_annotation,
    y = 0.14,
    label = label_sorensen,
    hjust = 0,
    vjust = 0.5,
    size = 3.55,
    color = color_values[["Sørensen"]]
  ) +
  annotate(
    "text",
    x = x_annotation,
    y = 0.06,
    label = label_contrast,
    hjust = 0,
    vjust = 0.5,
    size = 3.55,
    color = "black"
  ) +
  scale_color_manual(values = color_values, drop = FALSE) +
  scale_x_continuous(
    name = expression(paste("Geographic distance (km, ", log[10], " scale)")),
    breaks = log10(x_break_values),
    labels = scales::label_comma()(x_break_values),
    expand = expansion(mult = c(0.025, 0.025))
  ) +
  scale_y_continuous(
    name = "Community similarity (1 − dissimilarity)",
    breaks = seq(0, 1, by = 0.25),
    labels = scales::label_number(accuracy = 0.01),
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.015))
  ) +
  labs(
    tag = "c",
    title = "Taxonomic distance–decay"
  ) +
  guides(color = "none") +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    plot.tag = element_text(face = "bold", size = 17),
    plot.tag.position = c(-0.06, 1.03),
    axis.title = element_text(size = 11),
    axis.text = element_text(color = "black", size = 9.5),
    axis.line = element_line(linewidth = 0.45, color = "black"),
    axis.ticks = element_line(linewidth = 0.4, color = "black"),
    plot.margin = margin(t = 7, r = 7, b = 7, l = 14)
  )

ggsave(
  file.path(PLOT_DIR, "Figure_2c_taxonomic_distance_decay.png"),
  p_ddp,
  width = 7.2,
  height = 4.8,
  units = "in",
  dpi = 600,
  bg = "white"
)
ggsave(
  file.path(PLOT_DIR, "Figure_2c_taxonomic_distance_decay.pdf"),
  p_ddp,
  width = 7.2,
  height = 4.8,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)
ggsave(
  file.path(PLOT_DIR, "Figure_2c_taxonomic_distance_decay.tiff"),
  p_ddp,
  width = 7.2,
  height = 4.8,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

saveRDS(p_ddp, file.path(RDS_DIR, "06_Figure_2c_ggplot_object.rds"))
saveRDS(
  list(
    sample_table = sample_table,
    dist_bray = dist_bray,
    dist_sorensen = dist_sorensen,
    geo_matrix_km = geo_matrix,
    mantel_bray = mantel_bray,
    mantel_sorensen = mantel_sorensen,
    mantel_metric_contrast = mantel_metric_contrast,
    regression_statistics = regression_statistics,
    metric_contrast_statistics = metric_contrast_statistics,
    parameters = list(
      seed = SEED,
      n_permutations = N_PERM,
      n_site_balanced_replicates = N_SITE_BALANCED,
      rarefaction_depth = RAREFACTION_DEPTH,
      min_geo_km = MIN_GEO_KM
    )
  ),
  file.path(RDS_DIR, "06_panel_c_analysis_objects.rds")
)

## -----------------------------------------------------------------------------
## 10. Run summary and provenance
## -----------------------------------------------------------------------------

run_summary <- data.frame(
  status = "PASS",
  run_name = RUN_NAME,
  phyloseq_path = PHYLOSEQ_PATH,
  metadata_path = METADATA_PATH,
  cohort_path = COHORT_PATH,
  phyloseq_samples = nsamples(phy0),
  balanced_cohort_samples = length(cohort_ids),
  analysis_samples_with_coordinates = length(analysis_ids),
  coordinate_defined_sites = n_distinct(sample_table$site_id),
  cohort_taxa_nonzero_before_rarefaction = cohort_taxa_nonzero,
  taxa_nonzero_after_rarefaction = ncol(otu_rarefied),
  rarefaction_depth = RAREFACTION_DEPTH,
  eligible_pairs = nrow(pair_base),
  minimum_geographic_distance_km = MIN_GEO_KM,
  permutations = N_PERM,
  site_balanced_replicates = N_SITE_BALANCED,
  seed = SEED,
  latitude_column = lat_col,
  longitude_column = lon_col,
  study_column = ifelse(is.na(study_col), "NOT_FOUND", study_col),
  stringsAsFactors = FALSE
)

data.table::fwrite(run_summary, file.path(TABLE_DIR, "06_run_summary.csv"))

input_provenance <- data.frame(
  input = c("phyloseq", "metadata", "cohort_pointer", "cohort_table"),
  path = c(PHYLOSEQ_PATH, METADATA_PATH, COHORT_POINTER, COHORT_PATH),
  size_bytes = as.numeric(file.info(c(
    PHYLOSEQ_PATH,
    METADATA_PATH,
    COHORT_POINTER,
    COHORT_PATH
  ))$size),
  modified_time = as.character(file.info(c(
    PHYLOSEQ_PATH,
    METADATA_PATH,
    COHORT_POINTER,
    COHORT_PATH
  ))$mtime),
  stringsAsFactors = FALSE
)
data.table::fwrite(input_provenance, file.path(TABLE_DIR, "06_input_provenance.csv"))

writeLines(
  capture.output(sessionInfo()),
  file.path(LOG_DIR, "06_sessionInfo.txt")
)

log_msg("Bray-Curtis slope: ", sprintf("%.6f", fit_bray[["slope"]]))
log_msg("Sorensen slope: ", sprintf("%.6f", fit_sorensen[["slope"]]))
log_msg(
  "Slope difference (Bray-Curtis minus Sorensen): ",
  sprintf("%.6f", fit_bray[["slope"]] - fit_sorensen[["slope"]])
)
log_msg("Figure: ", file.path(PLOT_DIR, "Figure_2c_taxonomic_distance_decay.png"))
log_msg("=== SCRIPT COMPLETED SUCCESSFULLY ===")
