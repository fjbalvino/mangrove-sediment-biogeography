#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# ==============================================================================
# 05D_CORE_NONCORE_DISTANCE_DECAY.R
# Compare geographic and climatic distance-decay of core and non-core taxa.
#
# Uses core membership and pooled REALM || Locality counts from 05A. Bray-Curtis
# dissimilarity is calculated separately for each fraction after normalizing
# its counts within each locality. Geographic distance is Haversine distance
# between locality centroids; pairs <=0.1 km apart are excluded. Climatic distance
# is Euclidean distance across standardized BIO1, BIO2, BIO3, BIO13 and BIO15 for
# localities with complete observations.
#
# Mantel tests and paired slope comparisons use the same 9,999 locality-label
# permutations. Ordinary linear-model P values for sample pairs are not used
# for inference because the pairs are not independent observations.
# ==============================================================================

options(stringsAsFactors = FALSE, scipen = 999, bitmapType = "cairo")
Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
set.seed(120319)

suppressPackageStartupMessages({
  library(vegan)
  library(geosphere)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

# ============================== 1) Utilidades ================================

msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      ..., "\n", sep = "")
}

stopf <- function(...) stop(sprintf(...), call. = FALSE)

resolve_col <- function(df, candidates, label, required = TRUE) {
  nm <- names(df)
  hit <- candidates[candidates %in% nm]
  if (length(hit) == 0L) {
    idx <- match(tolower(candidates), tolower(nm), nomatch = 0L)
    hit <- nm[idx[idx > 0L]]
  }
  hit <- unique(hit)
  if (length(hit) == 0L) {
    if (required) {
      stopf("No se encontro la columna %s. Candidatas: %s",
            label, paste(candidates, collapse = ", "))
    }
    return(NA_character_)
  }
  hit[[1L]]
}

finite_median <- function(x) {
  x <- suppressWarnings(as.numeric(as.character(x)))
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else median(x)
}

n_distinct_finite <- function(x, digits = 10L) {
  x <- suppressWarnings(as.numeric(as.character(x)))
  x <- x[is.finite(x)]
  length(unique(round(x, digits = digits)))
}

fast_slope <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 3L || stats::var(x) <= 0) return(NA_real_)
  sum((x - mean(x)) * (y - mean(y))) /
    sum((x - mean(x))^2)
}

perm_p_two_sided <- function(observed, permuted) {
  permuted <- permuted[is.finite(permuted)]
  if (!is.finite(observed) || length(permuted) == 0L) return(NA_real_)
  (1 + sum(abs(permuted) >= abs(observed))) / (length(permuted) + 1)
}

fmt_p <- function(p) {
  if (!is.finite(p)) return("P = NA")
  if (p < 0.001) "P < 0.001" else sprintf("P = %.3f", p)
}

run_decay_test <- function(core_dissimilarity,
                           noncore_dissimilarity,
                           predictor_distance,
                           slope_x,
                           pair_keep,
                           n_permutations,
                           seed,
                           analysis_name) {
  core_mat <- as.matrix(core_dissimilarity)
  non_mat <- as.matrix(noncore_dissimilarity)
  pred_mat <- as.matrix(predictor_distance)
  x_mat <- as.matrix(slope_x)

  expected_names <- rownames(pred_mat)
  matrices <- list(core = core_mat, noncore = non_mat,
                   predictor = pred_mat, slope_predictor = x_mat)
  if (is.null(expected_names) || anyNA(expected_names) ||
      anyDuplicated(expected_names) || any(!nzchar(expected_names))) {
    stopf("%s requires unique, nonmissing locality labels.", analysis_name)
  }
  for (matrix_name in names(matrices)) {
    current <- matrices[[matrix_name]]
    if (!identical(dim(current), dim(pred_mat)) ||
        nrow(current) != ncol(current) ||
        !identical(rownames(current), expected_names) ||
        !identical(colnames(current), expected_names)) {
      stopf("The %s matrix is not aligned in %s.", matrix_name, analysis_name)
    }
  }
  if (!identical(dim(pair_keep), dim(pred_mat)) || anyNA(pair_keep)) {
    stopf("The pair mask is invalid in %s.", analysis_name)
  }
  if (any(!is.finite(pred_mat)) || any(pred_mat < 0)) {
    stopf("Predictor distances must be finite and nonnegative in %s.", analysis_name)
  }
  for (matrix_name in c("core", "noncore")) {
    current <- matrices[[matrix_name]]
    if (any(!is.finite(current)) || any(current < -1e-12 | current > 1 + 1e-12) ||
        max(abs(current - t(current))) > 1e-10 ||
        any(abs(diag(current)) > 1e-10)) {
      stopf("Invalid Bray-Curtis dissimilarities in %s (%s).", analysis_name, matrix_name)
    }
  }

  mask <- upper.tri(pred_mat) & pair_keep &
    is.finite(pred_mat) & is.finite(x_mat)

  pair_idx <- which(mask, arr.ind = TRUE)
  if (nrow(pair_idx) < 10L) {
    stopf("%s conserva menos de 10 pares validos.", analysis_name)
  }

  pred_vec <- pred_mat[mask]
  x_vec <- x_mat[mask]
  core_d_vec <- core_mat[mask]
  non_d_vec <- non_mat[mask]
  core_sim_vec <- 1 - core_d_vec
  non_sim_vec <- 1 - non_d_vec

  observed <- tibble(
    fraction = c("Core-associated", "Non-core"),
    mantel_r = c(
      cor(core_d_vec, pred_vec, method = "pearson"),
      cor(non_d_vec, pred_vec, method = "pearson")
    ),
    similarity_slope = c(
      fast_slope(x_vec, core_sim_vec),
      fast_slope(x_vec, non_sim_vec)
    )
  )

  observed_delta <- observed$similarity_slope[observed$fraction ==
                                                "Core-associated"] -
    observed$similarity_slope[observed$fraction == "Non-core"]

  set.seed(seed)
  n_units <- nrow(pred_mat)
  perm_core_r <- numeric(n_permutations)
  perm_non_r <- numeric(n_permutations)
  perm_delta_slope <- numeric(n_permutations)

  for (b in seq_len(n_permutations)) {
    idx <- sample.int(n_units, size = n_units, replace = FALSE)
    core_perm <- core_mat[idx, idx, drop = FALSE][mask]
    non_perm <- non_mat[idx, idx, drop = FALSE][mask]

    perm_core_r[b] <- cor(core_perm, pred_vec, method = "pearson")
    perm_non_r[b] <- cor(non_perm, pred_vec, method = "pearson")
    perm_delta_slope[b] <-
      fast_slope(x_vec, 1 - core_perm) -
      fast_slope(x_vec, 1 - non_perm)
  }

  observed <- observed %>%
    mutate(
      permutation_p = c(
        perm_p_two_sided(mantel_r[1], perm_core_r),
        perm_p_two_sided(mantel_r[2], perm_non_r)
      ),
      analysis = analysis_name,
      n_units = n_units,
      n_pairs = length(pred_vec),
      .before = 1
    )

  contrast <- tibble(
    analysis = analysis_name,
    contrast = "Core-associated minus Non-core similarity slope",
    delta_slope = observed_delta,
    permutation_p = perm_p_two_sided(observed_delta, perm_delta_slope),
    n_permutations = n_permutations,
    n_units = n_units,
    n_pairs = length(pred_vec)
  )

  permutation_table <- tibble(
    analysis = analysis_name,
    permutation = seq_len(n_permutations),
    core_mantel_r = perm_core_r,
    noncore_mantel_r = perm_non_r,
    delta_similarity_slope = perm_delta_slope
  )

  pair_table <- tibble(
    locality_1 = rownames(pred_mat)[pair_idx[, 1]],
    locality_2 = colnames(pred_mat)[pair_idx[, 2]],
    predictor_distance = pred_vec,
    slope_predictor = x_vec,
    core_dissimilarity = core_d_vec,
    noncore_dissimilarity = non_d_vec,
    core_similarity = core_sim_vec,
    noncore_similarity = non_sim_vec
  )

  list(
    statistics = observed,
    contrast = contrast,
    permutations = permutation_table,
    pairs = pair_table
  )
}

# ============================== 2) Rutas =====================================

PROJECT_DIR <- Sys.getenv("BIOGEOGRAPHY_BASE", unset = ".")
BASE_DIR <- PROJECT_DIR
RESULTS_DIR <- Sys.getenv("BIOGEOGRAPHY_RESULTS", unset = file.path(PROJECT_DIR, "results"))

METADATA_CSV <- Sys.getenv(
  "MAGENTA_METADATA_CSV",
  unset = file.path(PROJECT_DIR, "data", "environmental_metadata.csv")
)

SCRIPT_09A_TAG <- "05A_TAXONOMIC_CORE_ASSOCIATED"
LATEST_09A_PTR <- file.path(
  RESULTS_DIR, paste0("LATEST_", SCRIPT_09A_TAG, ".txt")
)

SCRIPT_TAG <- "05D_CORE_NONCORE_DISTANCE_DECAY"
RUN_STAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUT_DIR <- file.path(RESULTS_DIR, paste0(SCRIPT_TAG, "_", RUN_STAMP))
TABLE_DIR <- file.path(OUT_DIR, "tables")
PLOT_DIR <- file.path(OUT_DIR, "plots")
RDS_DIR <- file.path(OUT_DIR, "rds")

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)

LATEST_PTR <- file.path(RESULTS_DIR, paste0("LATEST_", SCRIPT_TAG, ".txt"))
STATUS_FILE <- file.path(OUT_DIR, "STATUS.txt")
writeLines("RUNNING", STATUS_FILE)

completed_successfully <- FALSE
on.exit({
  writeLines(
    if (completed_successfully) "PASS" else "FAILED",
    STATUS_FILE
  )
}, add = TRUE)

# ============================ 3) Parametros ==================================

N_PERMUTATIONS <- 9999L
MIN_GEOGRAPHIC_DISTANCE_KM <- 0.1
SEED_GEOGRAPHIC <- 120319L
SEED_ENVIRONMENTAL <- 120320L

CLIMATE_SPECS <- list(
  BIO1 = c("bio1", "BIO1", "bio1_mean", "BIO1_mean"),
  BIO2 = c("bio2", "BIO2", "bio2_mean", "BIO2_mean"),
  BIO3 = c("bio3", "BIO3", "bio3_mean", "BIO3_mean"),
  BIO13 = c("bio13", "BIO13", "bio13_mean", "BIO13_mean"),
  BIO15 = c("bio15", "BIO15", "bio15_mean", "BIO15_mean")
)

CORE_COLOUR <- "#B2182B"
NONCORE_COLOUR <- "#737373"
FRACTION_COLOURS <- c(
  "Core-associated" = CORE_COLOUR,
  "Non-core" = NONCORE_COLOUR
)

msg("Starting ", SCRIPT_TAG)
msg("Metadata: ", METADATA_CSV)
msg("Output: ", OUT_DIR)

if (!file.exists(METADATA_CSV)) stopf("No se encontro la metadata: %s", METADATA_CSV)
if (!file.exists(LATEST_09A_PTR)) {
  stopf("No se encontro el puntero 09A: %s", LATEST_09A_PTR)
}

RUN_09A <- trimws(readLines(LATEST_09A_PTR, warn = FALSE)[1])
if (!nzchar(RUN_09A) || !dir.exists(RUN_09A)) {
  stopf("El run 09A indicado por LATEST no existe: %s", RUN_09A)
}

INPUT_09A_RDS <- file.path(
  RUN_09A, "rds", "09A_Figure_3a_analysis_objects.rds"
)
if (!file.exists(INPUT_09A_RDS)) {
  stopf("No se encontro el objeto 09A: %s", INPUT_09A_RDS)
}

msg("Input 09A run: ", RUN_09A)
msg("Input 09A object: ", INPUT_09A_RDS)

# =========================== 4) Objeto 09A ===================================

obj09a <- readRDS(INPUT_09A_RDS)

required_09a <- c(
  "metadata_aligned",
  "locality_metadata",
  "locality_counts",
  "taxon_classification",
  "taxon_id_harmonization_audit",
  "parameters"
)
missing_09a <- setdiff(required_09a, names(obj09a))
if (length(missing_09a) > 0L) {
  stopf(
    "El objeto 09A no contiene: %s",
    paste(missing_09a, collapse = ", ")
  )
}

meta_09a <- as.data.frame(obj09a$metadata_aligned)
locality_meta_09a <- as.data.frame(obj09a$locality_metadata)
otu_locality <- as.matrix(obj09a$locality_counts)
taxon_class <- as.data.frame(obj09a$taxon_classification)
harmonization_audit <- as.data.frame(obj09a$taxon_id_harmonization_audit)
parameters_09a <- as.data.frame(obj09a$parameters)

if (!all(c("parameter", "value") %in% names(parameters_09a))) {
  stopf("El objeto 09A no contiene una tabla de parametros valida.")
}

parameter_value <- function(key) {
  value <- parameters_09a$value[parameters_09a$parameter == key]
  if (length(value) != 1L) {
    stopf("El parametro 09A '%s' esta ausente o duplicado.", key)
  }
  value[[1L]]
}

expected_09a <- c(
  minimum_global_reads = 50,
  minimum_global_localities = 5,
  minimum_realm_localities_to_support_core = 3,
  core_realm_prevalence = 1,
  core_mean_relative_abundance = 0.0001
)

observed_09a <- vapply(
  names(expected_09a),
  function(key) suppressWarnings(as.numeric(parameter_value(key))),
  numeric(1)
)

if (any(!is.finite(observed_09a)) ||
    any(abs(observed_09a - expected_09a) > 1e-12)) {
  stopf(
    paste0(
      "El ultimo 09A no usa la definicion final del core. Esperado: ",
      ">=50 lecturas, >=5 localidades, realms con >=3 localidades, ",
      "prevalencia=1 y abundancia media>=0.0001."
    )
  )
}

required_class <- c("taxon_id", "core_associated", "component")
if (!all(required_class %in% names(taxon_class))) {
  stopf(
    "Faltan columnas en taxon_classification: %s",
    paste(setdiff(required_class, names(taxon_class)), collapse = ", ")
  )
}
if (anyDuplicated(taxon_class$taxon_id)) {
  stopf("taxon_classification contiene taxon_id duplicados.")
}
if (is.null(rownames(otu_locality)) || is.null(colnames(otu_locality))) {
  stopf("locality_counts debe tener rownames y colnames.")
}
if (any(!is.finite(otu_locality)) || any(otu_locality < 0)) {
  stopf("locality_counts contiene valores no finitos o negativos.")
}
if (anyDuplicated(rownames(otu_locality)) || anyDuplicated(colnames(otu_locality))) {
  stopf("locality_counts contiene IDs duplicados.")
}

required_locality_meta <- c("locality_id", "Locality", "REALM")
if (!all(required_locality_meta %in% names(locality_meta_09a))) {
  stopf(
    "Faltan columnas en locality_metadata: %s",
    paste(setdiff(required_locality_meta, names(locality_meta_09a)),
          collapse = ", ")
  )
}

locality_meta_09a <- locality_meta_09a[
  match(rownames(otu_locality), locality_meta_09a$locality_id),
  , drop = FALSE
]
if (anyNA(locality_meta_09a$locality_id)) {
  stopf("Fallo la alineacion de locality_counts con locality_metadata.")
}

eligible_ids <- as.character(taxon_class$taxon_id)
core_ids <- as.character(
  taxon_class$taxon_id[taxon_class$core_associated %in% TRUE]
)
noncore_ids <- as.character(
  taxon_class$taxon_id[!(taxon_class$core_associated %in% TRUE)]
)

if (length(core_ids) == 0L || length(noncore_ids) == 0L) {
  stopf("La clasificacion 09A no contiene ambas fracciones.")
}
missing_eligible <- setdiff(eligible_ids, colnames(otu_locality))
if (length(missing_eligible) > 0L) {
  stopf(
    "Hay %d taxa elegibles ausentes de locality_counts. Primeros: %s",
    length(missing_eligible), paste(head(missing_eligible, 20), collapse = ", ")
  )
}

core_counts_all <- otu_locality[, core_ids, drop = FALSE]
non_counts_all <- otu_locality[, noncore_ids, drop = FALSE]

if (any(rowSums(core_counts_all) <= 0)) {
  stopf("Hay localidades sin lecturas core-associated.")
}
if (any(rowSums(non_counts_all) <= 0)) {
  stopf("Hay localidades sin lecturas non-core.")
}

msg("09A samples: ", nrow(meta_09a))
msg("09A realm-locality units: ", nrow(otu_locality))
msg("09A eligible taxa: ", length(eligible_ids))
msg("Core-associated taxa: ", length(core_ids))
msg("Non-core taxa: ", length(noncore_ids))

# ====================== 5) Numeric environmental metadata ==========================

meta_raw <- read_csv(METADATA_CSV, show_col_types = FALSE, progress = FALSE)

meta_key_col <- resolve_col(meta_raw, c("run_accession"), "run_accession")
lat_col <- resolve_col(
  meta_raw,
  c("lat", "latitude", "Latitude", "decimalLatitude", "decimal_latitude"),
  "latitude"
)
lon_col <- resolve_col(
  meta_raw,
  c("lon", "longitude", "Longitude", "decimalLongitude", "decimal_longitude"),
  "longitude"
)

climate_source_cols <- vapply(
  names(CLIMATE_SPECS),
  function(nm) resolve_col(meta_raw, CLIMATE_SPECS[[nm]], nm),
  character(1)
)

meta_numeric <- meta_raw %>%
  transmute(
    run_accession = trimws(as.character(.data[[meta_key_col]])),
    latitude = suppressWarnings(as.numeric(as.character(.data[[lat_col]]))),
    longitude = suppressWarnings(as.numeric(as.character(.data[[lon_col]]))),
    BIO1 = suppressWarnings(as.numeric(as.character(
      .data[[climate_source_cols[["BIO1"]]]]
    ))),
    BIO2 = suppressWarnings(as.numeric(as.character(
      .data[[climate_source_cols[["BIO2"]]]]
    ))),
    BIO3 = suppressWarnings(as.numeric(as.character(
      .data[[climate_source_cols[["BIO3"]]]]
    ))),
    BIO13 = suppressWarnings(as.numeric(as.character(
      .data[[climate_source_cols[["BIO13"]]]]
    ))),
    BIO15 = suppressWarnings(as.numeric(as.character(
      .data[[climate_source_cols[["BIO15"]]]]
    )))
  ) %>%
  filter(!is.na(run_accession), nzchar(run_accession))

numeric_vars <- c("latitude", "longitude", names(CLIMATE_SPECS))

duplicate_audit <- meta_numeric %>%
  group_by(run_accession) %>%
  summarise(
    n_rows = n(),
    across(
      all_of(numeric_vars),
      list(n_distinct_finite = ~ n_distinct_finite(.x)),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  mutate(
    any_numeric_discordance = if_any(
      ends_with("_n_distinct_finite"), ~ .x > 1L
    )
  )

write_tsv(
  duplicate_audit,
  file.path(TABLE_DIR, "12_metadata_duplicate_run_accession_audit.tsv")
)

n_discordant_runs <- sum(duplicate_audit$any_numeric_discordance)
if (n_discordant_runs > 0L) {
  msg(
    "WARNING: ", n_discordant_runs,
    " run_accession have discordant numeric metadata; medians will be used."
  )
}

invalid_coordinates <- meta_numeric %>%
  filter((is.finite(latitude) & (latitude < -90 | latitude > 90)) |
         (is.finite(longitude) & (longitude < -180 | longitude > 180)))
if (nrow(invalid_coordinates) > 0L) {
  write_tsv(invalid_coordinates,
            file.path(TABLE_DIR, "12_invalid_geographic_coordinates_STOP.tsv"))
  stopf("Metadata contain %d coordinate records outside geographic bounds; see the audit table.",
        nrow(invalid_coordinates))
}

meta_one <- meta_numeric %>%
  group_by(run_accession) %>%
  summarise(across(all_of(numeric_vars), finite_median), .groups = "drop")

required_meta09 <- c("run_accession", "Locality", "REALM")
if (!all(required_meta09 %in% names(meta_09a))) {
  stopf(
    "Faltan columnas en metadata_aligned de 09A: %s",
    paste(setdiff(required_meta09, names(meta_09a)), collapse = ", ")
  )
}

sample_meta <- as_tibble(meta_09a) %>%
  select(run_accession, Locality, REALM) %>%
  mutate(
    run_accession = trimws(as.character(run_accession)),
    Locality = trimws(as.character(Locality)),
    REALM = trimws(as.character(REALM)),
    locality_id = paste(REALM, Locality, sep = " || ")
  ) %>%
  left_join(meta_one, by = "run_accession")

missing_numeric_rows <- sample_meta %>%
  filter(if_all(all_of(numeric_vars), is.na))

write_tsv(
  missing_numeric_rows,
  file.path(
    TABLE_DIR,
    "12_samples_without_coordinates_or_climate_metadata.tsv"
  )
)

if (nrow(missing_numeric_rows) > 0L) {
  msg(
    "WARNING: coordinates and climate are completely missing for ",
    nrow(missing_numeric_rows),
    paste0(
      " run_accession. They will not stop the analysis; availability ",
      "will be evaluated after aggregation by REALM || Locality."
    )
  )
}

locality_numeric <- sample_meta %>%
  group_by(locality_id, Locality, REALM) %>%
  summarise(
    n_metagenomes = n(),
    across(all_of(numeric_vars), finite_median),
    .groups = "drop"
  )

locality_numeric <- locality_numeric[
  match(rownames(otu_locality), locality_numeric$locality_id),
  , drop = FALSE
]
if (anyNA(locality_numeric$locality_id)) {
  stopf("Fallo la alineacion de metadata numerica con locality_counts.")
}

locality_numeric <- locality_numeric %>%
  mutate(
    has_coordinates = is.finite(longitude) & is.finite(latitude) &
      longitude >= -180 & longitude <= 180 & latitude >= -90 & latitude <= 90,
    n_complete_climate_predictors = rowSums(
      is.finite(as.matrix(across(all_of(names(CLIMATE_SPECS)))))
    ),
    has_complete_climate =
      n_complete_climate_predictors == length(CLIMATE_SPECS),
    geographic_exclusion_reason = if_else(
      has_coordinates,
      NA_character_,
      "No finite locality coordinates"
    ),
    environmental_exclusion_reason = if_else(
      has_complete_climate,
      NA_character_,
      paste0(
        "Complete climate predictors: ",
        n_complete_climate_predictors,
        "/",
        length(CLIMATE_SPECS)
      )
    )
  )

locality_availability_audit <- locality_numeric %>%
  select(
    locality_id, Locality, REALM, n_metagenomes,
    has_coordinates, geographic_exclusion_reason,
    n_complete_climate_predictors, has_complete_climate,
    environmental_exclusion_reason
  )

write_tsv(
  locality_availability_audit,
  file.path(TABLE_DIR, "12_locality_analysis_availability_audit.tsv")
)

n_localities_without_coordinates <- sum(!locality_numeric$has_coordinates)
n_localities_without_complete_climate <- sum(
  !locality_numeric$has_complete_climate
)

if (n_localities_without_coordinates > 0L) {
  msg(
    "WARNING: ", n_localities_without_coordinates,
    " locality units lack finite coordinates and will be excluded only from panel 3c."
  )
}
if (n_localities_without_complete_climate > 0L) {
  msg(
    "WARNING: ", n_localities_without_complete_climate,
    " locality units lack one or more climate predictors and will be excluded only from panel 3d."
  )
}

coord_spread <- sample_meta %>%
  filter(is.finite(longitude), is.finite(latitude)) %>%
  group_by(locality_id, Locality, REALM) %>%
  group_modify(function(.x, .y) {
    centroid_lon <- finite_median(.x$longitude)
    centroid_lat <- finite_median(.x$latitude)
    distances <- geosphere::distHaversine(
      cbind(.x$longitude, .x$latitude),
      c(centroid_lon, centroid_lat)
    ) / 1000
    tibble(
      n_coordinate_records = nrow(.x),
      n_unique_coordinate_pairs = nrow(unique(.x[c("longitude", "latitude")])),
      centroid_longitude = centroid_lon,
      centroid_latitude = centroid_lat,
      maximum_distance_to_centroid_km = max(distances, na.rm = TRUE)
    )
  }) %>%
  ungroup()

write_tsv(
  locality_numeric,
  file.path(TABLE_DIR, "12_locality_metadata_for_distance_decay.tsv")
)
write_tsv(
  coord_spread,
  file.path(TABLE_DIR, "12_locality_coordinate_spread_audit.tsv")
)

# =================== 6) Distancias de comunidad =============================

compute_fraction_dissimilarity <- function(count_matrix) {
  totals <- rowSums(count_matrix)
  if (any(!is.finite(totals)) || any(totals <= 0)) {
    stopf("Una fraccion contiene localidades con cero lecturas.")
  }
  relative <- sweep(count_matrix, 1, totals, "/")
  vegan::vegdist(relative, method = "bray")
}

# -------------------------- Panel 3c: geografia ------------------------------

geo_keep <- locality_numeric$has_coordinates

geo_meta <- locality_numeric[geo_keep, , drop = FALSE]
if (nrow(geo_meta) < 4L) stopf("Menos de cuatro localidades con coordenadas.")

geo_ids <- geo_meta$locality_id
core_geo_counts <- core_counts_all[geo_ids, , drop = FALSE]
non_geo_counts <- non_counts_all[geo_ids, , drop = FALSE]

core_geo_dist <- compute_fraction_dissimilarity(core_geo_counts)
non_geo_dist <- compute_fraction_dissimilarity(non_geo_counts)

geo_matrix <- geosphere::distm(
  cbind(geo_meta$longitude, geo_meta$latitude),
  fun = geosphere::distHaversine
) / 1000
rownames(geo_matrix) <- geo_ids
colnames(geo_matrix) <- geo_ids

log_geo_matrix <- log10(geo_matrix)
diag(log_geo_matrix) <- NA_real_
geo_pair_keep <- geo_matrix > MIN_GEOGRAPHIC_DISTANCE_KM
diag(geo_pair_keep) <- FALSE

geo_result <- run_decay_test(
  core_dissimilarity = core_geo_dist,
  noncore_dissimilarity = non_geo_dist,
  predictor_distance = geo_matrix,
  slope_x = log_geo_matrix,
  pair_keep = geo_pair_keep,
  n_permutations = N_PERMUTATIONS,
  seed = SEED_GEOGRAPHIC,
  analysis_name = "Geographic"
)

geo_pairs <- geo_result$pairs %>%
  rename(
    geographic_distance_km = predictor_distance,
    log10_geographic_distance = slope_predictor
  )

# ------------------------- Panel 3d: ambiente --------------------------------

climate_vars <- names(CLIMATE_SPECS)
env_keep <- locality_numeric$has_complete_climate
env_meta <- locality_numeric[env_keep, , drop = FALSE]
if (nrow(env_meta) < 4L) {
  stopf("Menos de cuatro localidades con los cinco predictores climaticos.")
}

env_raw <- as.matrix(env_meta[, climate_vars, drop = FALSE])
storage.mode(env_raw) <- "double"

env_sd <- apply(env_raw, 2, sd)
if (any(!is.finite(env_sd) | env_sd <= 0)) {
  bad <- names(env_sd)[!is.finite(env_sd) | env_sd <= 0]
  stopf("Predictores climaticos sin varianza: %s", paste(bad, collapse = ", "))
}

env_scaled <- scale(env_raw)
rownames(env_scaled) <- env_meta$locality_id
env_matrix <- as.matrix(dist(env_scaled, method = "euclidean"))

env_ids <- env_meta$locality_id
core_env_counts <- core_counts_all[env_ids, , drop = FALSE]
non_env_counts <- non_counts_all[env_ids, , drop = FALSE]

core_env_dist <- compute_fraction_dissimilarity(core_env_counts)
non_env_dist <- compute_fraction_dissimilarity(non_env_counts)

env_pair_keep <- matrix(
  TRUE, nrow = nrow(env_matrix), ncol = ncol(env_matrix),
  dimnames = dimnames(env_matrix)
)
diag(env_pair_keep) <- FALSE

env_result <- run_decay_test(
  core_dissimilarity = core_env_dist,
  noncore_dissimilarity = non_env_dist,
  predictor_distance = env_matrix,
  slope_x = env_matrix,
  pair_keep = env_pair_keep,
  n_permutations = N_PERMUTATIONS,
  seed = SEED_ENVIRONMENTAL,
  analysis_name = "Environmental"
)

env_pairs <- env_result$pairs %>%
  rename(
    environmental_distance = predictor_distance,
    environmental_distance_for_slope = slope_predictor
  )

env_scaling <- tibble(
  predictor = climate_vars,
  source_column = unname(climate_source_cols[climate_vars]),
  locality_mean = as.numeric(attr(env_scaled, "scaled:center")),
  locality_sd = as.numeric(attr(env_scaled, "scaled:scale")),
  n_complete_localities = nrow(env_scaled)
)

# ============================== 7) Tablas ====================================

statistics_table <- bind_rows(
  geo_result$statistics,
  env_result$statistics
)

contrast_table <- bind_rows(
  geo_result$contrast,
  env_result$contrast
)

permutation_table <- bind_rows(
  geo_result$permutations,
  env_result$permutations
)

write_tsv(
  statistics_table,
  file.path(TABLE_DIR, "12_distance_decay_statistics.tsv")
)
write_tsv(
  contrast_table,
  file.path(TABLE_DIR, "12_core_noncore_slope_contrasts.tsv")
)
write_tsv(
  geo_pairs,
  file.path(TABLE_DIR, "12_geographic_pairwise_data.tsv")
)
write_tsv(
  env_pairs,
  file.path(TABLE_DIR, "12_environmental_pairwise_data.tsv")
)
write_tsv(
  env_scaling,
  file.path(TABLE_DIR, "12_environmental_predictor_scaling.tsv")
)
write_tsv(
  permutation_table,
  file.path(TABLE_DIR, "12_permutation_distributions.tsv.gz")
)

# ============================== 8) Graficos ==================================

geo_plot_data <- geo_pairs %>%
  transmute(
    locality_1,
    locality_2,
    geographic_distance_km,
    core_similarity,
    noncore_similarity
  ) %>%
  pivot_longer(
    cols = c(core_similarity, noncore_similarity),
    names_to = "fraction_code",
    values_to = "similarity"
  ) %>%
  mutate(
    fraction = recode(
      fraction_code,
      core_similarity = "Core-associated",
      noncore_similarity = "Non-core"
    ),
    fraction = factor(
      fraction,
      levels = c("Core-associated", "Non-core")
    )
  )

env_plot_data <- env_pairs %>%
  transmute(
    locality_1,
    locality_2,
    environmental_distance,
    core_similarity,
    noncore_similarity
  ) %>%
  pivot_longer(
    cols = c(core_similarity, noncore_similarity),
    names_to = "fraction_code",
    values_to = "similarity"
  ) %>%
  mutate(
    fraction = recode(
      fraction_code,
      core_similarity = "Core-associated",
      noncore_similarity = "Non-core"
    ),
    fraction = factor(
      fraction,
      levels = c("Core-associated", "Non-core")
    )
  )

make_stat_labels <- function(statistics, contrast) {
  core_row <- statistics %>% filter(fraction == "Core-associated")
  non_row <- statistics %>% filter(fraction == "Non-core")
  list(
    core = sprintf(
      "Core-associated: r = %.3f, slope = %.3f, %s",
      core_row$mantel_r, core_row$similarity_slope,
      fmt_p(core_row$permutation_p)
    ),
    non = sprintf(
      "Non-core: r = %.3f, slope = %.3f, %s",
      non_row$mantel_r, non_row$similarity_slope,
      fmt_p(non_row$permutation_p)
    ),
    contrast = sprintf(
      "Slope difference: Δ = %.3f, %s",
      contrast$delta_slope, fmt_p(contrast$permutation_p)
    )
  )
}

geo_labels <- make_stat_labels(
  geo_result$statistics,
  geo_result$contrast
)
env_labels <- make_stat_labels(
  env_result$statistics,
  env_result$contrast
)

base_decay_plot <- function(data, x_var, x_label, labels) {
  ggplot(
    data,
    aes(x = .data[[x_var]], y = similarity, colour = fraction)
  ) +
    geom_point(size = 1.25, alpha = 0.32, stroke = 0) +
    geom_smooth(
      method = "lm", formula = y ~ x, se = FALSE,
      linewidth = 1.15
    ) +
    annotate(
      "text", x = Inf, y = 0.245, label = labels$core,
      hjust = 1.03, vjust = 0.5,
      colour = CORE_COLOUR, size = 4.60
    ) +
    annotate(
      "text", x = Inf, y = 0.155, label = labels$non,
      hjust = 1.03, vjust = 0.5,
      colour = "#595959", size = 4.60
    ) +
    annotate(
      "text", x = Inf, y = 0.065, label = labels$contrast,
      hjust = 1.03, vjust = 0.5,
      colour = "black", size = 4.45
    ) +
    scale_colour_manual(
      values = FRACTION_COLOURS,
      breaks = c("Core-associated", "Non-core"),
      name = NULL,
      drop = FALSE
    ) +
    scale_y_continuous(
      name = "Community similarity (1 - dissimilarity)",
      limits = c(0, 1), breaks = seq(0, 1, 0.25),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(x = x_label) +
    guides(colour = guide_legend(nrow = 1, byrow = TRUE)) +
    theme_classic(base_size = 12) +
    theme(
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      legend.position = "bottom",
      legend.justification = "center",
      legend.text = element_text(size = 10),
      legend.key.width = grid::unit(0.65, "cm"),
      plot.margin = margin(8, 12, 5, 8)
    )
}

p_geo <- base_decay_plot(
  geo_plot_data,
  "geographic_distance_km",
  "Geographic distance (km, log10 scale)",
  geo_labels
) +
  scale_x_log10(
    breaks = c(0.1, 1, 10, 100, 1000, 10000),
    labels = c("0.1", "1", "10", "100", "1,000", "10,000")
  )

p_env <- base_decay_plot(
  env_plot_data,
  "environmental_distance",
  "Environmental distance",
  env_labels
)

p_env_combined <- p_env + labs(y = NULL)

p_combined <- (p_geo | p_env_combined) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

save_plot_triplet <- function(plot_object, stem, width, height) {
  png_file <- file.path(PLOT_DIR, paste0(stem, ".png"))
  pdf_file <- file.path(PLOT_DIR, paste0(stem, ".pdf"))
  tiff_file <- file.path(PLOT_DIR, paste0(stem, ".tiff"))

  ggsave(
    png_file, plot_object, width = width, height = height,
    units = "in", dpi = 600, bg = "white"
  )
  ggsave(
    pdf_file, plot_object, width = width, height = height,
    units = "in", device = cairo_pdf, bg = "white"
  )
  ggsave(
    tiff_file, plot_object, width = width, height = height,
    units = "in", dpi = 600, compression = "lzw", bg = "white"
  )
  png_file
}

png_geo <- save_plot_triplet(
  p_geo, "Figure_3c_geographic_core_noncore_decay", 7.2, 5.2
)
png_env <- save_plot_triplet(
  p_env, "Figure_3d_environmental_core_noncore_decay", 7.2, 5.2
)
png_combined <- save_plot_triplet(
  p_combined, "Figure_3cd_core_noncore_distance_decay_combined", 12.4, 5.2
)

# ============================ 9) Provenance ==================================

parameters <- tibble(
  parameter = c(
    "input_core_classification",
    "analysis_unit",
    "community_dissimilarity",
    "fraction_normalization",
    "geographic_distance",
    "minimum_geographic_distance_km",
    "geographic_slope_predictor",
    "environmental_distance",
    "environmental_predictors",
    "metadata_duplicate_resolution",
    "mantel_correlation",
    "n_permutations",
    "slope_contrast",
    "sample_key"
  ),
  value = c(
    "09A harmonized core-associated classification",
    "REALM || Locality; replicate counts summed upstream by 09A",
    "Bray-Curtis",
    "Relative abundance calculated separately within each fraction and locality",
    "Haversine distance between locality median coordinates",
    MIN_GEOGRAPHIC_DISTANCE_KM,
    "log10 geographic distance in km",
    "Euclidean distance between standardized locality climate predictors",
    paste(climate_vars, collapse = " | "),
    "Median by run_accession, then median by REALM || Locality; audit retained",
    "Pearson; two-sided locality-label permutations",
    N_PERMUTATIONS,
    "Paired core-associated minus non-core similarity-slope difference under the same locality-label permutations",
    "run_accession"
  )
)

git_commit <- tryCatch(
  trimws(system2(
    "git", c("-C", BASE_DIR, "rev-parse", "HEAD"),
    stdout = TRUE, stderr = FALSE
  )[1]),
  error = function(e) NA_character_
)

run_summary <- tibble(
  item = c(
    "script",
    "input_09A_run",
    "input_09A_rds",
    "metadata",
    "n_samples_in_09A",
    "n_realm_locality_units_in_09A",
    "n_eligible_taxa",
    "n_core_associated_taxa",
    "n_noncore_taxa",
    "n_geographic_localities",
    "n_geographic_pairs",
    "geographic_distance_min_km",
    "geographic_distance_max_km",
    "n_environmental_localities",
    "n_environmental_pairs",
    "n_samples_without_coordinates_or_climate",
    "n_localities_without_coordinates",
    "n_localities_without_complete_climate",
    "n_discordant_duplicate_run_accessions",
    "n_permutations",
    "git_commit",
    "output_directory"
  ),
  value = c(
    SCRIPT_TAG,
    RUN_09A,
    INPUT_09A_RDS,
    METADATA_CSV,
    nrow(meta_09a),
    nrow(otu_locality),
    length(eligible_ids),
    length(core_ids),
    length(noncore_ids),
    nrow(geo_meta),
    nrow(geo_pairs),
    min(geo_pairs$geographic_distance_km),
    max(geo_pairs$geographic_distance_km),
    nrow(env_meta),
    nrow(env_pairs),
    nrow(missing_numeric_rows),
    n_localities_without_coordinates,
    n_localities_without_complete_climate,
    n_discordant_runs,
    N_PERMUTATIONS,
    git_commit,
    OUT_DIR
  )
)

write_tsv(parameters, file.path(TABLE_DIR, "12_parameters.tsv"))
write_tsv(run_summary, file.path(TABLE_DIR, "12_run_summary.tsv"))

saveRDS(
  list(
    input_09A_run = RUN_09A,
    parameters = parameters,
    run_summary = run_summary,
    locality_metadata = locality_numeric,
    environmental_scaling = env_scaling,
    core_taxa = core_ids,
    noncore_taxa = noncore_ids,
    geographic = list(
      core_dissimilarity = core_geo_dist,
      noncore_dissimilarity = non_geo_dist,
      predictor_distance = geo_matrix,
      result = geo_result,
      plot = p_geo
    ),
    environmental = list(
      core_dissimilarity = core_env_dist,
      noncore_dissimilarity = non_env_dist,
      predictor_distance = env_matrix,
      result = env_result,
      plot = p_env
    ),
    combined_plot = p_combined
  ),
  file.path(RDS_DIR, "12_Figure_3cd_analysis_objects.rds")
)

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
writeLines(OUT_DIR, LATEST_PTR)

completed_successfully <- TRUE
writeLines("PASS", STATUS_FILE)

msg("Geographic figure: ", png_geo)
msg("Environmental figure: ", png_env)
msg("Combined figure: ", png_combined)
msg("=== SCRIPT COMPLETED SUCCESSFULLY ===")
