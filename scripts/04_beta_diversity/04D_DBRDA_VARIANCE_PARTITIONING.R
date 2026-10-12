#!/usr/bin/env Rscript
# -*- coding: UTF-8 -*-

# Supplementary dbRDA and variance partitioning for Figure 2.
#
# Response: site-centroid Aitchison distances reconstructed from the exact
# sample geometry used by Panel D and sampling-site identifiers in the metadata.
# Explanatory blocks: environmental PCA, spatial structure, and marine realm.
# Study/BioProject is evaluated as a conditioning sensitivity because its
# contribution may be inseparable from geography in public-data compilations.

options(stringsAsFactors = FALSE, scipen = 999, bitmapType = "cairo")

suppressPackageStartupMessages({
  library(vegan)
  library(adespatial)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(data.table)
  library(scales)
  library(patchwork)
})

## -----------------------------------------------------------------------------
## 1. Configuration
## -----------------------------------------------------------------------------

SEED <- 12345L
PERMUTATIONS <- 9999L
CLR_PSEUDOCOUNT_EXPECTED <- 1
MIN_COMPLETE_SITES <- 50L
ENV_PCA_CUMULATIVE_TARGET <- 0.90
MAX_ENV_PCS <- 8L
MAX_DBMEM <- 8L
MAX_ENV_ARROWS <- 4L

PROJECT_DIR <- Sys.getenv("BIOGEOGRAPHY_BASE", unset = normalizePath(".", mustWork = TRUE))
RESULTS_DIR <- Sys.getenv("BIOGEOGRAPHY_RESULTS", unset = file.path(PROJECT_DIR, "results"))
METADATA_PATH <- Sys.getenv(
  "MAGENTA_METADATA_CSV",
  unset = file.path(PROJECT_DIR, "data", "metadata_taxonomic.csv")
)

PANEL_D_POINTER <- file.path(
  RESULTS_DIR,
  "LATEST_TAXONOMIC_NMDS.txt"
)
PANEL_D_OBJECT_NAME <- file.path("rds", "05_panel_d_analysis_objects.rds")

SCRIPT_ID <- "08"
SCRIPT_BASENAME <- "TAXONOMIC_VARIATION_PARTITIONING"
RUN_NAME <- "vpa"

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
  paste0("LATEST_", SCRIPT_BASENAME, ".txt")
)

LOG_FILE <- file.path(LOG_DIR, "08_dbrda_vpa_supplementary.log")
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

clean_character <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | !nzchar(x) | tolower(x) %in% c("na", "nan", "null")] <-
    NA_character_
  x
}

unique_non_missing <- function(x) {
  x <- clean_character(x)
  sort(unique(x[!is.na(x)]))
}

single_character_or_na <- function(x) {
  values <- unique_non_missing(x)
  if (length(values) == 1L) values else NA_character_
}

median_or_na <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else stats::median(x)
}

range_or_na <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else diff(range(x))
}

resolve_column <- function(df, candidates, label) {
  nms <- names(df)
  hit <- match(tolower(candidates), tolower(nms), nomatch = 0L)
  hit <- hit[hit > 0L]
  if (length(hit) == 0L) {
    stopf(
      "No se encontró %s. Candidatas: %s",
      label,
      paste(candidates, collapse = ", ")
    )
  }
  nms[hit[1L]]
}

make_full_rank <- function(x, prefix) {
  x <- as.matrix(x)
  if (ncol(x) == 0L) {
    return(matrix(numeric(), nrow = nrow(x), ncol = 0L))
  }

  finite_columns <- apply(x, 2L, function(z) all(is.finite(z)))
  variable_columns <- apply(x, 2L, function(z) stats::sd(z) > 1e-10)
  x <- x[, finite_columns & variable_columns, drop = FALSE]
  if (ncol(x) == 0L) {
    return(matrix(numeric(), nrow = nrow(x), ncol = 0L))
  }

  qr_x <- qr(x, tol = 1e-8)
  keep <- sort(qr_x$pivot[seq_len(qr_x$rank)])
  x <- x[, keep, drop = FALSE]
  colnames(x) <- make.names(
    paste0(prefix, "_", colnames(x)),
    unique = TRUE
  )
  x
}

haversine_matrix_km <- function(longitude, latitude) {
  lon <- longitude * pi / 180
  lat <- latitude * pi / 180
  n <- length(lon)
  out <- matrix(0, nrow = n, ncol = n)

  for (i in seq_len(n - 1L)) {
    j <- seq.int(i + 1L, n)
    dlon <- lon[j] - lon[i]
    dlat <- lat[j] - lat[i]
    a <- sin(dlat / 2)^2 + cos(lat[i]) * cos(lat[j]) * sin(dlon / 2)^2
    a <- pmin(1, pmax(0, a))
    distance <- 6371.0088 * 2 * atan2(sqrt(a), sqrt(1 - a))
    out[i, j] <- distance
    out[j, i] <- distance
  }
  out
}

euclidean_coordinates_from_distance <- function(distance_matrix) {
  d <- as.matrix(distance_matrix)
  d2 <- d^2
  gram <- -0.5 * (
    d2 - rowMeans(d2) -
      matrix(colMeans(d2), nrow = nrow(d2), ncol = ncol(d2), byrow = TRUE) +
      mean(d2)
  )

  eig <- eigen(gram, symmetric = TRUE)
  tolerance <- max(abs(eig$values)) * 1e-10
  if (any(eig$values < -tolerance)) {
    stopf(
      "La distancia Aitchison no fue euclidiana: min eigenvalue = %.6g.",
      min(eig$values)
    )
  }

  positive <- eig$values > tolerance
  coordinates <- sweep(
    eig$vectors[, positive, drop = FALSE],
    2L,
    sqrt(eig$values[positive]),
    "*"
  )
  rownames(coordinates) <- rownames(d)

  reconstructed <- as.matrix(stats::dist(coordinates))
  reconstruction_error <- max(abs(reconstructed - d))
  allowed_error <- max(1e-8, max(d) * 1e-7)
  if (!is.finite(reconstruction_error) || reconstruction_error > allowed_error) {
    stopf(
      "Error de reconstrucción Aitchison %.6g > %.6g.",
      reconstruction_error,
      allowed_error
    )
  }

  list(
    coordinates = coordinates,
    reconstruction_error = reconstruction_error,
    n_positive_axes = sum(positive),
    min_eigenvalue = min(eig$values)
  )
}

fit_dbrda_block <- function(
    distance_object,
    explanatory,
    conditioning = NULL,
    model_name,
    permutations = PERMUTATIONS,
    seed = SEED) {

  explanatory <- as.matrix(explanatory)
  if (ncol(explanatory) == 0L) {
    stopf("El modelo %s no tiene predictores estimables.", model_name)
  }

  if (is.null(conditioning) || ncol(as.matrix(conditioning)) == 0L) {
    model <- vegan::dbrda(distance_object ~ explanatory)
    n_conditioning <- 0L
  } else {
    conditioning <- as.matrix(conditioning)
    model <- vegan::dbrda(
      distance_object ~ explanatory + Condition(conditioning)
    )
    n_conditioning <- ncol(conditioning)
  }

  adj_r2 <- suppressWarnings(vegan::RsquareAdj(model)$adj.r.squared)
  raw_r2 <- suppressWarnings(vegan::RsquareAdj(model)$r.squared)
  if (length(adj_r2) != 1L) adj_r2 <- NA_real_
  if (length(raw_r2) != 1L) raw_r2 <- NA_real_

  model_rank <- if (is.null(model$CCA)) 0L else model$CCA$rank
  condition_rank <- if (is.null(model$pCCA)) 0L else model$pCCA$rank
  residual_df <- attr(distance_object, "Size") - 1L - model_rank - condition_rank
  set.seed(seed)
  if (model_rank > 0L && residual_df > 0L) {
    permutation_test <- vegan::anova.cca(model, permutations = permutations)
    fit_status <- "PASS"
    performed_permutations <- permutations
  } else {
    permutation_test <- data.frame(F = NA_real_, check.names = FALSE)
    permutation_test[["Pr(>F)"]] <- NA_real_
    fit_status <- "NOT_ESTIMABLE_ZERO_CONSTRAINED_RANK"
    performed_permutations <- 0L
  }

  statistics <- data.frame(
    model = model_name,
    status = fit_status,
    model_df = model_rank,
    conditioning_df = condition_rank,
    residual_df = residual_df,
    n_sites = attr(distance_object, "Size"),
    n_explanatory_columns = ncol(explanatory),
    n_conditioning_columns = n_conditioning,
    raw_R2 = raw_r2,
    adjusted_R2 = adj_r2,
    pseudo_F = permutation_test$F[1L],
    p_value = permutation_test$`Pr(>F)`[1L],
    permutations = performed_permutations,
    requested_permutations = permutations,
    stringsAsFactors = FALSE
  )

  list(model = model, statistics = statistics)
}

safe_fit_dbrda_block <- function(...) {
  arguments <- list(...)
  requested_model_name <- arguments$model_name
  tryCatch(
    do.call(fit_dbrda_block, arguments),
    error = function(e) list(
      model = NULL,
      statistics = data.frame(
        model = requested_model_name,
        status = paste0("NOT_ESTIMABLE: ", conditionMessage(e)),
        n_sites = NA_integer_,
        n_explanatory_columns = NA_integer_,
        n_conditioning_columns = NA_integer_,
        raw_R2 = NA_real_,
        adjusted_R2 = NA_real_,
        pseudo_F = NA_real_,
        p_value = NA_real_,
        permutations = PERMUTATIONS,
        stringsAsFactors = FALSE
      )
    )
  )
}

## -----------------------------------------------------------------------------
## 2. Load synchronized Panel D objects
## -----------------------------------------------------------------------------

required_files <- c(PANEL_D_POINTER, METADATA_PATH)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stopf("Faltan inputs: %s", paste(missing_files, collapse = "; "))
}

panel_d_dir <- trimws(readLines(PANEL_D_POINTER, warn = FALSE, n = 1L))
if (!nzchar(panel_d_dir) || !dir.exists(panel_d_dir)) {
  stopf("El puntero del Panel D no contiene un directorio válido.")
}

PANEL_D_OBJECT_PATH <- file.path(panel_d_dir, PANEL_D_OBJECT_NAME)
if (!file.exists(PANEL_D_OBJECT_PATH)) {
  stopf("No existe %s", PANEL_D_OBJECT_PATH)
}

panel_d <- readRDS(PANEL_D_OBJECT_PATH)
required_objects <- c("sample_table", "dist_aitchison", "parameters")
missing_objects <- setdiff(required_objects, names(panel_d))
if (length(missing_objects) > 0L) {
  stopf(
    "El objeto Panel D no contiene: %s",
    paste(missing_objects, collapse = ", ")
  )
}

sample_table <- as.data.frame(panel_d$sample_table, stringsAsFactors = FALSE)
required_sample_columns <- c(
  "run_accession",
  "site_id",
  "latitude",
  "longitude",
  "REALM",
  "study_id"
)
missing_sample_columns <- setdiff(required_sample_columns, names(sample_table))
if (length(missing_sample_columns) > 0L) {
  stopf(
    "sample_table no contiene: %s",
    paste(missing_sample_columns, collapse = ", ")
  )
}

sample_table$run_accession <- clean_character(sample_table$run_accession)
sample_table$site_id <- clean_character(sample_table$site_id)
if (anyNA(sample_table$run_accession) || anyDuplicated(sample_table$run_accession)) {
  stopf("run_accession contiene faltantes o duplicados en Panel D.")
}
if (anyNA(sample_table$site_id)) stopf("Panel D contiene sitios no resueltos.")

if (!is.null(panel_d$parameters$clr_pseudocount) &&
    panel_d$parameters$clr_pseudocount != CLR_PSEUDOCOUNT_EXPECTED) {
  stopf(
    "El Panel D usa pseudoconteo %s; se esperaba %s.",
    panel_d$parameters$clr_pseudocount,
    CLR_PSEUDOCOUNT_EXPECTED
  )
}

analysis_ids <- sample_table$run_accession
sample_distance <- as.matrix(panel_d$dist_aitchison)
if (!setequal(rownames(sample_distance), analysis_ids)) {
  stopf("dist_aitchison no coincide con run_accession del Panel D.")
}
sample_distance <- sample_distance[analysis_ids, analysis_ids, drop = FALSE]

log_msg("Panel D object: ", PANEL_D_OBJECT_PATH)
log_msg("Panel D samples: ", nrow(sample_table))
log_msg("Panel D sites: ", dplyr::n_distinct(sample_table$site_id))

## -----------------------------------------------------------------------------
## 3. Resolve the current environmental metadata
## -----------------------------------------------------------------------------

metadata_raw <- data.table::fread(
  METADATA_PATH,
  data.table = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA", "NaN", "NULL")
)
if (!"run_accession" %in% names(metadata_raw)) {
  stopf("La metadata actual no contiene run_accession.")
}
identity_columns <- c("sampling_site_id", "study_id")
if (!all(identity_columns %in% names(metadata_raw))) {
  stopf("La metadata actual requiere sampling_site_id y study_id.")
}
if (anyDuplicated(metadata_raw$run_accession)) {
  stopf("La metadata actual contiene run_accession duplicados.")
}
metadata_rows <- match(sample_table$run_accession, metadata_raw$run_accession)
if (anyNA(metadata_rows)) stopf("Panel D contiene muestras ausentes de la metadata actual.")
current_site_id <- clean_character(metadata_raw$sampling_site_id[metadata_rows])
current_study_id <- clean_character(metadata_raw$study_id[metadata_rows])
if (anyNA(current_site_id) ||
    !identical(sample_table$site_id, current_site_id) ||
    !identical(clean_character(sample_table$study_id), current_study_id)) {
  stopf("Los sitios o estudios del Panel D no coinciden con la metadata actual; repita 04A.")
}

environment_aliases <- list(
  human_footprint_2020 = c(
    "human_footprint_2020",
    "hfp_2020",
    "human_footprint"
  ),
  bio1 = c("bio1", "BIO1"),
  bio2 = c("bio2", "BIO2"),
  bio3 = c("bio3", "BIO3"),
  bio13 = c("bio13", "BIO13"),
  bio15 = c("bio15", "BIO15"),
  gbif_shannon_H_2km = c(
    "gbif_shannon_H_2km",
    "gbif_shannon_2km"
  ),
  gbif_richness_species_2km = c(
    "gbif_richness_species_2km",
    "gbif_species_richness_2km"
  ),
  soc_mean = c("soc_mean", "soil_organic_carbon_mean"),
  phh2o_mean = c("phh2o_mean", "ph_mean", "soil_ph_mean"),
  sand_mean = c("sand_mean", "soil_sand_mean"),
  silt_mean = c("silt_mean", "soil_silt_mean"),
  clay_mean = c("clay_mean", "soil_clay_mean"),
  bdod_mean = c("bdod_mean", "soil_bulk_density_mean")
)

resolved_environment_columns <- vapply(
  names(environment_aliases),
  function(variable) {
    resolve_column(
      metadata_raw,
      environment_aliases[[variable]],
      paste0("la variable ambiental ", variable)
    )
  },
  character(1L)
)

environment_mapping <- data.frame(
  canonical_variable = names(resolved_environment_columns),
  metadata_column = unname(resolved_environment_columns),
  stringsAsFactors = FALSE
)
data.table::fwrite(
  environment_mapping,
  file.path(TABLE_DIR, "08_environmental_variable_mapping.csv")
)

metadata_environment <- data.frame(
  run_accession = clean_character(metadata_raw$run_accession),
  stringsAsFactors = FALSE
)
for (variable in names(resolved_environment_columns)) {
  metadata_environment[[variable]] <- suppressWarnings(as.numeric(
    metadata_raw[[resolved_environment_columns[[variable]]]]
  ))
}
metadata_environment <- metadata_environment[
  !is.na(metadata_environment$run_accession) &
    metadata_environment$run_accession %in% analysis_ids,
  ,
  drop = FALSE
]

metadata_dt <- data.table::as.data.table(metadata_environment)
environment_variables <- names(environment_aliases)

metadata_duplicate_audit <- metadata_dt[, c(
  list(n_metadata_rows = .N),
  lapply(.SD, function(x) {
    x <- x[is.finite(x)]
    data.table::uniqueN(round(x, 8L))
  })
), by = run_accession, .SDcols = environment_variables]

data.table::setnames(
  metadata_duplicate_audit,
  environment_variables,
  paste0("n_unique_", environment_variables)
)
data.table::fwrite(
  metadata_duplicate_audit,
  file.path(TABLE_DIR, "08_metadata_duplicate_environment_audit.csv")
)

conflict_columns <- paste0("n_unique_", environment_variables)
metadata_conflicts <- metadata_duplicate_audit[
  apply(metadata_duplicate_audit[, ..conflict_columns] > 1L, 1L, any)
]
if (nrow(metadata_conflicts) > 0L) {
  conflict_ids <- metadata_conflicts$run_accession

  data.table::fwrite(
    metadata_conflicts,
    file.path(TABLE_DIR, "08_metadata_environment_conflicts_resolved.csv")
  )

  data.table::fwrite(
    metadata_dt[run_accession %in% conflict_ids],
    file.path(TABLE_DIR, "08_metadata_environment_conflicting_rows.csv")
  )

  log_msg(
    "WARNING: ",
    nrow(metadata_conflicts),
    paste0(
      " run_accession tienen valores ambientales discordantes entre filas ",
      "duplicadas. Se resolverán mediante la mediana por run_accession; ",
      "las filas originales y el audit se conservaron en tables/."
    )
  )
} else {
  log_msg("Metadata environmental conflicts: 0")
}

metadata_one <- metadata_dt[, lapply(.SD, median_or_na),
  by = run_accession,
  .SDcols = environment_variables
]

if (anyDuplicated(metadata_one$run_accession)) {
  stopf(
    "La resolución de duplicados falló para %d run_accession.",
    nrow(metadata_conflicts)
  )
}

sample_environment <- sample_table %>%
  dplyr::left_join(
    as.data.frame(metadata_one),
    by = "run_accession"
  )

## -----------------------------------------------------------------------------
## 4. Aggregate environmental metadata and Aitchison geometry by site
## -----------------------------------------------------------------------------

site_identity <- sample_environment %>%
  dplyr::group_by(site_id) %>%
  dplyr::summarise(
    n_metagenomes = dplyr::n(),
    latitude = median_or_na(latitude),
    longitude = median_or_na(longitude),
    n_REALM_values = length(unique_non_missing(REALM)),
    n_studies = length(unique_non_missing(study_id)),
    REALM = single_character_or_na(REALM),
    study_id = single_character_or_na(study_id),
    .groups = "drop"
  )

if (any(site_identity$n_REALM_values != 1L) || anyNA(site_identity$REALM)) {
  data.table::fwrite(
    site_identity %>% dplyr::filter(n_REALM_values != 1L | is.na(REALM)),
    file.path(TABLE_DIR, "08_site_REALM_conflicts_STOP.csv")
  )
  stopf("Existen sitios sin un REALM único.")
}

sample_environment_dt <- data.table::as.data.table(sample_environment)
site_environment_values <- sample_environment_dt[, lapply(.SD, median_or_na),
  by = site_id,
  .SDcols = environment_variables
]
site_environment_ranges <- sample_environment_dt[, lapply(.SD, range_or_na),
  by = site_id,
  .SDcols = environment_variables
]
data.table::setnames(
  site_environment_ranges,
  environment_variables,
  paste0(environment_variables, "_within_site_range")
)

site_table <- site_identity %>%
  dplyr::left_join(as.data.frame(site_environment_values), by = "site_id") %>%
  dplyr::left_join(as.data.frame(site_environment_ranges), by = "site_id")

site_table$study_block <- ifelse(
  site_table$n_studies == 1L & !is.na(site_table$study_id),
  site_table$study_id,
  ifelse(site_table$n_studies > 1L, "MULTIPLE_STUDIES", "UNKNOWN_STUDY")
)
site_table$site_index_all <- seq_len(nrow(site_table))

aitchison_geometry <- euclidean_coordinates_from_distance(sample_distance)
sample_coordinates <- aitchison_geometry$coordinates

site_centroid_coordinates <- do.call(rbind, lapply(site_table$site_id, function(site) {
  rows <- sample_table$site_id == site
  colMeans(sample_coordinates[rows, , drop = FALSE])
}))
rownames(site_centroid_coordinates) <- site_table$site_id

site_distance_all <- as.matrix(stats::dist(site_centroid_coordinates))

complete_environment <- stats::complete.cases(
  site_table[, environment_variables, drop = FALSE]
)
complete_geography <- is.finite(site_table$latitude) &
  is.finite(site_table$longitude) &
  !is.na(site_table$REALM)
analysis_keep <- complete_environment & complete_geography

excluded_sites <- site_table[!analysis_keep, , drop = FALSE]
data.table::fwrite(
  excluded_sites,
  file.path(TABLE_DIR, "08_sites_excluded_missing_environment.csv")
)

analysis_site_table <- as.data.frame(
  site_table[analysis_keep, , drop = FALSE],
  stringsAsFactors = FALSE
)
if (nrow(analysis_site_table) < MIN_COMPLETE_SITES) {
  stopf(
    "Solo %d sitios tienen ambiente completo; se requieren al menos %d.",
    nrow(analysis_site_table),
    MIN_COMPLETE_SITES
  )
}

analysis_indices <- analysis_site_table$site_index_all
site_distance <- stats::as.dist(
  site_distance_all[analysis_indices, analysis_indices, drop = FALSE]
)
if (anyNA(analysis_site_table$site_id) ||
    anyDuplicated(analysis_site_table$site_id)) {
  stopf("site_id contiene faltantes o duplicados antes del dbRDA.")
}
rownames(analysis_site_table) <- analysis_site_table$site_id

data.table::fwrite(
  analysis_site_table,
  file.path(TABLE_DIR, "08_sites_used_dbrda_vpa.csv")
)
data.table::fwrite(
  data.frame(
    n_panel_D_samples = nrow(sample_table),
    n_panel_D_sites = nrow(site_table),
    n_complete_sites = nrow(analysis_site_table),
    n_positive_Aitchison_axes = aitchison_geometry$n_positive_axes,
    min_Aitchison_eigenvalue = aitchison_geometry$min_eigenvalue,
    max_Aitchison_reconstruction_error =
      aitchison_geometry$reconstruction_error
  ),
  file.path(TABLE_DIR, "08_analysis_geometry_audit.csv")
)

log_msg("Complete sites for dbRDA/VPA: ", nrow(analysis_site_table))

## -----------------------------------------------------------------------------
## 5. Environmental PCA and spatial predictors
## -----------------------------------------------------------------------------

environment_raw <- as.matrix(
  analysis_site_table[, environment_variables, drop = FALSE]
)
environment_sd <- apply(environment_raw, 2L, stats::sd)
environment_zero_variance <- names(environment_sd)[
  !is.finite(environment_sd) | environment_sd <= 1e-10
]
if (length(environment_zero_variance) > 0L) {
  data.table::fwrite(
    data.frame(variable = environment_zero_variance),
    file.path(TABLE_DIR, "08_environment_zero_variance_removed.csv")
  )
  environment_raw <- environment_raw[
    ,
    setdiff(colnames(environment_raw), environment_zero_variance),
    drop = FALSE
  ]
}
if (ncol(environment_raw) < 2L) {
  stopf("Quedaron menos de dos variables ambientales con varianza.")
}

environment_scaled <- scale(environment_raw)
environment_pca <- stats::prcomp(
  environment_scaled,
  center = FALSE,
  scale. = FALSE
)

pca_variance <- environment_pca$sdev^2 / sum(environment_pca$sdev^2)
pca_cumulative <- cumsum(pca_variance)
target_axis <- which(pca_cumulative >= ENV_PCA_CUMULATIVE_TARGET)[1L]
if (is.na(target_axis)) target_axis <- length(pca_variance)
n_environment_pcs <- min(MAX_ENV_PCS, target_axis)
n_environment_pcs <- max(1L, n_environment_pcs)

ENV <- environment_pca$x[, seq_len(n_environment_pcs), drop = FALSE]
colnames(ENV) <- paste0("ENV_PC", seq_len(ncol(ENV)))

pca_summary <- data.frame(
  axis = paste0("PC", seq_along(pca_variance)),
  proportion_variance = pca_variance,
  cumulative_variance = pca_cumulative,
  retained = seq_along(pca_variance) <= n_environment_pcs
)
data.table::fwrite(
  pca_summary,
  file.path(TABLE_DIR, "08_environment_PCA_variance.csv")
)

pca_loadings <- as.data.frame(environment_pca$rotation) %>%
  tibble::rownames_to_column("variable")
data.table::fwrite(
  pca_loadings,
  file.path(TABLE_DIR, "08_environment_PCA_loadings.csv")
)

longitude_radians <- analysis_site_table$longitude * pi / 180
geographic_trend <- cbind(
  latitude = as.numeric(scale(analysis_site_table$latitude)),
  longitude_sin = as.numeric(scale(sin(longitude_radians))),
  longitude_cos = as.numeric(scale(cos(longitude_radians)))
)

haversine_distances <- haversine_matrix_km(
  analysis_site_table$longitude,
  analysis_site_table$latitude
)
rownames(haversine_distances) <- analysis_site_table$site_id
colnames(haversine_distances) <- analysis_site_table$site_id

dbmem_raw <- tryCatch(
  adespatial::dbmem(
    stats::as.dist(haversine_distances),
    MEM.autocor = "positive",
    silent = TRUE
  ),
  error = function(e) e
)

if (inherits(dbmem_raw, "error")) {
  stopf("dbMEM not estimable; refusing to omit spatial MEMs: %s", conditionMessage(dbmem_raw))
} else {
  dbmem_eigenvalues <- attr(dbmem_raw, "values")
  data.table::fwrite(data.frame(
    MEM = seq_along(dbmem_eigenvalues), eigenvalue = dbmem_eigenvalues,
    retained = seq_along(dbmem_eigenvalues) <= MAX_DBMEM
  ), file.path(TABLE_DIR, "08_dbMEM_eigenvalues.csv"))
  dbmem_raw <- as.matrix(dbmem_raw)
  stopifnot(ncol(dbmem_raw) > 0L)
  data.table::fwrite(data.frame(
    n_positive_MEMs = ncol(dbmem_raw),
    max_abs_column_mean = max(abs(colMeans(dbmem_raw))),
    max_abs_weighted_norm_error = max(abs(colSums(dbmem_raw^2) / nrow(dbmem_raw) - 1)),
    threshold_km = adespatial::give.thresh(stats::as.dist(haversine_distances)),
    implementation = paste0("adespatial ", packageVersion("adespatial"), " official dbmem(dist, positive)")
  ), file.path(TABLE_DIR, "08_dbMEM_construction_audit.csv"))
  n_dbmem <- min(MAX_DBMEM, ncol(dbmem_raw))
  DBMEM <- dbmem_raw[, seq_len(n_dbmem), drop = FALSE]
  colnames(DBMEM) <- paste0("dbMEM", seq_len(ncol(DBMEM)))
}

SPACE <- make_full_rank(
  cbind(geographic_trend, DBMEM),
  prefix = "SPACE"
)
if (ncol(SPACE) == 0L) stopf("No se generaron predictores espaciales.")

realm_factor <- droplevels(factor(analysis_site_table$REALM))
if (nlevels(realm_factor) < 2L) stopf("REALM tiene menos de dos niveles.")
REALM <- stats::model.matrix(~ realm_factor)[, -1L, drop = FALSE]
REALM <- make_full_rank(REALM, prefix = "REALM")

study_factor <- droplevels(factor(analysis_site_table$study_block))
if (nlevels(study_factor) >= 2L) {
  STUDY <- stats::model.matrix(~ study_factor)[, -1L, drop = FALSE]
  STUDY <- make_full_rank(STUDY, prefix = "STUDY")
} else {
  STUDY <- matrix(numeric(), nrow = nrow(analysis_site_table), ncol = 0L)
}

rownames(ENV) <- analysis_site_table$site_id
rownames(SPACE) <- analysis_site_table$site_id
rownames(REALM) <- analysis_site_table$site_id
rownames(STUDY) <- analysis_site_table$site_id

predictor_block_summary <- data.frame(
  block = c("Environment_PCA", "Spatial", "Marine_realm", "Study_sensitivity"),
  n_columns = c(ncol(ENV), ncol(SPACE), ncol(REALM), ncol(STUDY)),
  stringsAsFactors = FALSE
)
data.table::fwrite(
  predictor_block_summary,
  file.path(TABLE_DIR, "08_predictor_block_summary.csv")
)

rank_diagnostics <- data.frame(
  matrix = c(
    "ENV",
    "SPACE",
    "REALM",
    "STUDY",
    "ENV_SPACE_REALM",
    "ENV_SPACE_REALM_STUDY"
  ),
  n_columns = c(
    ncol(ENV),
    ncol(SPACE),
    ncol(REALM),
    ncol(STUDY),
    ncol(cbind(ENV, SPACE, REALM)),
    ncol(cbind(ENV, SPACE, REALM, STUDY))
  ),
  rank = c(
    qr(cbind(1, ENV))$rank - 1L,
    qr(cbind(1, SPACE))$rank - 1L,
    qr(cbind(1, REALM))$rank - 1L,
    if (ncol(STUDY) > 0L) qr(cbind(1, STUDY))$rank - 1L else 0L,
    qr(cbind(1, ENV, SPACE, REALM))$rank - 1L,
    qr(cbind(1, ENV, SPACE, REALM, STUDY))$rank - 1L
  )
)
rank_diagnostics$rank_deficiency <-
  rank_diagnostics$n_columns - rank_diagnostics$rank
data.table::fwrite(
  rank_diagnostics,
  file.path(TABLE_DIR, "08_predictor_rank_and_confounding_audit.csv")
)

realm_study_table <- as.data.frame.matrix(table(
  REALM = analysis_site_table$REALM,
  Study = analysis_site_table$study_block
)) %>%
  tibble::rownames_to_column("REALM")
data.table::fwrite(
  realm_study_table,
  file.path(TABLE_DIR, "08_REALM_by_Study_site_counts.csv")
)

## -----------------------------------------------------------------------------
## 6. Formal three-block VPA and partial dbRDA tests
## -----------------------------------------------------------------------------

FULL <- cbind(ENV, SPACE, REALM)

full_fit <- fit_dbrda_block(
  site_distance,
  FULL,
  model_name = "Full_ENV_SPACE_REALM",
  seed = SEED + 1L
)
environment_unique_fit <- fit_dbrda_block(
  site_distance,
  ENV,
  conditioning = cbind(SPACE, REALM),
  model_name = "Environment_unique_given_SPACE_REALM",
  seed = SEED + 2L
)
space_unique_fit <- fit_dbrda_block(
  site_distance,
  SPACE,
  conditioning = cbind(ENV, REALM),
  model_name = "Space_unique_given_ENV_REALM",
  seed = SEED + 3L
)
realm_unique_fit <- fit_dbrda_block(
  site_distance,
  REALM,
  conditioning = cbind(ENV, SPACE),
  model_name = "REALM_unique_given_ENV_SPACE",
  seed = SEED + 4L
)

primary_model_statistics <- dplyr::bind_rows(
  full_fit$statistics,
  environment_unique_fit$statistics,
  space_unique_fit$statistics,
  realm_unique_fit$statistics
) %>%
  dplyr::mutate(
    p_BH_unique_components = c(
      NA_real_,
      stats::p.adjust(p_value[-1L], method = "BH")
    )
  )

data.table::fwrite(
  primary_model_statistics,
  file.path(TABLE_DIR, "08_dbrda_primary_and_unique_tests.csv")
)

variance_partition <- vegan::varpart(
  site_distance,
  ENV,
  SPACE,
  REALM
)

writeLines(
  capture.output(print(variance_partition, digits = 5)),
  file.path(LOG_DIR, "08_varpart_full_output.txt")
)

individual_fractions <- as.data.frame(
  variance_partition$part$indfract,
  check.names = FALSE
) %>%
  tibble::rownames_to_column("fraction_raw") %>%
  dplyr::mutate(
    fraction_code = sub("^(\\[[a-z]\\]).*$", "\\1", fraction_raw),
    fraction_label = dplyr::recode(
      fraction_code,
      "[a]" = "Unique environment",
      "[b]" = "Unique spatial structure",
      "[c]" = "Unique marine realm",
      "[d]" = "Shared environment + space",
      "[e]" = "Shared space + realm",
      "[f]" = "Shared environment + realm",
      "[g]" = "Shared environment + space + realm",
      "[h]" = "Unexplained",
      .default = fraction_raw
    )
  )

data.table::fwrite(
  individual_fractions,
  file.path(TABLE_DIR, "08_variance_partition_individual_fractions.csv")
)

## -----------------------------------------------------------------------------
## 7. Study/BioProject-conditioned sensitivity
## -----------------------------------------------------------------------------

study_sensitivity_fits <- list()
if (ncol(STUDY) > 0L) {
  study_sensitivity_fits <- list(
    safe_fit_dbrda_block(
      site_distance,
      FULL,
      conditioning = STUDY,
      model_name = "Full_ENV_SPACE_REALM_condition_STUDY",
      seed = SEED + 11L
    ),
    safe_fit_dbrda_block(
      site_distance,
      ENV,
      conditioning = cbind(SPACE, REALM, STUDY),
      model_name = "Environment_unique_condition_SPACE_REALM_STUDY",
      seed = SEED + 12L
    ),
    safe_fit_dbrda_block(
      site_distance,
      SPACE,
      conditioning = cbind(ENV, REALM, STUDY),
      model_name = "Space_unique_condition_ENV_REALM_STUDY",
      seed = SEED + 13L
    ),
    safe_fit_dbrda_block(
      site_distance,
      REALM,
      conditioning = cbind(ENV, SPACE, STUDY),
      model_name = "REALM_unique_condition_ENV_SPACE_STUDY",
      seed = SEED + 14L
    ),
    safe_fit_dbrda_block(
      site_distance,
      cbind(ENV, SPACE, STUDY),
      model_name = "Alternative_full_ENV_SPACE_STUDY_without_REALM",
      seed = SEED + 15L
    )
  )
}

if (length(study_sensitivity_fits) > 0L) {
  study_sensitivity_statistics <- dplyr::bind_rows(lapply(
    study_sensitivity_fits,
    `[[`,
    "statistics"
  ))
} else {
  study_sensitivity_statistics <- data.frame(
    model = "Study_not_available",
    status = "NOT_ESTIMABLE",
    n_sites = nrow(analysis_site_table),
    n_explanatory_columns = NA_integer_,
    n_conditioning_columns = 0L,
    raw_R2 = NA_real_,
    adjusted_R2 = NA_real_,
    pseudo_F = NA_real_,
    p_value = NA_real_,
    permutations = PERMUTATIONS,
    stringsAsFactors = FALSE
  )
}

data.table::fwrite(
  study_sensitivity_statistics,
  file.path(TABLE_DIR, "08_Study_conditioned_sensitivity.csv")
)

## -----------------------------------------------------------------------------
## 8. Supplementary dbRDA ordination without ellipses or internal titles
## -----------------------------------------------------------------------------

ordination_model <- full_fit$model
site_scores <- vegan::scores(
  ordination_model,
  display = "sites",
  choices = 1:2,
  scaling = 2
)
site_scores <- as.data.frame(site_scores) %>%
  tibble::rownames_to_column("site_id")
colnames(site_scores)[2:3] <- c("dbRDA1", "dbRDA2")
site_scores <- site_scores %>%
  dplyr::left_join(
    analysis_site_table %>%
      dplyr::select(
        site_id,
        REALM,
        latitude,
        longitude,
        study_block,
        n_metagenomes
      ),
    by = "site_id"
  )

constrained_eigenvalues <- vegan::eigenvals(
  ordination_model,
  model = "constrained"
)
axis_percent <- 100 * constrained_eigenvalues / ordination_model$tot.chi
axis1_percent <- if (length(axis_percent) >= 1L) axis_percent[1L] else NA_real_
axis2_percent <- if (length(axis_percent) >= 2L) axis_percent[2L] else NA_real_

set.seed(SEED + 20L)
environment_score_order <- match(
  site_scores$site_id,
  analysis_site_table$site_id
)
if (anyNA(environment_score_order)) {
  stopf(
    "No se pudieron alinear %d sitios del dbRDA con la matriz ambiental.",
    sum(is.na(environment_score_order))
  )
}
environment_for_envfit <- environment_scaled[
  environment_score_order,
  ,
  drop = FALSE
]
rownames(environment_for_envfit) <- site_scores$site_id

environment_fit <- vegan::envfit(
  as.matrix(site_scores[, c("dbRDA1", "dbRDA2")]),
  environment_for_envfit,
  permutations = PERMUTATIONS
)

environment_vector_scores <- vegan::scores(
  environment_fit,
  display = "vectors"
)
environment_vectors <- data.frame(
  variable = rownames(environment_vector_scores),
  dbRDA1 = environment_vector_scores[, 1L],
  dbRDA2 = environment_vector_scores[, 2L],
  r_squared = unname(environment_fit$vectors$r),
  p_value = unname(environment_fit$vectors$pvals),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(p_BH = stats::p.adjust(p_value, method = "BH"))

pretty_environment_labels <- c(
  human_footprint_2020 = "Human footprint",
  bio1 = "Annual mean temperature",
  bio2 = "Mean diurnal temperature range",
  bio3 = "Isothermality",
  bio13 = "Wettest-month precipitation",
  bio15 = "Precipitation seasonality",
  gbif_shannon_H_2km = "Plant Shannon diversity",
  gbif_richness_species_2km = "Plant species richness",
  soc_mean = "Soil organic carbon",
  phh2o_mean = "Soil pH",
  sand_mean = "Soil sand",
  silt_mean = "Soil silt",
  clay_mean = "Soil clay",
  bdod_mean = "Soil bulk density"
)

short_environment_labels <- c(
  human_footprint_2020 = "HFP",
  bio1 = "BIO1",
  bio2 = "BIO2",
  bio3 = "BIO3",
  bio13 = "BIO13",
  bio15 = "BIO15",
  gbif_shannon_H_2km = "Plant H'",
  gbif_richness_species_2km = "Plant richness",
  soc_mean = "SOC",
  phh2o_mean = "Soil pH",
  sand_mean = "Sand",
  silt_mean = "Silt",
  clay_mean = "Clay",
  bdod_mean = "BDOD"
)

environment_vectors$label_full <- unname(pretty_environment_labels[
  environment_vectors$variable
])
environment_vectors$label_full[is.na(environment_vectors$label_full)] <-
  environment_vectors$variable[is.na(environment_vectors$label_full)]
environment_vectors$label <- unname(short_environment_labels[
  environment_vectors$variable
])
environment_vectors$label[is.na(environment_vectors$label)] <-
  environment_vectors$variable[is.na(environment_vectors$label)]

environment_vectors_plot <- environment_vectors %>%
  dplyr::filter(p_BH <= 0.05) %>%
  dplyr::arrange(dplyr::desc(r_squared)) %>%
  dplyr::slice_head(n = MAX_ENV_ARROWS) %>%
  dplyr::mutate(
    arrow_x = dbRDA1,
    arrow_y = dbRDA2
  )

if (nrow(environment_vectors_plot) > 0L) {
  x_span <- diff(stats::quantile(
    site_scores$dbRDA1,
    probs = c(0.05, 0.95),
    na.rm = TRUE,
    names = FALSE
  ))
  y_span <- diff(stats::quantile(
    site_scores$dbRDA2,
    probs = c(0.05, 0.95),
    na.rm = TRUE,
    names = FALSE
  ))
  x_denominator <- max(abs(environment_vectors_plot$dbRDA1))
  y_denominator <- max(abs(environment_vectors_plot$dbRDA2))
  scale_candidates <- c(
    if (x_denominator > 0) 0.34 * x_span / x_denominator else Inf,
    if (y_denominator > 0) 0.34 * y_span / y_denominator else Inf
  )
  arrow_multiplier <- min(scale_candidates[is.finite(scale_candidates)])
  environment_vectors_plot <- environment_vectors_plot %>%
    dplyr::mutate(
      arrow_x = dbRDA1 * arrow_multiplier,
      arrow_y = dbRDA2 * arrow_multiplier
    )
}

data.table::fwrite(
  environment_vectors,
  file.path(TABLE_DIR, "08_dbRDA_environment_vectors_all.csv")
)
data.table::fwrite(
  environment_vectors_plot,
  file.path(TABLE_DIR, "08_dbRDA_environment_vectors_plotted.csv")
)

score_center <- c(
  dbRDA1 = stats::median(site_scores$dbRDA1, na.rm = TRUE),
  dbRDA2 = stats::median(site_scores$dbRDA2, na.rm = TRUE)
)
score_scale <- c(
  dbRDA1 = stats::mad(site_scores$dbRDA1, constant = 1, na.rm = TRUE),
  dbRDA2 = stats::mad(site_scores$dbRDA2, constant = 1, na.rm = TRUE)
)
score_scale[!is.finite(score_scale) | score_scale <= 1e-10] <- 1
site_scores$robust_ordination_distance <- sqrt(
  ((site_scores$dbRDA1 - score_center[["dbRDA1"]]) /
      score_scale[["dbRDA1"]])^2 +
    ((site_scores$dbRDA2 - score_center[["dbRDA2"]]) /
      score_scale[["dbRDA2"]])^2
)
outlier_cutoff <- unname(stats::quantile(
  site_scores$robust_ordination_distance,
  probs = 0.975,
  na.rm = TRUE
))
site_scores$ordination_outlier_top_2_5pct <-
  site_scores$robust_ordination_distance >= outlier_cutoff
data.table::fwrite(
  site_scores %>%
    dplyr::arrange(dplyr::desc(robust_ordination_distance)),
  file.path(TABLE_DIR, "08_dbRDA_site_score_outlier_audit.csv")
)

realm_colors <- c(
  "Central Indo-Pacific" = "#66C2A5",
  "Indo-Pacific" = "#FC8D62",
  "Temperate Northern Atlantic" = "#8DA0CB",
  "Tropical Atlantic" = "#E78AC3",
  "Tropical Eastern Pacific" = "#A6D854",
  "Western Indo-Pacific" = "#FFD92F"
)

# The eigenaxes are swapped only for display so that the larger score range is
# horizontal. This is a rigid 90-degree rotation; no scores are transformed.
p_dbrda <- ggplot(site_scores, aes(x = dbRDA2, y = dbRDA1)) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.3,
    linetype = "dashed",
    color = "grey80"
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.3,
    linetype = "dashed",
    color = "grey80"
  ) +
  geom_point(
    aes(fill = REALM, color = REALM),
    shape = 21,
    size = 2.35,
    stroke = 0.35,
    alpha = 0.76
  ) +
  geom_segment(
    data = environment_vectors_plot,
    aes(x = 0, y = 0, xend = arrow_y, yend = arrow_x),
    inherit.aes = FALSE,
    color = "grey20",
    linewidth = 0.4,
    arrow = grid::arrow(length = grid::unit(0.11, "cm"), type = "closed")
  ) +
  ggrepel::geom_text_repel(
    data = environment_vectors_plot,
    aes(x = arrow_y, y = arrow_x, label = label),
    inherit.aes = FALSE,
    size = 3.2,
    color = "grey15",
    min.segment.length = 0,
    seed = SEED,
    box.padding = 0.32,
    point.padding = 0.16,
    force = 2,
    max.overlaps = Inf
  ) +
  scale_fill_manual(
    values = realm_colors,
    drop = FALSE,
    name = "Marine realm",
    guide = guide_legend(
      nrow = 3,
      byrow = TRUE,
      override.aes = list(size = 2.8, alpha = 1)
    )
  ) +
  scale_color_manual(values = realm_colors, drop = FALSE, guide = "none") +
  labs(
    x = sprintf("dbRDA2 (%.1f%% of total inertia)", axis2_percent),
    y = sprintf("dbRDA1 (%.1f%% of total inertia)", axis1_percent)
  ) +
  coord_equal(clip = "off") +
  theme_classic(base_size = 11) +
  theme(
    axis.text = element_text(color = "black"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(face = "bold"),
    legend.box.margin = margin(t = 3, r = 0, b = 0, l = 0),
    plot.margin = margin(7, 9, 5, 7)
  )

## -----------------------------------------------------------------------------
## 9. Supplementary VPA plot without internal title or panel letter
## -----------------------------------------------------------------------------

adjusted_r2_column <- grep(
  "Adj\\.R\\.square|Adjusted",
  names(individual_fractions),
  value = TRUE
)[1L]
if (is.na(adjusted_r2_column)) {
  stopf("No se encontró la columna de R2 ajustado en varpart$indfract.")
}

vpa_plot_data <- individual_fractions %>%
  dplyr::filter(fraction_code %in% paste0("[", letters[1:7], "]")) %>%
  dplyr::mutate(
    adjusted_R2 = .data[[adjusted_r2_column]],
    fraction_label = factor(fraction_label, levels = rev(c(
      "Unique environment",
      "Unique spatial structure",
      "Unique marine realm",
      "Shared environment + space",
      "Shared environment + realm",
      "Shared space + realm",
      "Shared environment + space + realm"
    )))
  )

total_adjusted_explained <- sum(
  vpa_plot_data$adjusted_R2,
  na.rm = TRUE
)

stopifnot(abs(total_adjusted_explained - full_fit$statistics$adjusted_R2) < 1e-8)

fraction_colors <- c(
  "Unique environment" = "#D9A441",
  "Unique spatial structure" = "#4C9FBE",
  "Unique marine realm" = "#4FA36C",
  "Shared environment + space" = "#C8B65A",
  "Shared environment + realm" = "#94A95E",
  "Shared space + realm" = "#4C9B8F",
  "Shared environment + space + realm" = "#8C8C8C"
)

p_vpa <- ggplot(
  vpa_plot_data,
  aes(x = adjusted_R2, y = fraction_label, fill = fraction_label)
) +
  geom_vline(xintercept = 0, linewidth = 0.35, color = "grey45") +
  geom_col(width = 0.68, color = "grey20", linewidth = 0.25) +
  geom_text(
    aes(
      label = scales::percent(adjusted_R2, accuracy = 0.1),
      hjust = ifelse(adjusted_R2 >= 0, -0.08, 1.08)
    ),
    size = 3.0
  ) +
  scale_fill_manual(values = fraction_colors, guide = "none") +
  scale_x_continuous(
    labels = scales::label_percent(accuracy = 1),
    expand = expansion(mult = c(0.16, 0.18))
  ) +
  labs(
    x = "Adjusted variance fraction",
    y = NULL,
    caption = sprintf(
      "Total adjusted explained variance: %.1f%%",
      100 * total_adjusted_explained
    )
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 11) +
  theme(
    axis.text = element_text(color = "black"),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.caption = element_text(hjust = 1, color = "grey30"),
    plot.margin = margin(6, 15, 6, 7)
  )

## -----------------------------------------------------------------------------
## 10. Save supplementary figures and synchronized objects
## -----------------------------------------------------------------------------

ggsave(
  file.path(PLOT_DIR, "08_Supplementary_dbRDA_site_Aitchison.png"),
  p_dbrda,
  width = 6.8,
  height = 5.2,
  units = "in",
  dpi = 600,
  bg = "white"
)
ggsave(
  file.path(PLOT_DIR, "08_Supplementary_dbRDA_site_Aitchison.pdf"),
  p_dbrda,
  width = 6.8,
  height = 5.2,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)

ggsave(
  file.path(PLOT_DIR, "08_Supplementary_variance_partition.png"),
  p_vpa,
  width = 6.8,
  height = 4.4,
  units = "in",
  dpi = 600,
  bg = "white"
)
ggsave(
  file.path(PLOT_DIR, "08_Supplementary_variance_partition.pdf"),
  p_vpa,
  width = 6.8,
  height = 4.4,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)

p_combined <- p_dbrda / p_vpa +
  patchwork::plot_layout(heights = c(1.2, 1))

ggsave(
  file.path(PLOT_DIR, "08_Supplementary_dbRDA_VPA_combined.png"),
  p_combined,
  width = 6.8,
  height = 9.3,
  units = "in",
  dpi = 600,
  bg = "white"
)
ggsave(
  file.path(PLOT_DIR, "08_Supplementary_dbRDA_VPA_combined.pdf"),
  p_combined,
  width = 6.8,
  height = 9.3,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)

run_summary <- data.frame(
  status = "PASS",
  run_name = RUN_NAME,
  panel_D_object = PANEL_D_OBJECT_PATH,
  metadata_path = METADATA_PATH,
  n_panel_D_samples = nrow(sample_table),
  n_panel_D_sites = nrow(site_table),
  n_run_accession_with_environmental_conflicts = nrow(metadata_conflicts),
  duplicate_environment_resolution = "median_within_run_accession",
  n_complete_sites = nrow(analysis_site_table),
  n_environmental_variables = ncol(environment_raw),
  n_environment_PCs = ncol(ENV),
  environment_PCA_cumulative_variance =
    pca_cumulative[n_environment_pcs],
  n_spatial_predictors = ncol(SPACE),
  n_dbMEM = ncol(DBMEM),
  n_REALM_predictors = ncol(REALM),
  n_Study_predictors = ncol(STUDY),
  dbRDA_display_orientation = "x=dbRDA2; y=dbRDA1",
  permutations = PERMUTATIONS,
  seed = SEED,
  adespatial_version = as.character(packageVersion("adespatial")),
  vegan_version = as.character(packageVersion("vegan")),
  stringsAsFactors = FALSE
)
data.table::fwrite(
  run_summary,
  file.path(TABLE_DIR, "08_run_summary.csv")
)

saveRDS(
  list(
    sample_table = sample_table,
    site_table = site_table,
    analysis_site_table = analysis_site_table,
    site_Aitchison_distance = site_distance,
    environment_scaled = environment_scaled,
    environment_PCA = environment_pca,
    ENV = ENV,
    SPACE = SPACE,
    REALM = REALM,
    STUDY = STUDY,
    variance_partition = variance_partition,
    primary_models = list(
      full = full_fit$model,
      environment_unique = environment_unique_fit$model,
      space_unique = space_unique_fit$model,
      realm_unique = realm_unique_fit$model
    ),
    primary_model_statistics = primary_model_statistics,
    study_sensitivity_statistics = study_sensitivity_statistics,
    dbRDA_site_scores = site_scores,
    dbRDA_environment_vectors = environment_vectors,
    parameters = list(
      seed = SEED,
      permutations = PERMUTATIONS,
      environment_PCA_target = ENV_PCA_CUMULATIVE_TARGET,
      max_environment_PCs = MAX_ENV_PCS,
      max_dbMEM = MAX_DBMEM
    )
  ),
  file.path(RDS_DIR, "08_DBRDA_VPA_supplementary_objects.rds")
)

input_paths <- c(PANEL_D_POINTER, PANEL_D_OBJECT_PATH, METADATA_PATH)
input_provenance <- data.frame(
  input = c("Panel_D_pointer", "Panel_D_object", "Current_metadata"),
  path = input_paths,
  size_bytes = as.numeric(file.info(input_paths)$size),
  modified_time = as.character(file.info(input_paths)$mtime),
  stringsAsFactors = FALSE
)
data.table::fwrite(
  input_provenance,
  file.path(TABLE_DIR, "08_input_provenance.csv")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(LOG_DIR, "08_sessionInfo.txt")
)

writeLines(c(
  "Taxonomic variation partitioning across metadata-defined sampling sites.",
  "Site identity, geometry and samples read from the synchronized Panel D object.",
  "Official adespatial::dbmem retained; package version recorded in run summary and sessionInfo.",
  "Shared fractions follow vegan: [e] is space + realm and [f] is environment + realm.",
  "All negative adjusted variance fractions retained, without truncation.",
  "Study-conditioned terms with zero estimable constrained rank are explicitly marked not estimable.",
  "Permutation tests use site-level unrestricted permutations; spatial/study confounding restricts interpretation.",
  "Study-conditioned models are sensitivity analyses, not proof that all geographic confounding is removed."
), file.path(OUT_DIR, "README_analysis.txt"))

log_msg("Primary full adjusted R2: ", sprintf("%.6f", full_fit$statistics$adjusted_R2))
log_msg("Supplementary figures: ", PLOT_DIR)
writeLines(normalizePath(OUT_DIR, mustWork = FALSE), LATEST_POINTER)
log_msg("=== SCRIPT COMPLETED SUCCESSFULLY ===")
