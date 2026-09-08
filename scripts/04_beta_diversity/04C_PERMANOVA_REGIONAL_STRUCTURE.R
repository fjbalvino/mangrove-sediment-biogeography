#!/usr/bin/env Rscript
# -*- coding: UTF-8 -*-

# Corrected site-level regional PERMANOVA for Figure 2d and Supplementary results.
#
# The script inherits the exact samples and distance matrices produced by
# 05_NMDS_AITCHISON_PANEL_D_SET18.R. The primary unit is the coordinate-defined
# site. Site Aitchison centroids are reconstructed exactly from the Euclidean
# sample-distance geometry. A deterministic within-site medoid is used only for
# the Bray-Curtis sensitivity analysis. Aitchison is the primary metric.

options(stringsAsFactors = FALSE, scipen = 999, bitmapType = "cairo")

suppressPackageStartupMessages({
  library(vegan)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(data.table)
  library(scales)
})

## -----------------------------------------------------------------------------
## 1. Configuration
## -----------------------------------------------------------------------------

SEED <- 4242L
PRIMARY_PERMUTATIONS <- 9999L
N_BALANCED_ITERATIONS <- 199L
PERMUTATIONS_PER_ITERATION <- 999L
LOSO_PERMUTATIONS <- 999L
LOCO_PERMUTATIONS <- 999L
STUDY_BLOCK_PERMUTATIONS <- 9999L
MIN_SITES_PER_LEVEL <- 5L

EXPECTED_ANALYSIS_SAMPLES <- 314L
EXPECTED_COORDINATE_SITES <- 130L

PROJECT_DIR <- Sys.getenv("BIOGEOGRAPHY_BASE", unset = "/home/fjbalvino/Biogeografia_2026")
RESULTS_DIR <- file.path(PROJECT_DIR, "resultados")

PANEL_D_POINTER <- file.path(
  RESULTS_DIR,
  "LATEST_05_NMDS_AITCHISON_PANEL_D_SET18.txt"
)
PANEL_D_OBJECT_NAME <- file.path("rds", "05_panel_d_analysis_objects.rds")

SCRIPT_ID <- "07"
SCRIPT_BASENAME <- "PERMANOVA_REGIONAL_SET18_CORRECTED"
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

LOG_FILE <- file.path(LOG_DIR, "07_permanova_regional.log")
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
  x[is.na(x) | !nzchar(x) | tolower(x) %in% c("na", "nan", "null")] <- NA_character_
  x
}

string_seed <- function(...) {
  z <- paste(..., collapse = "|")
  SEED + sum(utf8ToInt(z))
}

as_dist_subset <- function(distance_matrix, indices) {
  stats::as.dist(distance_matrix[indices, indices, drop = FALSE])
}

unique_non_missing <- function(x) {
  sort(unique(clean_character(x)[!is.na(clean_character(x))]))
}

single_value_or_na <- function(x) {
  values <- unique_non_missing(x)
  if (length(values) == 1L) values else NA_character_
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
      "La distancia Aitchison no es euclidiana dentro de la tolerancia: min eigenvalue = %.6g.",
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
      "No se pudo reconstruir la geometría Aitchison: error máximo %.6g > %.6g.",
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

## -----------------------------------------------------------------------------
## 2. Load the synchronized panel D analysis objects
## -----------------------------------------------------------------------------

if (!file.exists(PANEL_D_POINTER)) {
  stopf(
    "No existe %s. Ejecute primero 05_NMDS_AITCHISON_PANEL_D_SET18.R.",
    PANEL_D_POINTER
  )
}

panel_d_dir <- trimws(readLines(PANEL_D_POINTER, warn = FALSE, n = 1L))
if (!nzchar(panel_d_dir) || !dir.exists(panel_d_dir)) {
  stopf("El puntero del panel D no contiene un directorio válido.")
}

PANEL_D_OBJECT_PATH <- file.path(panel_d_dir, PANEL_D_OBJECT_NAME)
if (!file.exists(PANEL_D_OBJECT_PATH)) {
  stopf("No existe el objeto sincronizado: %s", PANEL_D_OBJECT_PATH)
}

panel_d <- readRDS(PANEL_D_OBJECT_PATH)
required_objects <- c("sample_table", "dist_aitchison", "dist_bray", "parameters")
missing_objects <- setdiff(required_objects, names(panel_d))
if (length(missing_objects) > 0L) {
  stopf("Faltan elementos en el objeto panel D: %s", paste(missing_objects, collapse = ", "))
}

sample_table <- as.data.frame(panel_d$sample_table, stringsAsFactors = FALSE)
required_columns <- c(
  "run_accession",
  "site_id",
  "REALM",
  "PROVINCE",
  "Ecoregion",
  "Country",
  "study_id"
)
missing_columns <- setdiff(required_columns, names(sample_table))
if (length(missing_columns) > 0L) {
  stopf("Faltan columnas en sample_table: %s", paste(missing_columns, collapse = ", "))
}

sample_table$run_accession <- clean_character(sample_table$run_accession)
if (anyNA(sample_table$run_accession) || anyDuplicated(sample_table$run_accession)) {
  stopf("sample_table contiene run_accession faltantes o duplicados.")
}
if (anyNA(sample_table$site_id)) stopf("sample_table contiene site_id faltantes.")
if (nrow(sample_table) != EXPECTED_ANALYSIS_SAMPLES) {
  stopf(
    "Se esperaban %d muestras del panel D; se encontraron %d.",
    EXPECTED_ANALYSIS_SAMPLES,
    nrow(sample_table)
  )
}
if (dplyr::n_distinct(sample_table$site_id) != EXPECTED_COORDINATE_SITES) {
  stopf(
    "Se esperaban %d sitios del panel D; se encontraron %d.",
    EXPECTED_COORDINATE_SITES,
    dplyr::n_distinct(sample_table$site_id)
  )
}

analysis_ids <- sample_table$run_accession
sample_table$sample_index <- seq_len(nrow(sample_table))

aitchison_matrix <- as.matrix(panel_d$dist_aitchison)
bray_matrix <- as.matrix(panel_d$dist_bray)

for (matrix_name in c("aitchison_matrix", "bray_matrix")) {
  x <- get(matrix_name)
  if (!setequal(rownames(x), analysis_ids) || !setequal(colnames(x), analysis_ids)) {
    stopf("%s no contiene exactamente los run_accession del panel D.", matrix_name)
  }
  x <- x[analysis_ids, analysis_ids, drop = FALSE]
  assign(matrix_name, x)
}

log_msg("Panel D object: ", PANEL_D_OBJECT_PATH)
log_msg("Samples inherited from panel D: n = ", nrow(sample_table))
log_msg("Coordinate-defined localities: n = ", dplyr::n_distinct(sample_table$site_id))

## -----------------------------------------------------------------------------
## 3. Build independent site-level analysis units
## -----------------------------------------------------------------------------

site_table <- sample_table %>%
  dplyr::group_by(site_id) %>%
  dplyr::summarise(
    n_metagenomes = dplyr::n(),
    run_accessions = paste(sort(run_accession), collapse = " | "),
    n_REALM_values = length(unique_non_missing(REALM)),
    n_PROVINCE_values = length(unique_non_missing(PROVINCE)),
    n_Ecoregion_values = length(unique_non_missing(Ecoregion)),
    n_Country_values = length(unique_non_missing(Country)),
    n_studies = length(unique_non_missing(study_id)),
    REALM = single_value_or_na(REALM),
    PROVINCE = single_value_or_na(PROVINCE),
    Ecoregion = single_value_or_na(Ecoregion),
    Country = single_value_or_na(Country),
    study_id = single_value_or_na(study_id),
    .groups = "drop"
  )

site_geography_conflicts <- site_table %>%
  dplyr::filter(
    n_REALM_values > 1L |
      n_PROVINCE_values > 1L |
      n_Ecoregion_values > 1L |
      n_Country_values > 1L
  )
if (nrow(site_geography_conflicts) > 0L) {
  data.table::fwrite(
    site_geography_conflicts,
    file.path(TABLE_DIR, "07_site_geography_conflicts_STOP.csv")
  )
  stopf(
    "Hay %d sitios con asignaciones geográficas incompatibles.",
    nrow(site_geography_conflicts)
  )
}

site_table$site_index <- seq_len(nrow(site_table))

aitchison_geometry <- euclidean_coordinates_from_distance(aitchison_matrix)
sample_coordinates <- aitchison_geometry$coordinates

site_centroid_coordinates <- do.call(rbind, lapply(site_table$site_id, function(site) {
  sample_rows <- sample_table$site_id == site
  colMeans(sample_coordinates[sample_rows, , drop = FALSE])
}))
rownames(site_centroid_coordinates) <- site_table$site_id

site_aitchison_matrix <- as.matrix(stats::dist(site_centroid_coordinates))

bray_medoid_indices <- vapply(site_table$site_id, function(site) {
  candidates <- which(sample_table$site_id == site)
  if (length(candidates) == 1L) return(candidates)

  within_site <- bray_matrix[candidates, candidates, drop = FALSE]
  mean_distance <- rowMeans(within_site)
  best <- candidates[
    abs(mean_distance - min(mean_distance)) <= sqrt(.Machine$double.eps)
  ]
  best[order(sample_table$run_accession[best])][1L]
}, integer(1L))

site_table$bray_medoid_run_accession <-
  sample_table$run_accession[bray_medoid_indices]

site_bray_matrix <- bray_matrix[
  bray_medoid_indices,
  bray_medoid_indices,
  drop = FALSE
]
rownames(site_bray_matrix) <- site_table$site_id
colnames(site_bray_matrix) <- site_table$site_id

distance_matrices <- list(
  Aitchison = site_aitchison_matrix,
  `Bray-Curtis` = site_bray_matrix
)

data.table::fwrite(
  site_table,
  file.path(TABLE_DIR, "07_site_registry_and_Bray_medoids.csv")
)
data.table::fwrite(
  data.frame(
    n_samples = nrow(sample_table),
    n_sites = nrow(site_table),
    n_positive_axes = aitchison_geometry$n_positive_axes,
    min_eigenvalue = aitchison_geometry$min_eigenvalue,
    max_distance_reconstruction_error =
      aitchison_geometry$reconstruction_error
  ),
  file.path(TABLE_DIR, "07_Aitchison_geometry_reconstruction_audit.csv")
)

log_msg(
  "Aitchison geometry reconstruction error: ",
  sprintf("%.6g", aitchison_geometry$reconstruction_error)
)

## -----------------------------------------------------------------------------
## 4. Site-level factor designs
## -----------------------------------------------------------------------------

prepare_factor_design <- function(metadata, factor_name, fixed_levels = NULL) {
  if (!factor_name %in% names(metadata)) {
    return(list(status = "factor_missing", factor = factor_name))
  }

  factor_values <- clean_character(metadata[[factor_name]])
  df <- data.frame(
    site_index = metadata$site_index,
    site_id = metadata$site_id,
    study_id = metadata$study_id,
    n_studies = metadata$n_studies,
    factor_level = factor_values,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(!is.na(factor_level), !is.na(site_id))

  if (nrow(df) == 0L) {
    return(list(status = "no_complete_data", factor = factor_name))
  }

  site_counts <- df %>%
    dplyr::count(factor_level, name = "n_unique_sites") %>%
    dplyr::arrange(dplyr::desc(n_unique_sites), factor_level)

  if (is.null(fixed_levels)) {
    eligible_levels <- site_counts %>%
      dplyr::filter(n_unique_sites >= MIN_SITES_PER_LEVEL) %>%
      dplyr::pull(factor_level)
  } else {
    eligible_levels <- intersect(fixed_levels, site_counts$factor_level)
    eligible_levels <- eligible_levels[
      site_counts$n_unique_sites[match(eligible_levels, site_counts$factor_level)] >=
        MIN_SITES_PER_LEVEL
    ]
  }

  if (length(eligible_levels) < 2L) {
    return(list(
      status = "fewer_than_two_eligible_levels",
      factor = factor_name,
      site_counts = site_counts,
      eligible_levels = eligible_levels
    ))
  }

  eligible_counts <- site_counts %>%
    dplyr::filter(factor_level %in% eligible_levels)
  n_sites_per_level <- min(eligible_counts$n_unique_sites)

  list(
    status = "PASS",
    factor = factor_name,
    data = df %>% dplyr::filter(factor_level %in% eligible_levels),
    site_counts = site_counts,
    eligible_counts = eligible_counts,
    eligible_levels = sort(eligible_levels),
    n_sites_per_level = n_sites_per_level
  )
}

draw_balanced_indices <- function(design) {
  unlist(lapply(design$eligible_levels, function(level) {
    candidates <- design$data$site_index[
      design$data$factor_level == level
    ]
    sample(candidates, design$n_sites_per_level, replace = FALSE)
  }), use.names = FALSE)
}

test_one_distance <- function(
    distance_matrix,
    indices,
    group,
    permutations,
    permutation_seed) {

  group <- droplevels(as.factor(group))
  if (length(group) != length(indices)) stopf("group e indices no están alineados.")
  if (nlevels(group) < 2L) stopf("PERMANOVA requiere al menos dos niveles.")

  d <- as_dist_subset(distance_matrix, indices)

  set.seed(permutation_seed)
  ad <- vegan::adonis2(
    d ~ group,
    permutations = permutations,
    by = NULL
  )

  dispersion <- vegan::betadisper(
    d,
    group,
    type = "median",
    bias.adjust = TRUE
  )
  set.seed(permutation_seed + 1L)
  dispersion_test <- vegan::permutest(
    dispersion,
    permutations = permutations
  )

  data.frame(
    n_sites = length(indices),
    n_levels = nlevels(group),
    df = ad$Df[1L],
    F = ad$F[1L],
    R2 = ad$R2[1L],
    p_value = ad$`Pr(>F)`[1L],
    dispersion_F = dispersion_test$tab$F[1L],
    dispersion_p = dispersion_test$tab$`Pr(>F)`[1L]
  )
}

factors_to_test <- c("REALM", "PROVINCE", "Ecoregion", "Country")
primary_designs <- setNames(
  lapply(factors_to_test, function(f) prepare_factor_design(site_table, f)),
  factors_to_test
)

factor_level_counts <- dplyr::bind_rows(lapply(primary_designs, function(design) {
  if (is.null(design$site_counts)) return(NULL)
  design$site_counts %>%
    dplyr::mutate(
      factor = design$factor,
      eligible_primary = factor_level %in% design$eligible_levels,
      n_balanced_sites_per_level = if (design$status == "PASS") {
        design$n_sites_per_level
      } else {
        NA_integer_
      }
    )
})) %>%
  dplyr::select(
    factor,
    factor_level,
    n_unique_sites,
    eligible_primary,
    n_balanced_sites_per_level
  )

data.table::fwrite(
  factor_level_counts,
  file.path(TABLE_DIR, "07_factor_level_unique_site_counts.csv")
)

design_summary <- dplyr::bind_rows(lapply(primary_designs, function(design) {
  data.frame(
    factor = design$factor,
    status = design$status,
    n_eligible_levels = length(design$eligible_levels),
    eligible_levels = paste(design$eligible_levels, collapse = " | "),
    n_sites_per_level = if (design$status == "PASS") {
      design$n_sites_per_level
    } else {
      NA_integer_
    },
    n_sites_per_iteration = if (design$status == "PASS") {
      design$n_sites_per_level * length(design$eligible_levels)
    } else {
      NA_integer_
    },
    stringsAsFactors = FALSE
  )
}))
data.table::fwrite(design_summary, file.path(TABLE_DIR, "07_balanced_design_summary.csv"))

## -----------------------------------------------------------------------------
## 5. Primary PERMANOVA using all eligible independent sites
## -----------------------------------------------------------------------------

primary_results <- list()
primary_counter <- 0L

for (factor_index in seq_along(primary_designs)) {
  factor_name <- names(primary_designs)[factor_index]
  design <- primary_designs[[factor_index]]
  if (design$status != "PASS") next

  indices <- design$data$site_index
  group <- design$data$factor_level

  for (metric_index in seq_along(distance_matrices)) {
    metric_name <- names(distance_matrices)[metric_index]
    primary_counter <- primary_counter + 1L

    primary_results[[primary_counter]] <- test_one_distance(
      distance_matrix = distance_matrices[[metric_index]],
      indices = indices,
      group = group,
      permutations = PRIMARY_PERMUTATIONS,
      permutation_seed = SEED + factor_index * 10000L + metric_index
    ) %>%
      dplyr::mutate(
        factor = factor_name,
        metric = metric_name,
        eligible_levels = paste(design$eligible_levels, collapse = " | "),
        permutation_scheme = "free_among_independent_sites"
      )
  }
}

primary_results <- dplyr::bind_rows(primary_results) %>%
  dplyr::group_by(metric) %>%
  dplyr::mutate(
    p_BH = stats::p.adjust(p_value, method = "BH"),
    dispersion_p_BH = stats::p.adjust(dispersion_p, method = "BH")
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    factor,
    metric,
    n_sites,
    n_levels,
    eligible_levels,
    df,
    F,
    R2,
    p_value,
    p_BH,
    dispersion_F,
    dispersion_p,
    dispersion_p_BH,
    permutation_scheme
  )

if (nrow(primary_results) == 0L) {
  stopf("Ningún factor produjo una PERMANOVA primaria válida.")
}

data.table::fwrite(
  primary_results,
  file.path(TABLE_DIR, "07_primary_site_level_PERMANOVA.csv")
)

## -----------------------------------------------------------------------------
## 6. Repeated balanced sensitivity and multivariate dispersion
## -----------------------------------------------------------------------------

balanced_results <- list()
result_counter <- 0L

for (factor_index in seq_along(primary_designs)) {
  factor_name <- names(primary_designs)[factor_index]
  design <- primary_designs[[factor_index]]

  if (design$status != "PASS") {
    log_msg("Skipping ", factor_name, ": ", design$status)
    next
  }

  log_msg(
    factor_name,
    ": ",
    length(design$eligible_levels),
    " levels; ",
    design$n_sites_per_level,
    " sites per level; ",
    N_BALANCED_ITERATIONS,
    " iterations"
  )

  factor_results <- vector(
    "list",
    length = N_BALANCED_ITERATIONS * length(distance_matrices)
  )
  local_counter <- 0L

  for (iteration in seq_len(N_BALANCED_ITERATIONS)) {
    set.seed(SEED + factor_index * 100000L + iteration)
    indices <- draw_balanced_indices(design)
    group <- clean_character(site_table[[factor_name]][indices])

    for (metric_index in seq_along(distance_matrices)) {
      metric_name <- names(distance_matrices)[metric_index]
      local_counter <- local_counter + 1L

      test_result <- test_one_distance(
        distance_matrix = distance_matrices[[metric_index]],
        indices = indices,
        group = group,
        permutations = PERMUTATIONS_PER_ITERATION,
        permutation_seed = SEED + factor_index * 1000000L +
          iteration * 10L + metric_index
      )

      factor_results[[local_counter]] <- test_result %>%
        dplyr::mutate(
          factor = factor_name,
          metric = metric_name,
          iteration = iteration,
          n_sites_per_level = design$n_sites_per_level,
          eligible_levels = paste(design$eligible_levels, collapse = " | ")
        )
    }
  }

  result_counter <- result_counter + 1L
  balanced_results[[result_counter]] <- dplyr::bind_rows(factor_results)
}

balanced_results <- dplyr::bind_rows(balanced_results) %>%
  dplyr::select(
    factor,
    metric,
    iteration,
    n_sites,
    n_levels,
    n_sites_per_level,
    eligible_levels,
    df,
    F,
    R2,
    p_value,
    dispersion_F,
    dispersion_p
  )

if (nrow(balanced_results) == 0L) {
  stopf("Ningún factor produjo un diseño balanceado válido.")
}

data.table::fwrite(
  balanced_results,
  file.path(TABLE_DIR, "07_balanced_PERMANOVA_all_iterations.csv.gz"),
  compress = "gzip"
)

balanced_summary <- balanced_results %>%
  dplyr::group_by(factor, metric) %>%
  dplyr::summarise(
    n_iterations = dplyr::n(),
    n_levels = dplyr::first(n_levels),
    n_sites_per_level = dplyr::first(n_sites_per_level),
    n_sites_per_iteration = dplyr::first(n_sites),
    mean_R2 = mean(R2, na.rm = TRUE),
    sd_R2 = stats::sd(R2, na.rm = TRUE),
    median_R2 = stats::median(R2, na.rm = TRUE),
    q025_R2 = stats::quantile(R2, 0.025, na.rm = TRUE),
    q975_R2 = stats::quantile(R2, 0.975, na.rm = TRUE),
    mean_F = mean(F, na.rm = TRUE),
    median_F = stats::median(F, na.rm = TRUE),
    median_p = stats::median(p_value, na.rm = TRUE),
    proportion_PERMANOVA_p_le_0_05 = mean(p_value <= 0.05, na.rm = TRUE),
    median_dispersion_p = stats::median(dispersion_p, na.rm = TRUE),
    proportion_dispersion_p_le_0_05 = mean(dispersion_p <= 0.05, na.rm = TRUE),
    .groups = "drop"
  )

data.table::fwrite(
  balanced_summary,
  file.path(TABLE_DIR, "07_Supplementary_Table_balanced_PERMANOVA.csv")
)

## -----------------------------------------------------------------------------
## 7. Study-blocked Aitchison sensitivity when estimable
## -----------------------------------------------------------------------------

study_block_results <- list()
study_block_counter <- 0L

for (factor_name in names(primary_designs)) {
  design <- primary_designs[[factor_name]]
  if (design$status != "PASS") next

  block_data <- design$data %>%
    dplyr::filter(n_studies == 1L, !is.na(study_id))

  study_overlap <- block_data %>%
    dplyr::group_by(study_id) %>%
    dplyr::summarise(
      n_sites = dplyr::n(),
      n_factor_levels = dplyr::n_distinct(factor_level),
      .groups = "drop"
    )

  cross_level_studies <- study_overlap %>%
    dplyr::filter(n_factor_levels >= 2L)

  level_counts <- block_data %>%
    dplyr::count(factor_level, name = "n_sites")

  block_status <- "PASS"
  if (!setequal(level_counts$factor_level, design$eligible_levels) ||
      any(level_counts$n_sites < MIN_SITES_PER_LEVEL)) {
    block_status <- "NOT_ESTIMABLE_PRIMARY_LEVEL_LOST"
  } else if (nrow(cross_level_studies) < 2L) {
    block_status <- "NOT_ESTIMABLE_STUDY_NESTED_WITH_FACTOR"
  }

  study_block_counter <- study_block_counter + 1L

  if (block_status != "PASS") {
    study_block_results[[study_block_counter]] <- data.frame(
      factor = factor_name,
      status = block_status,
      n_sites = nrow(block_data),
      n_levels = dplyr::n_distinct(block_data$factor_level),
      n_studies = dplyr::n_distinct(block_data$study_id),
      n_cross_level_studies = nrow(cross_level_studies),
      F = NA_real_,
      R2 = NA_real_,
      p_value = NA_real_,
      dispersion_F = NA_real_,
      dispersion_p = NA_real_,
      stringsAsFactors = FALSE
    )
    next
  }

  permutation_control <- permute::how(
    nperm = STUDY_BLOCK_PERMUTATIONS,
    blocks = factor(block_data$study_id)
  )

  blocked_test <- tryCatch(
    test_one_distance(
      distance_matrix = site_aitchison_matrix,
      indices = block_data$site_index,
      group = block_data$factor_level,
      permutations = permutation_control,
      permutation_seed = string_seed("STUDY_BLOCK", factor_name)
    ),
    error = function(e) e
  )

  if (inherits(blocked_test, "error")) {
    study_block_results[[study_block_counter]] <- data.frame(
      factor = factor_name,
      status = paste0("ERROR: ", conditionMessage(blocked_test)),
      n_sites = nrow(block_data),
      n_levels = dplyr::n_distinct(block_data$factor_level),
      n_studies = dplyr::n_distinct(block_data$study_id),
      n_cross_level_studies = nrow(cross_level_studies),
      F = NA_real_,
      R2 = NA_real_,
      p_value = NA_real_,
      dispersion_F = NA_real_,
      dispersion_p = NA_real_,
      stringsAsFactors = FALSE
    )
  } else {
    study_block_results[[study_block_counter]] <- blocked_test %>%
      dplyr::mutate(
        factor = factor_name,
        status = "PASS",
        n_studies = dplyr::n_distinct(block_data$study_id),
        n_cross_level_studies = nrow(cross_level_studies)
      ) %>%
      dplyr::select(
        factor,
        status,
        n_sites,
        n_levels,
        n_studies,
        n_cross_level_studies,
        F,
        R2,
        p_value,
        dispersion_F,
        dispersion_p
      )
  }
}

study_block_results <- dplyr::bind_rows(study_block_results)
if (ncol(study_block_results) == 0L) {
  study_block_results <- data.frame(
    factor = character(),
    status = character(),
    n_sites = integer(),
    n_levels = integer(),
    n_studies = integer(),
    n_cross_level_studies = integer(),
    F = numeric(),
    R2 = numeric(),
    p_value = numeric(),
    dispersion_F = numeric(),
    dispersion_p = numeric()
  )
}
data.table::fwrite(
  study_block_results,
  file.path(TABLE_DIR, "07_study_blocked_Aitchison_PERMANOVA.csv")
)

## -----------------------------------------------------------------------------
## 8. Leave-one-study-out Aitchison sensitivity
## -----------------------------------------------------------------------------

study_values <- site_table %>%
  dplyr::filter(n_studies == 1L, !is.na(study_id)) %>%
  dplyr::pull(study_id) %>%
  unique() %>%
  sort()

loso_results <- list()
loso_counter <- 0L

for (factor_name in names(primary_designs)) {
  primary_design <- primary_designs[[factor_name]]
  if (primary_design$status != "PASS") next

  for (study in study_values) {
    sites_removed <- site_table$n_studies == 1L &
      !is.na(site_table$study_id) & site_table$study_id == study
    meta_excluded <- site_table[!sites_removed, , drop = FALSE]

    design_excluded <- prepare_factor_design(
      meta_excluded,
      factor_name,
      fixed_levels = primary_design$eligible_levels
    )

    loso_counter <- loso_counter + 1L
    if (design_excluded$status != "PASS" ||
        !setequal(design_excluded$eligible_levels, primary_design$eligible_levels)) {
      loso_results[[loso_counter]] <- data.frame(
        factor = factor_name,
        study_excluded = study,
        status = "SKIPPED_PRIMARY_LEVEL_LOST",
        n_sites_excluded = sum(sites_removed),
        n_sites = NA_integer_,
        n_levels = length(design_excluded$eligible_levels),
        F = NA_real_,
        R2 = NA_real_,
        p_value = NA_real_,
        dispersion_F = NA_real_,
        dispersion_p = NA_real_
      )
      next
    }

    test_result <- test_one_distance(
      distance_matrix = site_aitchison_matrix,
      indices = design_excluded$data$site_index,
      group = design_excluded$data$factor_level,
      permutations = LOSO_PERMUTATIONS,
      permutation_seed = string_seed("LOSO_PERM", factor_name, study)
    )

    loso_results[[loso_counter]] <- test_result %>%
      dplyr::mutate(
        factor = factor_name,
        study_excluded = study,
        status = "PASS",
        n_sites_excluded = sum(sites_removed)
      ) %>%
      dplyr::select(
        factor,
        study_excluded,
        status,
        n_sites_excluded,
        n_sites,
        n_levels,
        F,
        R2,
        p_value,
        dispersion_F,
        dispersion_p
      )
  }
}

loso_results <- dplyr::bind_rows(loso_results)
if (ncol(loso_results) == 0L) {
  loso_results <- data.frame(
    factor = character(),
    study_excluded = character(),
    status = character(),
    n_sites_excluded = integer(),
    n_sites = integer(),
    n_levels = integer(),
    F = numeric(),
    R2 = numeric(),
    p_value = numeric(),
    dispersion_F = numeric(),
    dispersion_p = numeric()
  )
}
data.table::fwrite(
  loso_results,
  file.path(TABLE_DIR, "07_leave_one_study_out_Aitchison.csv")
)

loso_summary <- loso_results %>%
  dplyr::filter(status == "PASS") %>%
  dplyr::group_by(factor) %>%
  dplyr::summarise(
    n_evaluable_studies = dplyr::n(),
    min_R2 = min(R2, na.rm = TRUE),
    median_R2 = stats::median(R2, na.rm = TRUE),
    max_R2 = max(R2, na.rm = TRUE),
    proportion_p_le_0_05 = mean(p_value <= 0.05, na.rm = TRUE),
    proportion_dispersion_p_le_0_05 = mean(dispersion_p <= 0.05, na.rm = TRUE),
    .groups = "drop"
  )

data.table::fwrite(
  loso_summary,
  file.path(TABLE_DIR, "07_leave_one_study_out_Aitchison_summary.csv")
)

## -----------------------------------------------------------------------------
## 9. Leave-one-country-out sensitivity for marine realm
## -----------------------------------------------------------------------------

loco_results <- list()
loco_counter <- 0L
realm_design <- primary_designs[["REALM"]]
country_values <- sort(unique_non_missing(site_table$Country))

if (realm_design$status == "PASS") {
  for (country in country_values) {
    sites_removed <- !is.na(site_table$Country) & site_table$Country == country
    meta_excluded <- site_table[!sites_removed, , drop = FALSE]
    design_excluded <- prepare_factor_design(
      meta_excluded,
      "REALM",
      fixed_levels = realm_design$eligible_levels
    )

    loco_counter <- loco_counter + 1L
    if (design_excluded$status != "PASS" ||
        !setequal(design_excluded$eligible_levels, realm_design$eligible_levels)) {
      loco_results[[loco_counter]] <- data.frame(
        country_excluded = country,
        status = "SKIPPED_PRIMARY_LEVEL_LOST",
        n_sites_excluded = sum(sites_removed),
        n_sites = NA_integer_,
        n_levels = length(design_excluded$eligible_levels),
        F = NA_real_,
        R2 = NA_real_,
        p_value = NA_real_,
        dispersion_F = NA_real_,
        dispersion_p = NA_real_
      )
      next
    }

    test_result <- test_one_distance(
      distance_matrix = site_aitchison_matrix,
      indices = design_excluded$data$site_index,
      group = design_excluded$data$factor_level,
      permutations = LOCO_PERMUTATIONS,
      permutation_seed = string_seed("LOCO_REALM", country)
    )

    loco_results[[loco_counter]] <- test_result %>%
      dplyr::mutate(
        country_excluded = country,
        status = "PASS",
        n_sites_excluded = sum(sites_removed)
      ) %>%
      dplyr::select(
        country_excluded,
        status,
        n_sites_excluded,
        n_sites,
        n_levels,
        F,
        R2,
        p_value,
        dispersion_F,
        dispersion_p
      )
  }
}

loco_results <- dplyr::bind_rows(loco_results)
if (ncol(loco_results) == 0L) {
  loco_results <- data.frame(
    country_excluded = character(),
    status = character(),
    n_sites_excluded = integer(),
    n_sites = integer(),
    n_levels = integer(),
    F = numeric(),
    R2 = numeric(),
    p_value = numeric(),
    dispersion_F = numeric(),
    dispersion_p = numeric()
  )
}
data.table::fwrite(
  loco_results,
  file.path(TABLE_DIR, "07_leave_one_country_out_REALM_Aitchison.csv")
)

## -----------------------------------------------------------------------------
## 10. Supplementary figure
## -----------------------------------------------------------------------------

factor_labels <- c(
  "REALM" = "Marine realm",
  "PROVINCE" = "Marine province",
  "Ecoregion" = "Marine ecoregion",
  "Country" = "Country"
)

plot_data <- balanced_results %>%
  dplyr::mutate(
    factor_label = factor(factor_labels[factor], levels = rev(factor_labels)),
    metric = factor(metric, levels = c("Aitchison", "Bray-Curtis"))
  )

p_supp <- ggplot(
  plot_data,
  aes(x = factor_label, y = R2, fill = metric)
) +
  geom_violin(
    alpha = 0.32,
    color = NA,
    scale = "width",
    position = position_dodge(width = 0.78)
  ) +
  geom_boxplot(
    width = 0.16,
    outlier.shape = NA,
    linewidth = 0.35,
    alpha = 0.78,
    position = position_dodge(width = 0.78)
  ) +
  coord_flip() +
  scale_y_continuous(
    labels = scales::label_percent(accuracy = 1),
    expand = expansion(mult = c(0, 0.04))
  ) +
  scale_fill_manual(
    values = c("Aitchison" = "#3B82A0", "Bray-Curtis" = "#D97706"),
    name = "Distance"
  ) +
  labs(
    x = NULL,
    y = expression(paste("Variance explained (", R^2, ")"))
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text = element_text(color = "black"),
    legend.position = "top",
    legend.title = element_text(face = "bold")
  )

ggsave(
  file.path(PLOT_DIR, "07_Supplementary_balanced_PERMANOVA.png"),
  p_supp,
  width = 8.0,
  height = 5.2,
  units = "in",
  dpi = 600,
  bg = "white"
)
ggsave(
  file.path(PLOT_DIR, "07_Supplementary_balanced_PERMANOVA.pdf"),
  p_supp,
  width = 8.0,
  height = 5.2,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)

## -----------------------------------------------------------------------------
## 11. Save objects and provenance
## -----------------------------------------------------------------------------

run_summary <- data.frame(
  status = "PASS",
  run_name = RUN_NAME,
  panel_d_object = PANEL_D_OBJECT_PATH,
  n_samples = nrow(sample_table),
  n_coordinate_sites = nrow(site_table),
  factors_requested = paste(factors_to_test, collapse = " | "),
  factors_evaluated = paste(unique(primary_results$factor), collapse = " | "),
  primary_metric = "Aitchison",
  sensitivity_metric = "Bray-Curtis",
  primary_unit = "coordinate_site_Aitchison_centroid",
  min_sites_per_level = MIN_SITES_PER_LEVEL,
  primary_permutations = PRIMARY_PERMUTATIONS,
  balanced_iterations = N_BALANCED_ITERATIONS,
  permutations_per_iteration = PERMUTATIONS_PER_ITERATION,
  study_block_permutations = STUDY_BLOCK_PERMUTATIONS,
  leave_one_study_out_permutations = LOSO_PERMUTATIONS,
  leave_one_country_out_permutations = LOCO_PERMUTATIONS,
  Aitchison_reconstruction_error = aitchison_geometry$reconstruction_error,
  seed = SEED,
  stringsAsFactors = FALSE
)
data.table::fwrite(run_summary, file.path(TABLE_DIR, "07_run_summary.csv"))

saveRDS(
  list(
    site_table = site_table,
    site_centroid_coordinates = site_centroid_coordinates,
    site_aitchison_distance = stats::as.dist(site_aitchison_matrix),
    site_bray_medoid_distance = stats::as.dist(site_bray_matrix),
    primary_designs = primary_designs,
    primary_results = primary_results,
    balanced_results = balanced_results,
    balanced_summary = balanced_summary,
    study_block_results = study_block_results,
    leave_one_study_out = loso_results,
    leave_one_study_out_summary = loso_summary,
    leave_one_country_out_REALM = loco_results,
    source_panel_d_object = PANEL_D_OBJECT_PATH,
    parameters = list(
      seed = SEED,
      min_sites_per_level = MIN_SITES_PER_LEVEL,
      primary_permutations = PRIMARY_PERMUTATIONS,
      balanced_iterations = N_BALANCED_ITERATIONS,
      permutations_per_iteration = PERMUTATIONS_PER_ITERATION,
      study_block_permutations = STUDY_BLOCK_PERMUTATIONS,
      loso_permutations = LOSO_PERMUTATIONS,
      loco_permutations = LOCO_PERMUTATIONS
    )
  ),
  file.path(RDS_DIR, "07_corrected_PERMANOVA_objects.rds")
)

writeLines(capture.output(sessionInfo()), file.path(LOG_DIR, "07_sessionInfo.txt"))

log_msg("Primary site-level PERMANOVA rows: ", nrow(primary_results))
log_msg("Balanced sensitivity rows: ", nrow(balanced_results))
log_msg("Study-block rows: ", nrow(study_block_results))
log_msg("LOSO rows: ", nrow(loso_results))
log_msg("LOCO rows: ", nrow(loco_results))
log_msg("Primary results: ", file.path(TABLE_DIR, "07_primary_site_level_PERMANOVA.csv"))
log_msg("=== SCRIPT COMPLETED SUCCESSFULLY ===")
