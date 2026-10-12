#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# Fit a descriptive Gaussian GAM for log1p-transformed Shannon entropy.
# Fourteen environmental smooths are fitted jointly with a locality random effect.
# Complete cases are retained, followed by a four-standard-deviation response rule.
# The fitted model, numerical checks, residual diagnostics, smooth statistics and
# partial-effect curves are exported. Convergence does not establish model adequacy.

options(stringsAsFactors = FALSE, scipen = 999, bitmapType = "cairo")

suppressPackageStartupMessages({
  library(mgcv)
  library(dplyr)
  library(readr)
})

# =============================================================================
# 1. Configuration and helpers
# =============================================================================

BASE_DIR <- Sys.getenv("BIOGEOGRAPHY_BASE", unset = ".")
OUTPUT_ROOT <- Sys.getenv("BIOGEOGRAPHY_RESULTS", unset = file.path(BASE_DIR, "results"))
ALPHA_DIR <- Sys.getenv("ALPHA_INPUT_DIR", unset = file.path(OUTPUT_ROOT, "alpha"))
METADATA_CSV <- Sys.getenv("MAGENTA_METADATA_CSV", unset = file.path(BASE_DIR, "data", "metadata_taxonomic.csv"))
SCRIPT_ID <- "SHANNON_GAM"
OUT_DIR <- Sys.getenv("GAM_OUTPUT_DIR", unset = file.path(OUTPUT_ROOT, "gam"))
for (subdir in c("plots", "models", "logs", "tables")) {
  dir.create(file.path(OUT_DIR, subdir), recursive = TRUE, showWarnings = FALSE)
}
OUT_DIR <- normalizePath(OUT_DIR, mustWork = TRUE)
LOG_FILE <- file.path(OUT_DIR, "logs", "SHANNON_GAM.log")
writeLines(character(), LOG_FILE)
LATEST_POINTER <- file.path(OUTPUT_ROOT, "LATEST_SHANNON_GAM.txt")
if (file.exists(LATEST_POINTER)) unlink(LATEST_POINTER)
SAMPLE_MANIFEST <- Sys.getenv("ALPHA_SAMPLE_MANIFEST", unset = file.path(BASE_DIR, "config", "samples", "taxonomic_accessions.txt"))
if (!file.exists(SAMPLE_MANIFEST)) stop("Sample manifest is missing: ", SAMPLE_MANIFEST, call. = FALSE)
EXPECTED_ACCESSIONS <- trimws(readLines(SAMPLE_MANIFEST, warn = FALSE))
EXPECTED_ACCESSIONS <- EXPECTED_ACCESSIONS[nzchar(EXPECTED_ACCESSIONS) & !startsWith(EXPECTED_ACCESSIONS, "#")]
if (!length(EXPECTED_ACCESSIONS) || anyDuplicated(EXPECTED_ACCESSIONS)) stop("Invalid sample manifest.", call. = FALSE)
EXPECTED_ALPHA_N <- length(EXPECTED_ACCESSIONS)
K_ENVIRONMENTAL <- 6L
GAMMA_VALUE <- 1.4
MIN_COMPLETE_CASES <- 150L
MIN_UNIQUE_VALUES <- 6L
EFFECT_GRID_N <- 240L
CONVERGENCE_GRADIENT_LIMIT <- 0.001

ENV_VARIABLES <- c(
  "human_footprint_2020",
  "bio1",
  "bio2",
  "bio3",
  "bio13",
  "bio15",
  "gbif_shannon_H_2km",
  "gbif_richness_species_2km",
  "soc_mean",
  "ph_mean",
  "sand_mean",
  "silt_mean",
  "clay_mean",
  "bdod_mean"
)

ENV_LABELS <- c(
  human_footprint_2020 = "Human footprint index (2020)",
  bio1 = "Annual mean temperature (BIO1)",
  bio2 = "Mean diurnal temperature range (BIO2)",
  bio3 = "Isothermality (BIO3)",
  bio13 = "Precipitation of wettest month (BIO13)",
  bio15 = "Precipitation seasonality (BIO15)",
  gbif_shannon_H_2km = "GBIF Shannon diversity (2 km)",
  gbif_richness_species_2km = "GBIF species richness (2 km)",
  soc_mean = "Soil organic carbon",
  ph_mean = "Soil pH",
  sand_mean = "Soil sand content",
  silt_mean = "Soil silt content",
  clay_mean = "Soil clay content",
  bdod_mean = "Soil bulk density"
)

msg <- function(...) {
  line <- paste0(
    sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0(..., collapse = "")
  )
  cat(line, "\n", sep = "")
  cat(line, "\n", file = LOG_FILE, append = TRUE, sep = "")
}

stop_with <- function(...) stop(paste0(..., collapse = ""), call. = FALSE)

normalize_text <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "null", "NULL")] <- NA_character_
  x
}

numeric_or_na <- function(x) suppressWarnings(as.numeric(x))

median_or_na <- function(x) {
  x <- numeric_or_na(x)
  if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
}

single_text_or_na <- function(x) {
  values <- unique(normalize_text(x))
  values <- values[!is.na(values)]
  if (length(values) == 0L) NA_character_ else values[1L]
}

msg("Starting full-variable Shannon GAM")
msg("Explicit alpha-diversity input directory: ", ALPHA_DIR)
msg("Environmental metadata: ", METADATA_CSV)
msg("Environmental predictors: ", length(ENV_VARIABLES))
msg("Output: ", OUT_DIR)

if (!nzchar(ALPHA_DIR) || !dir.exists(ALPHA_DIR)) {
  stop_with("Invalid explicit alpha-diversity input directory: ", ALPHA_DIR)
}

if (!file.exists(METADATA_CSV)) {
  stop_with("Required metadata file not found: ", METADATA_CSV)
}

# =============================================================================
# 2. Import the Figure 2a sample cohort
# =============================================================================

ALPHA_CSV <- file.path(
  ALPHA_DIR,
  "04a_alpha_diversity_filtered_samplewise.csv"
)

if (!file.exists(ALPHA_CSV)) {
  stop_with("Figure 2a samplewise table not found: ", ALPHA_CSV)
}

alpha_raw <- readr::read_csv(
  ALPHA_CSV,
  show_col_types = FALSE,
  na = c("", "NA", "N/A", "null", "NULL")
)

required_alpha <- c("run_accession", "Locality", "q1")
missing_alpha <- setdiff(required_alpha, names(alpha_raw))

if (length(missing_alpha) > 0L) {
  stop_with(
    "Missing required columns in the Figure 2a table: ",
    paste(missing_alpha, collapse = ", ")
  )
}

alpha <- alpha_raw %>%
  transmute(
    run_accession = normalize_text(run_accession),
    Locality_alpha = normalize_text(Locality),
    q1 = numeric_or_na(q1)
  )

if (nrow(alpha) != EXPECTED_ALPHA_N) {
  stop_with(
    "Figure 2a cohort contains ", nrow(alpha),
    " rows; expected ", EXPECTED_ALPHA_N, "."
  )
}

if (!setequal(alpha$run_accession, EXPECTED_ACCESSIONS)) stop_with("Alpha table does not match the sample manifest.")

if (anyNA(alpha$run_accession) || anyDuplicated(alpha$run_accession)) {
  stop_with("Figure 2a run_accession values are missing or duplicated.")
}

if (any(!is.finite(alpha$q1)) || any(alpha$q1 < 1)) {
  stop_with("Hill q = 1 must contain finite values greater than or equal to 1.")
}

alpha <- alpha %>%
  mutate(
    shannon = log(q1),
    y_log = log1p(shannon)
  )

msg("Inherited Figure 2a cohort: n = ", nrow(alpha))

# =============================================================================
# 3. Import and deterministically collapse duplicated metadata records
# =============================================================================

metadata_raw <- readr::read_csv(
  METADATA_CSV,
  show_col_types = FALSE,
  na = c("", "NA", "N/A", "null", "NULL")
)

names(metadata_raw) <- trimws(names(metadata_raw))

if ("phh2o_mean" %in% names(metadata_raw) &&
    !"ph_mean" %in% names(metadata_raw)) {
  metadata_raw <- metadata_raw %>% rename(ph_mean = phh2o_mean)
}

required_metadata <- c("run_accession", "Locality", ENV_VARIABLES)
missing_metadata <- setdiff(required_metadata, names(metadata_raw))

if (length(missing_metadata) > 0L) {
  stop_with(
    "Missing required columns in environmental metadata: ",
    paste(missing_metadata, collapse = ", ")
  )
}

metadata_selected <- metadata_raw %>%
  transmute(
    run_accession = normalize_text(run_accession),
    Locality_meta = normalize_text(Locality),
    across(all_of(ENV_VARIABLES), numeric_or_na)
  ) %>%
  filter(!is.na(run_accession))

duplicate_ids <- metadata_selected %>%
  count(run_accession, name = "n_rows") %>%
  filter(n_rows > 1L) %>%
  pull(run_accession)

duplicate_records <- metadata_selected %>%
  filter(run_accession %in% duplicate_ids) %>%
  arrange(run_accession)

write.csv(
  duplicate_records,
  file.path(OUT_DIR, "tables", "17A_duplicate_metadata_records_audit.csv"),
  row.names = FALSE
)

audit_variables <- c("Locality_meta", ENV_VARIABLES)

if (length(duplicate_ids) == 0L) {
  duplicate_conflict_audit <- data.frame(
    run_accession = character(),
    variable = character(),
    n_distinct_nonmissing = integer(),
    observed_values = character(),
    stringsAsFactors = FALSE
  )
} else {
  duplicate_conflict_audit <- bind_rows(lapply(duplicate_ids, function(id) {
    current <- metadata_selected %>% filter(run_accession == id)

    bind_rows(lapply(audit_variables, function(variable) {
      values <- current[[variable]]
      values <- values[!is.na(values)]
      distinct_values <- unique(values)

      data.frame(
        run_accession = id,
        variable = variable,
        n_distinct_nonmissing = length(distinct_values),
        observed_values = paste(distinct_values, collapse = " | "),
        stringsAsFactors = FALSE
      )
    }))
  }))
}

write.csv(
  duplicate_conflict_audit,
  file.path(OUT_DIR, "tables", "17A_duplicate_metadata_conflict_audit.csv"),
  row.names = FALSE
)

locality_duplicate_conflicts <- duplicate_conflict_audit %>%
  filter(variable == "Locality_meta", n_distinct_nonmissing > 1L)

if (nrow(locality_duplicate_conflicts) > 0L) {
  stop_with(
    "Duplicated run_accession records contain conflicting Locality values. ",
    "See 17A_duplicate_metadata_conflict_audit.csv."
  )
}

numeric_duplicate_conflicts <- duplicate_conflict_audit %>%
  filter(variable %in% ENV_VARIABLES, n_distinct_nonmissing > 1L)

if (nrow(numeric_duplicate_conflicts) > 0L) {
  msg(
    "WARNING: ", nrow(numeric_duplicate_conflicts),
    " duplicated accession-variable combinations contain conflicting numeric ",
    "values. Their within-accession medians will be used; see the audit table."
  )
}

metadata <- metadata_selected %>%
  group_by(run_accession) %>%
  summarise(
    Locality_meta = single_text_or_na(Locality_meta),
    across(all_of(ENV_VARIABLES), median_or_na),
    .groups = "drop"
  )

if (anyDuplicated(metadata$run_accession)) {
  stop_with("Metadata collapsing did not yield one row per run_accession.")
}

# =============================================================================
# 4. Build the full-model analysis table
# =============================================================================

datos_raw <- alpha %>%
  left_join(metadata, by = "run_accession") %>%
  mutate(
    Locality = Locality_meta
  )

if (nrow(datos_raw) != nrow(alpha) ||
    n_distinct(datos_raw$run_accession) != nrow(alpha)) {
  stop_with("The metadata join changed the size or uniqueness of the cohort.")
}

locality_source_conflicts <- datos_raw %>%
  filter(
    !is.na(Locality_alpha),
    !is.na(Locality_meta),
    Locality_alpha != Locality_meta
  ) %>%
  select(run_accession, Locality_alpha, Locality_meta)

write.csv(
  locality_source_conflicts,
  file.path(OUT_DIR, "tables", "17A_locality_source_conflicts.csv"),
  row.names = FALSE
)

if (nrow(locality_source_conflicts) > 0L) {
  msg(
    "WARNING: ", nrow(locality_source_conflicts),
    " locality labels differ between Figure 2a and the metadata CSV. ",
    "The external metadata locality label is retained."
  )
}

has_locality <- !is.na(datos_raw$Locality) & nzchar(datos_raw$Locality)

finite_numeric <- Reduce(
  `&`,
  lapply(
    c("y_log", ENV_VARIABLES),
    function(variable) is.finite(datos_raw[[variable]])
  )
)

complete_case_mask <- has_locality & finite_numeric

complete_case_audit <- datos_raw %>%
  mutate(
    complete_full_model_case = complete_case_mask,
    missing_or_nonfinite_variables = vapply(
      seq_len(n()),
      function(i) {
        problem_variables <- c("y_log", ENV_VARIABLES)[
          !vapply(
            c("y_log", ENV_VARIABLES),
            function(variable) is.finite(datos_raw[[variable]][i]),
            logical(1)
          )
        ]

        if (!has_locality[i]) {
          problem_variables <- c(problem_variables, "Locality")
        }

        paste(problem_variables, collapse = ";")
      },
      character(1)
    )
  )

write.csv(
  complete_case_audit,
  file.path(OUT_DIR, "tables", "17A_complete_case_audit.csv"),
  row.names = FALSE
)

datos_complete <- datos_raw[complete_case_mask, , drop = FALSE]

if (nrow(datos_complete) < MIN_COMPLETE_CASES) {
  stop_with(
    "Only ", nrow(datos_complete),
    " complete cases remain for the full model; at least ",
    MIN_COMPLETE_CASES, " are required."
  )
}

y_mean <- mean(datos_complete$y_log)
y_sd <- sd(datos_complete$y_log)

if (!is.finite(y_sd) || y_sd <= 0) {
  stop_with("The standard deviation of y_log is not positive.")
}

datos_complete <- datos_complete %>%
  mutate(
    response_outlier_4sd = abs(y_log - y_mean) >= 4 * y_sd
  )

write.csv(
  datos_complete %>% filter(response_outlier_4sd),
  file.path(OUT_DIR, "tables", "17A_response_4SD_exclusions.csv"),
  row.names = FALSE
)

datos <- datos_complete %>%
  filter(!response_outlier_4sd) %>%
  transmute(
    run_accession,
    q1,
    shannon,
    y_log,
    Locality = droplevels(factor(Locality)),
    across(all_of(ENV_VARIABLES))
  )

unique_values <- vapply(
  datos[ENV_VARIABLES],
  function(x) length(unique(x)),
  integer(1)
)

if (any(unique_values < MIN_UNIQUE_VALUES)) {
  stop_with(
    "Insufficient unique values for full-model smooths: ",
    paste(
      names(unique_values)[unique_values < MIN_UNIQUE_VALUES],
      collapse = ", "
    )
  )
}

MODEL_INPUT_CSV <- file.path(
  OUT_DIR,
  "tables",
  "17A_full_model_input.csv"
)

write.csv(datos, MODEL_INPUT_CSV, row.names = FALSE)

msg("Complete cases for all 14 predictors: n = ", nrow(datos_complete))
msg("Retained after 4-SD response rule: n = ", nrow(datos))
msg("Locality levels: ", nlevels(datos$Locality))
msg("Response retained: y_log = log1p(Shannon entropy)")

# =============================================================================
# 5. Fit the full-variable Gaussian mean GAM
# =============================================================================

environmental_smooths <- sprintf(
  "s(%s, k = %d, bs = 'tp')",
  ENV_VARIABLES,
  K_ENVIRONMENTAL
)

formula_text <- paste(
  "y_log ~",
  paste(
    c(environmental_smooths, "s(Locality, bs = 're')"),
    collapse = " + "
  )
)

model_formula <- as.formula(formula_text)

writeLines(
  formula_text,
  file.path(OUT_DIR, "logs", "17A_model_formula.txt")
)

msg("Fitting full Gaussian GAM with REML and outer/BFGS optimization")

set.seed(42)

m_full <- mgcv::gam(
  formula = model_formula,
  family = gaussian(link = "identity"),
  method = "REML",
  optimizer = c("outer", "bfgs"),
  select = TRUE,
  gamma = GAMMA_VALUE,
  control = mgcv::gam.control(maxit = 1000, trace = FALSE),
  na.action = na.fail,
  data = datos
)

MODEL_RDS <- file.path(
  OUT_DIR,
  "models",
  "17A_full_environmental_GAM_Shannon.rds"
)

saveRDS(m_full, MODEL_RDS)

summ <- summary(m_full)

outer_conv <- if (is.null(m_full$outer.info$conv)) {
  NA_character_
} else {
  as.character(m_full$outer.info$conv)[1L]
}

max_gradient <- if (is.null(m_full$outer.info$grad)) {
  NA_real_
} else {
  max(abs(m_full$outer.info$grad), na.rm = TRUE)
}

min_hessian_eigenvalue <- NA_real_

if (!is.null(m_full$outer.info$hess) &&
    is.matrix(m_full$outer.info$hess) &&
    all(is.finite(m_full$outer.info$hess))) {
  min_hessian_eigenvalue <- min(
    eigen(
      m_full$outer.info$hess,
      symmetric = TRUE,
      only.values = TRUE
    )$values
  )
}

hessian_ok <-
  is.finite(min_hessian_eigenvalue) &&
  min_hessian_eigenvalue > -0.000001

convergence_ok <-
  is.character(outer_conv) &&
  grepl("full convergence", outer_conv, ignore.case = TRUE) &&
  is.finite(max_gradient) &&
  max_gradient < CONVERGENCE_GRADIENT_LIMIT &&
  hessian_ok

model_metrics <- data.frame(
  n_inherited_cohort = nrow(alpha),
  n_complete_all_predictors = nrow(datos_complete),
  n_after_4sd_rule = nrow(datos),
  n_localities = nlevels(datos$Locality),
  n_environmental_predictors = length(ENV_VARIABLES),
  response = "log1p(Shannon entropy)",
  family = "Gaussian identity",
  estimation = "REML",
  optimizer = "outer/bfgs",
  environmental_basis_k = K_ENVIRONMENTAL,
  select = TRUE,
  gamma = GAMMA_VALUE,
  deviance_explained = summ$dev.expl,
  adjusted_r_squared = summ$r.sq,
  AIC = AIC(m_full),
  convergence_message = outer_conv,
  maximum_absolute_gradient = max_gradient,
  minimum_hessian_eigenvalue = min_hessian_eigenvalue,
  convergence_accepted = convergence_ok,
  stringsAsFactors = FALSE
)

write.csv(
  model_metrics,
  file.path(OUT_DIR, "tables", "17A_model_metrics.csv"),
  row.names = FALSE
)

capture.output(
  summ,
  file = file.path(OUT_DIR, "logs", "17A_model_summary.txt")
)

capture.output(
  m_full$outer.info,
  file = file.path(OUT_DIR, "logs", "17A_outer_optimizer_info.txt")
)

msg("Convergence message: ", outer_conv)
msg("Maximum absolute gradient: ", format(max_gradient, digits = 10))
msg("Minimum Hessian eigenvalue: ", format(min_hessian_eigenvalue, digits = 10))
msg("Deviance explained: ", round(100 * summ$dev.expl, 2), "%")

if (!convergence_ok) {
  stop_with(
    "The full environmental GAM did not meet convergence criteria. ",
    "No effect ranking was produced. See 17A_model_metrics.csv and ",
    "17A_outer_optimizer_info.txt."
  )
}

# =============================================================================
# 6. Export model statistics and diagnostics
# =============================================================================

smooth_table <- data.frame(
  term = rownames(summ$s.table),
  as.data.frame(summ$s.table, check.names = FALSE),
  row.names = NULL,
  check.names = FALSE
)

write.csv(
  smooth_table,
  file.path(OUT_DIR, "tables", "17A_all_smooth_terms.csv"),
  row.names = FALSE
)

parametric_table <- data.frame(
  term = rownames(summ$p.table),
  as.data.frame(summ$p.table, check.names = FALSE),
  row.names = NULL,
  check.names = FALSE
)

write.csv(
  parametric_table,
  file.path(OUT_DIR, "tables", "17A_parametric_terms.csv"),
  row.names = FALSE
)

set.seed(42)
k_check <- mgcv::k.check(m_full, n.rep = 400)

k_check_table <- data.frame(
  term = rownames(k_check),
  as.data.frame(k_check, check.names = FALSE),
  row.names = NULL,
  check.names = FALSE
)

write.csv(
  k_check_table,
  file.path(OUT_DIR, "tables", "17A_k_check.csv"),
  row.names = FALSE
)

concurvity_result <- mgcv::concurvity(m_full, full = TRUE)

concurvity_table <- data.frame(
  metric = rownames(concurvity_result),
  as.data.frame(concurvity_result, check.names = FALSE),
  row.names = NULL,
  check.names = FALSE
)

write.csv(
  concurvity_table,
  file.path(OUT_DIR, "tables", "17A_concurvity.csv"),
  row.names = FALSE
)

tryCatch({
  png(
    file.path(OUT_DIR, "plots", "17A_gam_check.png"),
    width = 1800,
    height = 1400,
    res = 180
  )
  mgcv::gam.check(m_full)
  dev.off()
}, error = function(error) {
  try(dev.off(), silent = TRUE)
  msg("WARNING: gam.check plot failed: ", conditionMessage(error))
})

capture.output(
  sessionInfo(),
  file = file.path(OUT_DIR, "logs", "17A_sessionInfo.txt")
)

# =============================================================================
# 7. Derive model-based partial-effect curves and rank all 14 effects
# =============================================================================

statistic_column <- intersect(c("F", "Chi.sq"), names(smooth_table))
p_column <- grep(
  "^p[-._ ]?value$|^p[-._ ]?val$",
  names(smooth_table),
  ignore.case = TRUE,
  value = TRUE
)

if (length(statistic_column) != 1L || length(p_column) != 1L) {
  stop_with(
    "Could not identify the smooth test statistic or P-value column."
  )
}

reference_row <- datos[1L, , drop = FALSE]

for (variable in ENV_VARIABLES) {
  reference_row[[variable]] <- median(datos[[variable]], na.rm = TRUE)
}

reference_row$Locality <- factor(
  levels(datos$Locality)[1L],
  levels = levels(datos$Locality)
)

extract_effect_curve <- function(variable) {
  observed <- datos[[variable]]
  display_range <- as.numeric(
    quantile(
      observed,
      probs = c(0.025, 0.975),
      na.rm = TRUE,
      names = FALSE,
      type = 7
    )
  )

  if (!all(is.finite(display_range)) || display_range[1L] >= display_range[2L]) {
    stop_with("Invalid central-95% range for ", variable, ".")
  }

  newdata <- reference_row[rep(1L, EFFECT_GRID_N), , drop = FALSE]
  newdata[[variable]] <- seq(
    display_range[1L],
    display_range[2L],
    length.out = EFFECT_GRID_N
  )

  prediction <- predict(
    m_full,
    newdata = newdata,
    type = "terms",
    se.fit = TRUE,
    exclude = "s(Locality)"
  )

  fit_matrix <- as.matrix(prediction$fit)
  se_matrix <- as.matrix(prediction$se.fit)
  term_name <- paste0("s(", variable, ")")
  term_index <- which(colnames(fit_matrix) == term_name)

  if (length(term_index) != 1L) {
    stop_with("Could not identify prediction term ", term_name, ".")
  }

  effect <- fit_matrix[, term_index]
  standard_error <- se_matrix[, term_index]

  data.frame(
    variable = variable,
    panel_label = unname(ENV_LABELS[variable]),
    x = newdata[[variable]],
    effect = effect,
    standard_error = standard_error,
    lower = effect - 1.96 * standard_error,
    upper = effect + 1.96 * standard_error,
    display_min = display_range[1L],
    display_max = display_range[2L],
    stringsAsFactors = FALSE
  )
}

effect_curves <- bind_rows(lapply(ENV_VARIABLES, extract_effect_curve))

smooth_statistics_for <- function(variable) {
  term_name <- paste0("s(", variable, ")")
  row <- smooth_table[smooth_table$term == term_name, , drop = FALSE]

  if (nrow(row) != 1L) {
    stop_with("Missing smooth statistics for ", term_name, ".")
  }

  data.frame(
    variable = variable,
    edf = as.numeric(row$edf[1L]),
    reference_df = as.numeric(row$Ref.df[1L]),
    test = statistic_column,
    statistic = as.numeric(row[[statistic_column]][1L]),
    p_value = as.numeric(row[[p_column]][1L]),
    stringsAsFactors = FALSE
  )
}

smooth_statistics <- bind_rows(lapply(ENV_VARIABLES, smooth_statistics_for))

effect_ranking <- effect_curves %>%
  group_by(variable, panel_label) %>%
  summarise(
    predictor_central_95_min = first(display_min),
    predictor_central_95_max = first(display_max),
    partial_effect_min = min(effect, na.rm = TRUE),
    partial_effect_max = max(effect, na.rm = TRUE),
    central_95_effect_amplitude = partial_effect_max - partial_effect_min,
    maximum_absolute_effect = max(abs(effect), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(smooth_statistics, by = "variable") %>%
  mutate(
    n_unique_predictor_values = unname(unique_values[variable])
  ) %>%
  arrange(desc(central_95_effect_amplitude), variable) %>%
  mutate(
    effect_rank = row_number(),
    selected_top4 = effect_rank <= 4L
  )

if (nrow(effect_ranking) != length(ENV_VARIABLES) ||
    sum(effect_ranking$selected_top4) != 4L) {
  stop_with("Effect ranking did not return 14 terms and exactly four selections.")
}

effect_curves <- effect_curves %>%
  left_join(
    effect_ranking %>%
      select(variable, effect_rank, selected_top4),
    by = "variable"
  )

RANKING_CSV <- file.path(
  OUT_DIR,
  "tables",
  "17A_environmental_effect_ranking.csv"
)

CURVES_CSV <- file.path(
  OUT_DIR,
  "tables",
  "17A_all_environmental_partial_effect_curves.csv"
)

write.csv(effect_ranking, RANKING_CSV, row.names = FALSE)
write.csv(effect_curves, CURVES_CSV, row.names = FALSE)

write.csv(
  effect_ranking %>% filter(selected_top4),
  file.path(OUT_DIR, "tables", "17A_selected_top4_effects.csv"),
  row.names = FALSE
)

writeLines(
  effect_ranking %>%
    filter(selected_top4) %>%
    arrange(effect_rank) %>%
    pull(variable),
  file.path(OUT_DIR, "17A_selected_top4_variables.txt")
)

msg(
  "Selected top four effects: ",
  paste(
    effect_ranking %>%
      filter(selected_top4) %>%
      arrange(effect_rank) %>%
      pull(variable),
    collapse = ", "
  )
)

# =============================================================================
# 8. Provenance, summary and final pointer
# =============================================================================

input_md5 <- data.frame(
  input = c("Figure_2a_samplewise", "Metadata_2026_copia.csv"),
  path = c(ALPHA_CSV, METADATA_CSV),
  md5 = unname(tools::md5sum(c(ALPHA_CSV, METADATA_CSV))),
  stringsAsFactors = FALSE
)

write.csv(
  input_md5,
  file.path(OUT_DIR, "tables", "17A_input_md5.csv"),
  row.names = FALSE
)

top4_summary <- effect_ranking %>%
  filter(selected_top4) %>%
  arrange(effect_rank) %>%
  transmute(
    text = paste0(
      effect_rank, ". ", variable,
      " (amplitude = ",
      format(central_95_effect_amplitude, digits = 6),
      ")"
    )
  ) %>%
  pull(text)

run_summary <- c(
  "Shannon diversity environmental GAM",
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("Figure 2a source: ", ALPHA_CSV),
  paste0("Metadata source: ", METADATA_CSV),
  "Join key: run_accession",
  paste0("Inherited Figure 2a samples: ", nrow(alpha)),
  paste0("Complete cases for all 14 predictors: ", nrow(datos_complete)),
  paste0(
    "Samples excluded by 4-SD response rule: ",
    sum(datos_complete$response_outlier_4sd)
  ),
  paste0("Analysis samples: ", nrow(datos)),
  paste0("Localities: ", nlevels(datos$Locality)),
  "Shannon entropy: log(Hill q = 1)",
  "Response: y_log = log1p(Shannon entropy)",
  paste0("Environmental predictors: ", paste(ENV_VARIABLES, collapse = ", ")),
  "Random effect: Locality",
  "Family: Gaussian identity GAM",
  "Estimation: REML with outer/BFGS optimization",
  "Smooth selection: double-penalty selection (select = TRUE)",
  paste0("Environmental basis dimension: k = ", K_ENVIRONMENTAL),
  paste0("Gamma: ", GAMMA_VALUE),
  paste0("Convergence message: ", outer_conv),
  paste0("Maximum absolute gradient: ", format(max_gradient, digits = 12)),
  paste0(
    "Minimum Hessian eigenvalue: ",
    format(min_hessian_eigenvalue, digits = 12)
  ),
  paste0("Deviance explained (%): ", round(100 * summ$dev.expl, 4)),
  paste0("Adjusted R-squared: ", round(summ$r.sq, 6)),
  paste0("AIC: ", format(AIC(m_full), digits = 12)),
  paste0(
    "Numeric duplicate conflicts collapsed by within-accession median: ",
    nrow(numeric_duplicate_conflicts)
  ),
  paste0(
    "Effect ranking criterion: partial-effect range across the central 95% ",
    "of each predictor distribution"
  ),
  "Selected top four effects:",
  top4_summary,
  paste0("Model RDS: ", MODEL_RDS),
  paste0("Ranking table: ", RANKING_CSV),
  paste0("All effect curves: ", CURVES_CSV),
  paste0("Output directory: ", OUT_DIR)
)

writeLines(run_summary, file.path(OUT_DIR, "17A_run_summary.txt"))
writeLines(OUT_DIR, LATEST_POINTER)

msg("Completed successfully: ", OUT_DIR)

