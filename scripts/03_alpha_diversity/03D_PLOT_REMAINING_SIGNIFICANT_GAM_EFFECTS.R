#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# Supplementary figure — plot-only stage.
# Reads the converged 14-predictor GAM produced by script 17A and plots every
# significant environmental smooth (P < 0.05) that was not selected among the
# four largest effects shown in Figure 2b. This script DOES NOT fit or refit a
# GAM.

options(stringsAsFactors = FALSE, scipen = 999, bitmapType = "cairo")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
})

BASE_DIR <- Sys.getenv("BIOGEOGRAPHY_BASE", unset = "/home/fjbalvino/Biogeografia_2026")
RESULT_DIR <- file.path(BASE_DIR, "resultados")

FULL_MODEL_POINTER <- file.path(
  RESULT_DIR,
  "LATEST_17A_GAM_SHANNON_ALL_VARIABLES_STABLE.txt"
)

SIGNIFICANCE_LEVEL <- 0.05
SCRIPT_ID <- "17C_PLOT_REMAINING_SIGNIFICANT_FROM_FULL_GAM"
RUN_STAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUT_DIR <- file.path(RESULT_DIR, paste0(SCRIPT_ID, "_", RUN_STAMP))

dir.create(file.path(OUT_DIR, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "logs"), recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(OUT_DIR, "logs", paste0(SCRIPT_ID, ".log"))
LATEST_POINTER <- file.path(RESULT_DIR, paste0("LATEST_", SCRIPT_ID, ".txt"))

writeLines(character(), LOG_FILE)
if (file.exists(LATEST_POINTER)) unlink(LATEST_POINTER)

msg <- function(...) {
  line <- paste0(
    sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0(..., collapse = "")
  )
  cat(line, "\n", sep = "")
  cat(line, "\n", file = LOG_FILE, append = TRUE, sep = "")
}

stop_with <- function(...) stop(paste0(..., collapse = ""), call. = FALSE)

as_flag <- function(x) {
  if (is.logical(x)) return(replace(x, is.na(x), FALSE))
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes")
}

format_p <- function(p) {
  p <- suppressWarnings(as.numeric(p))[1L]
  if (!is.finite(p)) return("P = NA")
  if (p < 0.001) return("P < 0.001")
  paste0("P = ", formatC(p, format = "f", digits = 3))
}

msg(
  "Starting supplementary plot-only workflow for the remaining significant ",
  "effects from the full 14-predictor GAM"
)
msg("Full-model pointer: ", FULL_MODEL_POINTER)
msg("Significance threshold: P < ", SIGNIFICANCE_LEVEL)
msg("Output: ", OUT_DIR)

if (!file.exists(FULL_MODEL_POINTER)) {
  stop_with(
    "Missing full-model pointer: ", FULL_MODEL_POINTER,
    ". Run 17A_GAM_SHANNON_ALL_VARIABLES_STABLE.R first."
  )
}

FULL_MODEL_DIR <- trimws(
  readLines(FULL_MODEL_POINTER, n = 1L, warn = FALSE)
)

if (!nzchar(FULL_MODEL_DIR) || !dir.exists(FULL_MODEL_DIR)) {
  stop_with("Invalid full-model directory: ", FULL_MODEL_DIR)
}

RANKING_CSV <- file.path(
  FULL_MODEL_DIR,
  "tables",
  "17A_environmental_effect_ranking.csv"
)

CURVES_CSV <- file.path(
  FULL_MODEL_DIR,
  "tables",
  "17A_all_environmental_partial_effect_curves.csv"
)

MODEL_INPUT_CSV <- file.path(
  FULL_MODEL_DIR,
  "tables",
  "17A_full_model_input.csv"
)

MODEL_METRICS_CSV <- file.path(
  FULL_MODEL_DIR,
  "tables",
  "17A_model_metrics.csv"
)

for (required_file in c(
  RANKING_CSV,
  CURVES_CSV,
  MODEL_INPUT_CSV,
  MODEL_METRICS_CSV
)) {
  if (!file.exists(required_file)) {
    stop_with("Missing required 17A output: ", required_file)
  }
}

ranking <- readr::read_csv(RANKING_CSV, show_col_types = FALSE)
curves <- readr::read_csv(CURVES_CSV, show_col_types = FALSE)
model_input <- readr::read_csv(MODEL_INPUT_CSV, show_col_types = FALSE)
metrics <- readr::read_csv(MODEL_METRICS_CSV, show_col_types = FALSE)

required_ranking <- c(
  "variable", "panel_label", "central_95_effect_amplitude",
  "edf", "reference_df", "statistic", "p_value", "effect_rank",
  "selected_top4"
)

required_curves <- c(
  "variable", "x", "effect", "lower", "upper",
  "display_min", "display_max", "effect_rank"
)

missing_ranking <- setdiff(required_ranking, names(ranking))
missing_curves <- setdiff(required_curves, names(curves))

if (length(missing_ranking) > 0L) {
  stop_with(
    "Ranking table lacks columns: ",
    paste(missing_ranking, collapse = ", ")
  )
}

if (length(missing_curves) > 0L) {
  stop_with(
    "Curve table lacks columns: ",
    paste(missing_curves, collapse = ", ")
  )
}

if (nrow(ranking) != 14L) {
  stop_with("Expected 14 ranked environmental effects; found ", nrow(ranking), ".")
}

if (nrow(metrics) != 1L ||
    !"convergence_accepted" %in% names(metrics) ||
    !isTRUE(as_flag(metrics$convergence_accepted)[1L])) {
  stop_with("The source full model is not marked as converged.")
}

ranking <- ranking %>%
  mutate(
    variable = as.character(variable),
    effect_rank = as.integer(effect_rank),
    p_value = suppressWarnings(as.numeric(p_value)),
    selected_top4 = as_flag(selected_top4)
  ) %>%
  arrange(effect_rank)

top4_variables <- ranking %>%
  filter(selected_top4) %>%
  arrange(effect_rank) %>%
  pull(variable)

if (length(top4_variables) != 4L || anyDuplicated(top4_variables)) {
  stop_with("The 17A ranking does not contain exactly four unique top effects.")
}

remaining_significant <- ranking %>%
  filter(
    !selected_top4,
    is.finite(p_value),
    p_value < SIGNIFICANCE_LEVEL
  ) %>%
  arrange(effect_rank)

if (nrow(remaining_significant) == 0L) {
  stop_with(
    "No significant environmental effects remain after excluding the top four."
  )
}

remaining_variables <- remaining_significant$variable

if (any(remaining_variables %in% top4_variables)) {
  stop_with("A selected supplementary term overlaps with the Figure 2b top four.")
}

missing_predictors <- setdiff(remaining_variables, names(model_input))

if (length(missing_predictors) > 0L) {
  stop_with(
    "Selected variables are absent from model input: ",
    paste(missing_predictors, collapse = ", ")
  )
}

effect_data <- curves %>%
  filter(variable %in% remaining_variables) %>%
  mutate(
    variable = factor(variable, levels = remaining_variables),
    effect_rank = as.integer(effect_rank)
  ) %>%
  arrange(effect_rank, x)

if (n_distinct(effect_data$variable) != length(remaining_variables)) {
  stop_with("At least one selected partial-effect curve is missing.")
}

msg(
  "Figure 2b top four excluded: ",
  paste(top4_variables, collapse = ", ")
)
msg(
  "Remaining significant effects selected: ",
  paste(remaining_variables, collapse = ", ")
)

write.csv(
  ranking,
  file.path(OUT_DIR, "tables", "17C_complete_effect_ranking.csv"),
  row.names = FALSE
)

write.csv(
  remaining_significant,
  file.path(
    OUT_DIR,
    "tables",
    "17C_remaining_significant_effects_statistics.csv"
  ),
  row.names = FALSE
)

write.csv(
  effect_data,
  file.path(
    OUT_DIR,
    "tables",
    "17C_remaining_significant_partial_effect_curves.csv"
  ),
  row.names = FALSE
)

# A common y-axis permits direct visual comparison of effect magnitude.
y_limits <- range(c(effect_data$lower, effect_data$upper), finite = TRUE)
y_span <- diff(y_limits)
if (!is.finite(y_span) || y_span <= 0) y_span <- 0.1
y_limits <- y_limits + c(-1, 1) * max(0.02, 0.07 * y_span)

make_panel <- function(variable, show_y_title) {
  curve_data <- effect_data %>%
    filter(as.character(.data$variable) == .env$variable) %>%
    arrange(x)

  term_stats <- remaining_significant %>%
    filter(.data$variable == .env$variable)

  if (nrow(curve_data) == 0L || nrow(term_stats) != 1L) {
    stop_with("Invalid panel data for ", variable, ".")
  }

  xmin <- suppressWarnings(as.numeric(curve_data$display_min[1L]))
  xmax <- suppressWarnings(as.numeric(curve_data$display_max[1L]))

  if (!is.finite(xmin) || !is.finite(xmax) || xmin >= xmax) {
    stop_with("Invalid display range for ", variable, ".")
  }

  rug_data <- data.frame(
    x = suppressWarnings(as.numeric(model_input[[variable]]))
  ) %>%
    filter(is.finite(x), x >= xmin, x <= xmax)

  annotation <- paste0(
    "F = ",
    formatC(as.numeric(term_stats$statistic[1L]), format = "f", digits = 1),
    ", ",
    format_p(term_stats$p_value[1L])
  )

  ggplot(curve_data, aes(x = x, y = effect)) +
    geom_hline(yintercept = 0, colour = "#9A9A9A", linewidth = 0.35) +
    geom_ribbon(
      aes(ymin = lower, ymax = upper),
      fill = "#9CA3AF",
      alpha = 0.34,
      colour = NA
    ) +
    geom_line(colour = "#111111", linewidth = 0.82) +
    geom_rug(
      data = rug_data,
      aes(x = x),
      inherit.aes = FALSE,
      sides = "b",
      alpha = 0.30,
      linewidth = 0.26
    ) +
    annotate(
      "text",
      x = min(curve_data$x, na.rm = TRUE),
      y = y_limits[1L] + 0.04 * diff(y_limits),
      label = annotation,
      hjust = 0,
      vjust = 0,
      size = 3.05,
      colour = "#222222"
    ) +
    coord_cartesian(ylim = y_limits, clip = "off") +
    labs(
      title = as.character(term_stats$panel_label[1L]),
      x = NULL,
      y = if (show_y_title) {
        "Partial effect on Shannon diversity"
      } else {
        NULL
      }
    ) +
    theme_minimal(base_size = 10, base_family = "sans") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#E5E5E5", linewidth = 0.35),
      axis.line = element_line(colour = "#555555", linewidth = 0.35),
      axis.ticks = element_line(colour = "#555555", linewidth = 0.30),
      axis.text = element_text(colour = "#333333", size = 8.5),
      axis.title.y = element_text(size = 9.3, margin = margin(r = 7)),
      plot.title = element_text(
        face = "bold",
        size = 9.6,
        hjust = 0.5,
        margin = margin(b = 5)
      ),
      plot.margin = margin(t = 6, r = 7, b = 5, l = 7)
    )
}

panels <- lapply(seq_along(remaining_variables), function(index) {
  make_panel(
    remaining_variables[index],
    show_y_title = ((index - 1L) %% 2L) == 0L
  )
})

# No composite title by design.
figure <- patchwork::wrap_plots(panels, ncol = 2)

n_rows <- ceiling(length(panels) / 2)
figure_width <- 8.8
figure_height <- max(3.2, 3.1 * n_rows)

figure_base <- file.path(
  OUT_DIR,
  "plots",
  "Figure_S_remaining_significant_effects_from_full_GAM"
)

ggsave(
  paste0(figure_base, ".png"), figure,
  width = figure_width, height = figure_height,
  units = "in", dpi = 600, bg = "white"
)

ggsave(
  paste0(figure_base, ".tiff"), figure,
  width = figure_width, height = figure_height,
  units = "in", dpi = 600, compression = "lzw", bg = "white"
)

ggsave(
  paste0(figure_base, ".pdf"), figure,
  width = figure_width, height = figure_height,
  units = "in", device = cairo_pdf, bg = "white"
)

summary_lines <- c(
  "Supplementary figure — remaining significant environmental effects",
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("Full-model source: ", FULL_MODEL_DIR),
  "No GAM was fitted or refitted by this script",
  paste0("Significance rule: P < ", SIGNIFICANCE_LEVEL),
  "Exclusion rule: selected_top4 = TRUE in the 17A ranking",
  paste0("Figure 2b variables excluded: ", paste(top4_variables, collapse = ", ")),
  paste0(
    "Remaining significant variables: ",
    paste(remaining_variables, collapse = ", ")
  ),
  paste0("Number of panels: ", length(remaining_variables)),
  "Panel order: descending partial-effect magnitude according to effect_rank",
  "Common y-axis across panels",
  "Y-axis: Partial effect on Shannon diversity",
  "Scale note: effects remain on the log1p(Shannon entropy) model scale",
  paste0("Output: ", OUT_DIR)
)

writeLines(summary_lines, file.path(OUT_DIR, "17C_run_summary.txt"))
capture.output(sessionInfo(), file = file.path(OUT_DIR, "logs", "sessionInfo.txt"))
writeLines(OUT_DIR, LATEST_POINTER)

msg("Completed successfully: ", OUT_DIR)
