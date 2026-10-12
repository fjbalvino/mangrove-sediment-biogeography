#!/usr/bin/env Rscript
# Supplementary effects with BH-adjusted P < 0.05. Uses the saved full model and does not refit it.
GAM_PLOT_SCOPE <- "additional"
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_arg) != 1L) stop("Run this script with Rscript.", call. = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE))
source(file.path(script_dir, "_gam_outputs.R"), chdir = FALSE)
