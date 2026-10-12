#!/usr/bin/env Rscript
# Small regression checks using expressions extracted from the analysis scripts.
# No community data or external R packages are loaded.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) != 1L) stop("Run this test with Rscript.")
repo <- dirname(dirname(normalizePath(sub("^--file=", "", file_arg))))
beta_dir <- file.path(repo, "scripts", "04_beta_diversity")
decay <- parse(file.path(beta_dir, "04B_DISTANCE_DECAY_AND_MANTEL.R"))
vpa <- parse(file.path(beta_dir, "04D_DBRDA_VARIANCE_PARTITIONING.R"))

is_call <- function(x, name) is.call(x) && identical(x[[1L]], as.name(name))
assignment <- function(expressions, name) {
  matches <- Filter(function(x) {
    is_call(x, "<-") && identical(x[[2L]], as.name(name))
  }, as.list(expressions))
  if (length(matches) != 1L) stop("Expected one assignment for ", name)
  matches[[1L]][[3L]]
}
# Index expressions may have missing arguments, and function formals need not
# be traversed. All targets below are ordinary calls or function bodies.
find_calls <- function(x, predicate) {
  found <- if (is.call(x) && predicate(x)) list(x) else list()
  if (!is.call(x)) return(found)
  if (is_call(x, "[") || is_call(x, "[[")) return(found)
  children <- if (is_call(x, "function")) list(x[[3L]]) else as.list(x)[-1L]
  for (child in children) found <- c(found, find_calls(child, predicate))
  found
}
only <- function(x, label) {
  if (length(x) != 1L) stop("Expected exactly one ", label)
  x[[1L]]
}

# Execute the actual selection and identity checks inside the distance-decay
# resampling callback, including singleton sites with numeric indices > 1.
resampling <- assignment(decay, "site_balanced_results")
callback <- only(find_calls(resampling, function(x) {
  is_call(x, "function") && "selected" %in% all.vars(x[[3L]])
}), "site-resampling callback")
callback_body <- as.list(callback[[3L]])[-1L]
selection_index <- which(vapply(callback_body, function(x) {
  is_call(x, "<-") && identical(x[[2L]], as.name("selected"))
}, logical(1L)))
stopifnot(length(selection_index) == 1L)
selection <- callback_body[[selection_index]]
identity_checks <- callback_body[[selection_index + 1L]]
stopifnot(is_call(identity_checks, "stopifnot"))

sampling_env <- new.env(parent = baseenv())
sampling_env$site_groups <- list(
  singleton_a = 17L, singleton_b = 30L,
  multiple_a = c(2L, 5L, 9L), multiple_b = c(11L, 12L)
)
sampling_env$sample_table <- data.frame(site_id = rep(NA_character_, 30L))
for (site in names(sampling_env$site_groups)) {
  sampling_env$sample_table$site_id[sampling_env$site_groups[[site]]] <- site
}
set.seed(3721L)
for (i in seq_len(1000L)) {
  eval(selection, sampling_env)
  eval(identity_checks, sampling_env)
  selected <- sampling_env$selected
  stopifnot(identical(names(selected), names(sampling_env$site_groups)))
  stopifnot(all(vapply(names(selected), function(site) {
    selected[[site]] %in% sampling_env$site_groups[[site]]
  }, logical(1L))))
}
cat("PASS: actual site-resampling code preserves membership in 1,000 draws.\n")

# Exercise the actual base-R predictor helper with dependent, constant and
# nonfinite columns, then with no usable predictors.
make_full_rank <- eval(assignment(vpa, "make_full_rank"), envir = baseenv())
i <- seq_len(6L)
x <- cbind(linear = i, duplicate = 2 * i, quadratic = i^2,
           constant = 1, invalid = NA_real_)
reduced <- make_full_rank(x, "test")
stopifnot(nrow(reduced) == 6L, ncol(reduced) == 2L,
          qr(reduced)$rank == 2L, all(is.finite(reduced)))
empty <- make_full_rank(cbind(constant = rep(1, 6), invalid = rep(NA_real_, 6)), "test")
stopifnot(identical(dim(empty), c(6L, 0L)))
empty_input <- make_full_rank(matrix(numeric(), nrow = 6L, ncol = 0L), "test")
stopifnot(identical(dim(empty_input), c(6L, 0L)))

# Evaluate the rank formula used by the real diagnostics on explanatory
# columns completely spanned by the study design, including an intercept.
rank_table <- assignment(vpa, "rank_diagnostics")
rank_expression <- as.list(rank_table)$rank
stopifnot(!is.null(rank_expression))
rank_env <- list2env(list(ENV = matrix(i), SPACE = matrix(i^2),
                         REALM = matrix(i), STUDY = cbind(i, i^2)), parent = baseenv())
ranks <- eval(rank_expression, rank_env)
stopifnot(length(ranks) == 6L, ranks[[4L]] == 2L,
          ranks[[6L]] - ranks[[4L]] == 0L)
cat("PASS: actual rank helper and diagnostics handle redundant predictors.\n")

# Check the actual inference guard without claiming to fit or validate a
# vegan model: zero constrained rank or no residual df must prevent testing.
fit_function <- assignment(vpa, "fit_dbrda_block")
rank_guard <- only(find_calls(fit_function, function(x) {
  is_call(x, "if") && setequal(all.vars(x[[2L]]), c("model_rank", "residual_df"))
}), "estimability guard")
cases <- data.frame(model_rank = c(0L, 2L, 0L, 2L),
                    residual_df = c(20L, 0L, 0L, 20L),
                    expected = c(FALSE, FALSE, FALSE, TRUE))
for (j in seq_len(nrow(cases))) {
  actual <- eval(rank_guard[[2L]], list2env(as.list(cases[j, ]), parent = baseenv()))
  stopifnot(identical(actual, cases$expected[[j]]))
}
# The non-estimable branch must not report performed permutations.
no_test_branch <- rank_guard[[4L]]
stopifnot(identical(assignment(as.list(no_test_branch)[-1L], "performed_permutations"), 0L))

# Static AST check for exported shared-fraction labels.
# This validates exported labels, not the numerical variance partition.
recode_call <- only(find_calls(assignment(vpa, "individual_fractions"), function(x) {
  identical(x[[1L]], quote(dplyr::recode))
}), "variance-fraction recode call")
labels <- as.list(recode_call)
stopifnot(identical(labels[["[e]"]], "Shared space + realm"),
          identical(labels[["[f]"]], "Shared environment + realm"))
cat("PASS: estimability guard and shared-fraction label contracts.\n")
