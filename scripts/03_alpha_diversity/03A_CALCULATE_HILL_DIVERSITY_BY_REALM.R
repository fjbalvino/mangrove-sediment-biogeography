#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# Taxonomic Hill diversity across marine realms
#
# Uses an explicit sample-accession manifest and metadata joined by run_accession.
# Standardizes profiles to a common read count, calculates Hill q0/q1/q2,
# and plots metagenome-level q1 distributions. These realm comparisons are
# descriptive; primary coordinate-site inference is implemented in script 03E.

options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(phyloseq)
  library(ggplot2)
  library(dplyr)
  library(vegan)
})

msg <- function(...) {
  cat(
    sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    ..., "\n", sep = ""
  )
}

stop_with <- function(...) stop(paste0(...), call. = FALSE)

normalize_text <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "null", "NULL")] <- NA_character_
  x
}

normalize_numeric <- function(x) {
  suppressWarnings(as.numeric(normalize_text(x)))
}

format_p <- function(p, prefix = "P") {
  if (is.na(p)) return(paste0(prefix, " = NA"))
  if (p < 0.001) return(paste0(prefix, " < 0.001"))
  paste0(prefix, " = ", formatC(p, format = "f", digits = 3))
}

read_accession_file <- function(path, label) {
  if (!file.exists(path)) {
    stop_with(label, " file not found: ", path)
  }

  values <- normalize_text(readLines(path, warn = FALSE))
  values <- values[!is.na(values) & !grepl("^#", values)]

  if (length(values) == 0L) {
    stop_with(label, " file contains no accessions: ", path)
  }

  if (anyDuplicated(values)) {
    stop_with(label, " file contains duplicated accessions: ", path)
  }

  values
}

read_integer_env <- function(name, default = NA_integer_) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  parsed <- suppressWarnings(as.integer(value))
  if (length(parsed) != 1L || is.na(parsed)) {
    stop_with("Environment variable ", name, " is not an integer: ", value)
  }
  parsed
}

# ==============================================================================
# 1. Configuration
# ==============================================================================

BASE_DIR <- Sys.getenv("BIOGEOGRAPHY_BASE", unset = ".")
OUTPUT_ROOT <- Sys.getenv("BIOGEOGRAPHY_RESULTS", unset = file.path(BASE_DIR, "results"))
INPUT_RDS <- Sys.getenv("MAGENTA_PHYLOSEQ_405", unset = file.path(BASE_DIR, "data", "taxonomic_counts.rds"))
METADATA_CSV <- Sys.getenv("MAGENTA_METADATA_CSV", unset = file.path(BASE_DIR, "data", "environmental_metadata.csv"))
SAMPLE_MANIFEST <- Sys.getenv("ALPHA_SAMPLE_MANIFEST", unset = file.path(BASE_DIR, "config", "samples", "taxonomic_diversity_accessions.txt"))
ANALYSIS_ACCESSIONS <- read_accession_file(SAMPLE_MANIFEST, "Taxonomic diversity sample manifest")
EXPECTED_ANALYSIS_N <- length(ANALYSIS_ACCESSIONS)
SEED_RAREFACTION <- 42L
SEED_POINTS <- 3030L
ALPHA <- 0.05

EXPECTED_REALMS <- c(
  "Central Indo-Pacific",
  "Indo-Pacific",
  "Temperate Northern Atlantic",
  "Tropical Atlantic",
  "Tropical Eastern Pacific",
  "Western Indo-Pacific"
)

REALM_CODES <- c(
  "Central Indo-Pacific" = "CIP",
  "Indo-Pacific" = "IP",
  "Temperate Northern Atlantic" = "TNA",
  "Tropical Atlantic" = "TA",
  "Tropical Eastern Pacific" = "TEP",
  "Western Indo-Pacific" = "WIP"
)

REALM_COLORS <- c(
  "Central Indo-Pacific" = "#66C2A5",
  "Indo-Pacific" = "#FC8D62",
  "Temperate Northern Atlantic" = "#8DA0CB",
  "Tropical Atlantic" = "#E78AC3",
  "Tropical Eastern Pacific" = "#A6D854",
  "Western Indo-Pacific" = "#FFD92F"
)

run_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUT_DIR <- Sys.getenv("ALPHA_OUTPUT_DIR", unset = file.path(OUTPUT_ROOT, paste0("TAXONOMIC_ALPHA_", run_stamp)))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_DIR <- normalizePath(OUT_DIR, mustWork = TRUE)
msg("Starting taxonomic Hill diversity analysis")
msg("Phyloseq: ", INPUT_RDS)
msg("Metadata: ", METADATA_CSV)
msg("Sample manifest: ", SAMPLE_MANIFEST)
msg("Output: ", OUT_DIR)

if (!file.exists(INPUT_RDS)) {
  stop_with("Phyloseq RDS not found: ", INPUT_RDS)
}

if (!file.exists(METADATA_CSV)) {
  stop_with("Metadata CSV not found: ", METADATA_CSV)
}

# ==============================================================================
# 2. Read and validate external metadata
# ==============================================================================

metadata_raw <- read.csv(
  METADATA_CSV,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA", "N/A", "null", "NULL")
)

names(metadata_raw) <- trimws(names(metadata_raw))

required_metadata <- c("run_accession", "Country", "REALM")
missing_metadata <- setdiff(required_metadata, names(metadata_raw))

if (length(missing_metadata) > 0L) {
  stop_with(
    "Metadata are missing required columns: ",
    paste(missing_metadata, collapse = ", ")
  )
}

for (column_name in c("Locality", "ECOREGION", "PROVINCE", "lat", "lon")) {
  if (!column_name %in% names(metadata_raw)) {
    metadata_raw[[column_name]] <- NA
  }
}

metadata_raw <- metadata_raw %>%
  mutate(
    run_accession = normalize_text(run_accession),
    Country = normalize_text(Country),
    REALM = normalize_text(REALM),
    Locality = normalize_text(Locality),
    ECOREGION = normalize_text(ECOREGION),
    PROVINCE = normalize_text(PROVINCE),
    lat = normalize_numeric(lat),
    lon = normalize_numeric(lon)
  ) %>%
  filter(!is.na(run_accession))

duplicate_ids <- metadata_raw %>%
  count(run_accession, name = "metadata_rows") %>%
  filter(metadata_rows > 1L)

write.csv(
  duplicate_ids,
  file.path(OUT_DIR, "04a_metadata_duplicate_run_accessions.csv"),
  row.names = FALSE
)

if (nrow(duplicate_ids) > 0L) {
  duplicate_conflicts <- metadata_raw %>%
    filter(run_accession %in% duplicate_ids$run_accession) %>%
    group_by(run_accession) %>%
    summarise(
      n_country = n_distinct(Country[!is.na(Country)]),
      n_realm = n_distinct(REALM[!is.na(REALM)]),
      n_locality = n_distinct(Locality[!is.na(Locality)]),
      .groups = "drop"
    ) %>%
    filter(n_country > 1L | n_realm > 1L | n_locality > 1L)

  write.csv(
    duplicate_conflicts,
    file.path(OUT_DIR, "04a_metadata_duplicate_conflicts.csv"),
    row.names = FALSE
  )

  if (nrow(duplicate_conflicts) > 0L) {
    stop_with(
      "Conflicting metadata were found for duplicated run_accession IDs. ",
      "See 04a_metadata_duplicate_conflicts.csv."
    )
  }
}

metadata <- metadata_raw %>%
  distinct(run_accession, .keep_all = TRUE)

# ==============================================================================
# 3. Resolve metadata against phyloseq and select the analysis samples
# ==============================================================================

ps <- readRDS(INPUT_RDS)

if (!inherits(ps, "phyloseq")) {
  stop_with("INPUT_RDS is not a phyloseq object.")
}


ps_metadata <- as(sample_data(ps), "data.frame")
ps_metadata$.ps_sample_id <- rownames(ps_metadata)

if (!"run_accession" %in% names(ps_metadata)) {
  stop_with(
    "The phyloseq sample_data table has no run_accession column. ",
    "No alternative join key will be inferred."
  )
}

ps_text_column <- function(column_name) {
  if (column_name %in% names(ps_metadata)) {
    normalize_text(ps_metadata[[column_name]])
  } else {
    rep(NA_character_, nrow(ps_metadata))
  }
}

ps_numeric_column <- function(column_name) {
  if (column_name %in% names(ps_metadata)) {
    normalize_numeric(ps_metadata[[column_name]])
  } else {
    rep(NA_real_, nrow(ps_metadata))
  }
}

ps_fallback <- data.frame(
  .ps_sample_id = ps_metadata$.ps_sample_id,
  run_accession = normalize_text(ps_metadata$run_accession),
  Country_phyloseq = ps_text_column("Country"),
  REALM_phyloseq = ps_text_column("REALM"),
  Locality_phyloseq = ps_text_column("Locality"),
  ECOREGION_phyloseq = ps_text_column("ECOREGION"),
  PROVINCE_phyloseq = ps_text_column("PROVINCE"),
  lat_phyloseq = ps_numeric_column("lat"),
  lon_phyloseq = ps_numeric_column("lon"),
  stringsAsFactors = FALSE
)

duplicated_ps_key <- ps_fallback %>%
  count(run_accession, name = "phyloseq_samples") %>%
  filter(is.na(run_accession) | phyloseq_samples > 1L)

if (nrow(duplicated_ps_key) > 0L) {
  write.csv(
    duplicated_ps_key,
    file.path(OUT_DIR, "04a_phyloseq_invalid_run_accessions.csv"),
    row.names = FALSE
  )
  stop_with("run_accession is missing or duplicated in phyloseq sample_data.")
}

metadata_not_in_phyloseq <- metadata %>%
  anti_join(ps_fallback %>% select(run_accession), by = "run_accession")

phyloseq_not_in_metadata <- ps_fallback %>%
  anti_join(metadata %>% select(run_accession), by = "run_accession")

write.csv(
  metadata_not_in_phyloseq,
  file.path(OUT_DIR, "04a_metadata_only_run_accessions.csv"),
  row.names = FALSE
)

write.csv(
  phyloseq_not_in_metadata,
  file.path(OUT_DIR, "04a_phyloseq_only_run_accessions.csv"),
  row.names = FALSE
)

sample_map <- ps_fallback %>%
  left_join(metadata, by = "run_accession") %>%
  mutate(
    Country_external = normalize_text(Country),
    REALM_external = normalize_text(REALM),
    Locality_external = normalize_text(Locality),
    ECOREGION_external = normalize_text(ECOREGION),
    PROVINCE_external = normalize_text(PROVINCE),
    lat_external = normalize_numeric(lat),
    lon_external = normalize_numeric(lon),
    Country = coalesce(Country_external, Country_phyloseq),
    REALM = coalesce(REALM_external, REALM_phyloseq),
    Locality = coalesce(Locality_external, Locality_phyloseq),
    ECOREGION = coalesce(ECOREGION_external, ECOREGION_phyloseq),
    PROVINCE = coalesce(PROVINCE_external, PROVINCE_phyloseq),
    lat = coalesce(lat_external, lat_phyloseq),
    lon = coalesce(lon_external, lon_phyloseq),
    metadata_resolution = case_when(
      !is.na(Country_external) & !is.na(REALM_external) &
        is.finite(lon_external) ~ "external_metadata_complete",
      !is.na(REALM_external) ~ "external_realm_partially_complete",
      is.na(REALM_external) & !is.na(REALM_phyloseq) ~
        "phyloseq_fallback_for_realm",
      TRUE ~ "unresolved"
    )
  )

write.csv(
  sample_map %>%
    filter(metadata_resolution != "external_metadata_complete") %>%
    select(
      .ps_sample_id, run_accession,
      Country_external, REALM_external, Locality_external,
      lat_external, lon_external,
      Country_phyloseq, REALM_phyloseq, Locality_phyloseq,
      lat_phyloseq, lon_phyloseq,
      Country, REALM, Locality, lat, lon,
      metadata_resolution
    ),
  file.path(OUT_DIR, "04a_metadata_resolution_audit.csv"),
  row.names = FALSE
)

unresolved_samples <- sample_map %>%
  filter(is.na(REALM))

write.csv(
  unresolved_samples,
  file.path(OUT_DIR, "04a_unresolved_phyloseq_samples.csv"),
  row.names = FALSE
)

if (nrow(unresolved_samples) > 0L) {
  stop_with(
    nrow(unresolved_samples),
    " phyloseq samples still lack REALM after metadata resolution."
  )
}

missing_accessions <- setdiff(ANALYSIS_ACCESSIONS, sample_map$run_accession)
if (length(missing_accessions)) {
  stop_with("Manifest accessions absent from phyloseq: ", paste(missing_accessions, collapse = ", "))
}
sample_map <- sample_map %>% mutate(
  keep_sample = run_accession %in% ANALYSIS_ACCESSIONS,
  filter_status = ifelse(keep_sample, "included_analysis_manifest", "not_in_analysis_manifest")
)
write.csv(sample_map %>% select(.ps_sample_id, run_accession, Country, Locality, REALM, keep_sample, filter_status),
          file.path(OUT_DIR, "04a_sample_selection_manifest.csv"), row.names = FALSE)
n_before <- nsamples(ps)
keep_sample_names <- sample_map$.ps_sample_id[sample_map$keep_sample]
# Pruning preserves the source matrix order used for the published calculations.
ps_filtered <- prune_samples(sample_names(ps) %in% keep_sample_names, ps)
ps_filtered <- prune_taxa(taxa_sums(ps_filtered) > 0, ps_filtered)
n_after <- nsamples(ps_filtered)
if (n_after != EXPECTED_ANALYSIS_N) stop_with("Selected sample count does not match the accession manifest.")
msg("Analysis samples: ", n_after, " of ", n_before, " supplied profiles")

# ==============================================================================
# 4. Rarefy after filtering and calculate Hill numbers
# ==============================================================================

library_sizes_before <- data.frame(
  .ps_sample_id = sample_names(ps_filtered),
  reads_before_rarefaction = as.numeric(sample_sums(ps_filtered)),
  stringsAsFactors = FALSE
)

min_depth <- min(library_sizes_before$reads_before_rarefaction)

if (!is.finite(min_depth) || min_depth <= 0) {
  stop_with("Invalid minimum library size: ", min_depth)
}

set.seed(SEED_RAREFACTION)

ps_rare <- rarefy_even_depth(
  ps_filtered,
  sample.size = min_depth,
  replace = FALSE,
  rngseed = SEED_RAREFACTION,
  verbose = FALSE
)

otu <- as(otu_table(ps_rare), "matrix")

if (taxa_are_rows(ps_rare)) {
  otu <- t(otu)
}

alpha_df <- data.frame(
  .ps_sample_id = rownames(otu),
  q0 = vegan::specnumber(otu),
  q1 = exp(vegan::diversity(otu, index = "shannon")),
  q2 = vegan::diversity(otu, index = "invsimpson"),
  stringsAsFactors = FALSE,
  row.names = NULL
)

plot_df <- alpha_df %>%
  left_join(
    sample_map %>% filter(keep_sample),
    by = ".ps_sample_id"
  )

if (any(!is.finite(plot_df$q1))) {
  stop_with("Non-finite Hill q=1 values were generated.")
}

unexpected_realms <- setdiff(unique(plot_df$REALM), EXPECTED_REALMS)
missing_realms <- setdiff(EXPECTED_REALMS, unique(plot_df$REALM))

if (length(unexpected_realms) > 0L) {
  stop_with(
    "Unexpected REALM labels: ",
    paste(unexpected_realms, collapse = ", ")
  )
}

if (length(missing_realms) > 0L) {
  stop_with(
    "Expected REALM labels absent after filtering: ",
    paste(missing_realms, collapse = ", ")
  )
}

# ==============================================================================
# 5. Order realms west-to-east using unique sampling sites
# ==============================================================================

site_longitudes <- plot_df %>%
  filter(is.finite(lat), is.finite(lon)) %>%
  mutate(
    site_key = paste(
      coalesce(Locality, "unknown_locality"),
      round(lat, 6),
      round(lon, 6),
      sep = "|"
    )
  ) %>%
  distinct(REALM, site_key, .keep_all = TRUE)

realm_longitude_order <- site_longitudes %>%
  group_by(REALM) %>%
  summarise(
    n_unique_sites_with_coordinates = n(),
    mean_longitude = mean(lon),
    median_longitude = median(lon),
    min_longitude = min(lon),
    max_longitude = max(lon),
    .groups = "drop"
  ) %>%
  arrange(median_longitude) %>%
  mutate(west_to_east_position = row_number())

missing_longitude_realms <- setdiff(
  EXPECTED_REALMS,
  realm_longitude_order$REALM
)

if (length(missing_longitude_realms) > 0L) {
  stop_with(
    "Cannot order realms because these lack valid coordinates: ",
    paste(missing_longitude_realms, collapse = ", ")
  )
}

REALM_LEVELS <- realm_longitude_order$REALM

plot_df <- plot_df %>%
  mutate(
    REALM = factor(REALM, levels = REALM_LEVELS),
    realm_index = as.numeric(REALM)
  )

set.seed(SEED_POINTS)
plot_df$point_x <- plot_df$realm_index + runif(nrow(plot_df), 0.09, 0.28)

realm_counts <- plot_df %>%
  count(REALM, name = "n") %>%
  arrange(REALM)

x_labels <- paste0(
  unname(REALM_CODES[as.character(realm_counts$REALM)]),
  "\nn = ", realm_counts$n
)

write.csv(
  realm_longitude_order,
  file.path(OUT_DIR, "04a_realm_longitude_order_west_to_east.csv"),
  row.names = FALSE
)

# ==============================================================================
# 6. Global Kruskal-Wallis and Dunn-Holm post hoc tests
# ==============================================================================

kw_test <- kruskal.test(q1 ~ REALM, data = plot_df)

# Rank-based epsilon-squared effect size for the Kruskal-Wallis comparison:
# epsilon^2_H = (H - k + 1) / (n - k). It is exported explicitly because the
# manuscript audit found that an earlier text mentioned epsilon-squared without
# a traceable table.
kw_epsilon_squared <- (
  unname(kw_test$statistic) - nlevels(plot_df$REALM) + 1
) / (
  nrow(plot_df) - nlevels(plot_df$REALM)
)

kw_table <- data.frame(
  method = unname(kw_test$method),
  statistic = unname(kw_test$statistic),
  df = unname(kw_test$parameter),
  p_value = kw_test$p.value,
  epsilon_squared_H = kw_epsilon_squared,
  n_samples = nrow(plot_df),
  n_realms = nlevels(plot_df$REALM),
  stringsAsFactors = FALSE
)

dunn_holm_test <- function(values, groups) {
  keep <- is.finite(values) & !is.na(groups)
  values <- values[keep]
  groups <- droplevels(factor(groups[keep]))

  n_total <- length(values)

  if (n_total < 2L || nlevels(groups) < 2L) {
    stop_with("Dunn test requires at least two nonempty groups.")
  }

  pooled_ranks <- rank(values, ties.method = "average")
  tie_sizes <- as.numeric(table(values))
  tie_correction <- sum(tie_sizes^3 - tie_sizes)

  variance_base <-
    n_total * (n_total + 1) / 12 -
    tie_correction / (12 * (n_total - 1))

  group_statistics <- data.frame(
    group = levels(groups),
    n = as.numeric(table(groups)[levels(groups)]),
    mean_rank = as.numeric(
      tapply(pooled_ranks, groups, mean)[levels(groups)]
    ),
    stringsAsFactors = FALSE
  )

  group_pairs <- combn(levels(groups), 2L, simplify = FALSE)

  results <- bind_rows(lapply(group_pairs, function(pair) {
    first <- group_statistics[group_statistics$group == pair[1], ]
    second <- group_statistics[group_statistics$group == pair[2], ]

    standard_error <- sqrt(
      variance_base * (1 / first$n + 1 / second$n)
    )

    z_value <- (first$mean_rank - second$mean_rank) / standard_error
    p_raw <- 2 * pnorm(-abs(z_value))

    data.frame(
      group_1 = pair[1],
      group_2 = pair[2],
      n_1 = first$n,
      n_2 = second$n,
      mean_rank_1 = first$mean_rank,
      mean_rank_2 = second$mean_rank,
      z = z_value,
      p_raw = p_raw,
      stringsAsFactors = FALSE
    )
  }))

  results %>%
    mutate(
      p_adjust_holm = p.adjust(p_raw, method = "holm"),
      significant_005 = p_adjust_holm <= ALPHA,
      significance = case_when(
        p_adjust_holm <= 0.001 ~ "***",
        p_adjust_holm <= 0.01 ~ "**",
        p_adjust_holm <= 0.05 ~ "*",
        TRUE ~ "ns"
      ),
      group_1_code = unname(REALM_CODES[group_1]),
      group_2_code = unname(REALM_CODES[group_2])
    ) %>%
    arrange(p_adjust_holm)
}

dunn_results <- dunn_holm_test(plot_df$q1, plot_df$REALM)

compact_letters <- function(pairwise_results, group_order, alpha = 0.05) {
  n_groups <- length(group_order)

  p_matrix <- matrix(
    1,
    nrow = n_groups,
    ncol = n_groups,
    dimnames = list(group_order, group_order)
  )

  for (row_index in seq_len(nrow(pairwise_results))) {
    group_1 <- pairwise_results$group_1[row_index]
    group_2 <- pairwise_results$group_2[row_index]
    p_value <- pairwise_results$p_adjust_holm[row_index]
    p_matrix[group_1, group_2] <- p_value
    p_matrix[group_2, group_1] <- p_value
  }

  nonsignificant <- p_matrix > alpha
  diag(nonsignificant) <- TRUE

  candidate_cliques <- list()

  for (mask in seq_len(2^n_groups - 1L)) {
    included <- as.logical(intToBits(mask)[seq_len(n_groups)])
    members <- group_order[included]
    submatrix <- nonsignificant[members, members, drop = FALSE]

    if (all(submatrix)) {
      candidate_cliques[[length(candidate_cliques) + 1L]] <- members
    }
  }

  maximal_cliques <- Filter(
    function(candidate) {
      !any(vapply(
        candidate_cliques,
        function(other) {
          length(other) > length(candidate) && all(candidate %in% other)
        },
        logical(1)
      ))
    },
    candidate_cliques
  )

  clique_key <- vapply(
    maximal_cliques,
    function(x) paste(sort(x), collapse = "|"),
    character(1)
  )
  maximal_cliques <- maximal_cliques[!duplicated(clique_key)]

  clique_order <- order(vapply(
    maximal_cliques,
    function(x) min(match(x, group_order)),
    numeric(1)
  ))
  maximal_cliques <- maximal_cliques[clique_order]

  letter_pool <- c(letters, paste0("a", seq_len(100)))
  group_letters <- setNames(rep("", n_groups), group_order)

  for (clique_index in seq_along(maximal_cliques)) {
    label <- letter_pool[clique_index]
    for (group_name in maximal_cliques[[clique_index]]) {
      group_letters[group_name] <- paste0(group_letters[group_name], label)
    }
  }

  data.frame(
    REALM = names(group_letters),
    dunn_holm_letters = unname(group_letters),
    stringsAsFactors = FALSE
  )
}

realm_descriptives <- plot_df %>%
  group_by(REALM) %>%
  summarise(
    n = n(),
    median_q1 = median(q1),
    mean_q1 = mean(q1),
    sd_q1 = sd(q1),
    q1_25 = quantile(q1, 0.25),
    q1_75 = quantile(q1, 0.75),
    max_q1 = max(q1),
    .groups = "drop"
  )

letter_group_order <- realm_descriptives %>%
  arrange(desc(median_q1)) %>%
  pull(REALM) %>%
  as.character()

letter_table <- compact_letters(
  dunn_results,
  group_order = letter_group_order,
  alpha = ALPHA
)

realm_descriptives <- realm_descriptives %>%
  left_join(letter_table, by = "REALM") %>%
  mutate(
    REALM = factor(REALM, levels = REALM_LEVELS),
    realm_index = as.numeric(REALM)
  ) %>%
  arrange(REALM)

# Coordinate-site inference is calculated by 03E using valid geographic
# coordinates only, with a minimum of five sites per realm.

# ==============================================================================
# 7. Raincloud plot
# ==============================================================================

GeomHalfViolin <- ggproto(
  "GeomHalfViolin",
  GeomViolin,
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data <- transform(
      data,
      xminv = x - violinwidth * (x - xmin),
      xmaxv = x
    )

    polygon_data <- rbind(
      transform(data, x = xminv)[order(data$y), ],
      transform(data, x = xmaxv)[order(data$y, decreasing = TRUE), ]
    )

    polygon_data <- rbind(polygon_data, polygon_data[1, ])
    GeomPolygon$draw_panel(polygon_data, ...)
  }
)

geom_half_violin <- function(
    mapping = NULL,
    data = NULL,
    stat = "ydensity",
    position = "identity",
    ...,
    trim = TRUE,
    scale = "area",
    na.rm = FALSE,
    show.legend = NA,
    inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomHalfViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      na.rm = na.rm,
      ...
    )
  )
}

y_data_max <- max(plot_df$q1)
realm_descriptives$letter_y <-
  realm_descriptives$max_q1 + 0.035 * y_data_max

global_label_y <- y_data_max * 1.16
y_limit <- max(global_label_y * 1.04, max(realm_descriptives$letter_y) * 1.05)

p <- ggplot(
  plot_df,
  aes(x = realm_index, y = q1, group = REALM, fill = REALM)
) +
  geom_half_violin(
    width = 0.72,
    trim = TRUE,
    scale = "width",
    alpha = 0.55,
    colour = NA
  ) +
  geom_boxplot(
    aes(x = realm_index, y = q1, group = REALM, fill = REALM),
    width = 0.15,
    outlier.shape = NA,
    alpha = 0.82,
    colour = "#333333",
    linewidth = 0.38
  ) +
  geom_point(
    aes(x = point_x, y = q1, fill = REALM),
    inherit.aes = FALSE,
    shape = 21,
    size = 1.55,
    stroke = 0.22,
    colour = "white",
    alpha = 0.86
  ) +
  geom_text(
    data = realm_descriptives,
    aes(
      x = realm_index,
      y = letter_y,
      label = dunn_holm_letters
    ),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 3.7,
    colour = "#222222"
  ) +
  annotate(
    "text",
    x = length(REALM_LEVELS) + 0.34,
    y = global_label_y,
    label = paste0(
      "Kruskal-Wallis: H = ",
      formatC(unname(kw_test$statistic), format = "f", digits = 1),
      ", ",
      format_p(kw_test$p.value)
    ),
    hjust = 1,
    vjust = 1,
    size = 3.15,
    colour = "#222222"
  ) +
  scale_fill_manual(values = REALM_COLORS, drop = FALSE) +
  scale_x_continuous(
    breaks = seq_along(REALM_LEVELS),
    labels = x_labels,
    limits = c(0.48, length(REALM_LEVELS) + 0.48),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    limits = c(0, y_limit),
    expand = expansion(mult = c(0, 0)),
    labels = scales::label_comma()
  ) +
  labs(
    tag = "a",
    title = "Taxonomic diversity across marine realms",
    x = "Marine realms ordered west to east",
    y = "Hill diversity (q = 1)"
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 10.5, base_family = "sans") +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(
      colour = "#E8E8E8",
      linewidth = 0.35
    ),
    axis.line = element_line(colour = "#333333", linewidth = 0.42),
    axis.ticks = element_line(colour = "#333333", linewidth = 0.35),
    axis.ticks.length = grid::unit(1.8, "mm"),
    axis.text.x = element_text(
      colour = "#222222",
      size = 9.1,
      lineheight = 1.05,
      margin = margin(t = 6)
    ),
    axis.text.y = element_text(colour = "#333333", size = 9.2),
    axis.title.x = element_text(
      colour = "#333333",
      size = 9.4,
      margin = margin(t = 9)
    ),
    axis.title.y = element_text(
      colour = "#222222",
      size = 10.7,
      margin = margin(r = 8)
    ),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 11.3,
      margin = margin(b = 10)
    ),
    plot.tag = element_text(face = "bold", size = 15),
    plot.tag.position = c(0.012, 0.992),
    plot.margin = margin(t = 9, r = 10, b = 6, l = 8)
  )

# ==============================================================================
# 8. Export figure and audit tables
# ==============================================================================

png_file <- file.path(
  OUT_DIR,
  "Figure_2a_Hill_q1_by_marine_realm.png"
)
tif_file <- file.path(
  OUT_DIR,
  "Figure_2a_Hill_q1_by_marine_realm.tiff"
)
pdf_file <- file.path(
  OUT_DIR,
  "Figure_2a_Hill_q1_by_marine_realm.pdf"
)

ggsave(
  png_file,
  p,
  width = 8.2,
  height = 5.35,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  tif_file,
  p,
  width = 8.2,
  height = 5.35,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

ggsave(
  pdf_file,
  p,
  width = 8.2,
  height = 5.35,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

write.csv(
  plot_df %>%
    select(
      .ps_sample_id,
      run_accession,
      Country,
      Locality,
      ECOREGION,
      PROVINCE,
      REALM,
      lat,
      lon,
      metadata_resolution,
      q0,
      q1,
      q2
    ),
  file.path(OUT_DIR, "04a_alpha_diversity_filtered_samplewise.csv"),
  row.names = FALSE
)

write.csv(
  realm_counts,
  file.path(OUT_DIR, "04a_sample_counts_by_realm.csv"),
  row.names = FALSE
)

write.csv(
  realm_descriptives,
  file.path(OUT_DIR, "04a_realm_q1_descriptive_statistics.csv"),
  row.names = FALSE
)

write.csv(
  kw_table,
  file.path(OUT_DIR, "04a_Kruskal_Wallis_q1.csv"),
  row.names = FALSE
)

write.csv(
  dunn_results,
  file.path(OUT_DIR, "04a_Dunn_Holm_q1.csv"),
  row.names = FALSE
)

write.csv(
  letter_table,
  file.path(OUT_DIR, "04a_compact_letters_Dunn_Holm_q1.csv"),
  row.names = FALSE
)

write.csv(
  library_sizes_before,
  file.path(OUT_DIR, "04a_library_sizes_before_rarefaction.csv"),
  row.names = FALSE
)

summary_lines <- c(
  "Taxonomic Hill diversity across marine realms",
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("Phyloseq object: ", INPUT_RDS),
  paste0("Metadata file: ", METADATA_CSV),
  paste0("Sample manifest: ", SAMPLE_MANIFEST),
  "Join key: run_accession",
  paste0("Supplied samples: ", n_before),
  paste0("Analysis samples: ", n_after),
  paste0("Nonzero taxonomic features: ", ntaxa(ps_filtered)),
  paste0("Standardized read depth: ", min_depth),
  "Input counts may already be standardized; supplied library sizes are exported.",
  "Metagenome-level tests describe the plotted distributions; use 03E for primary site-level inference.",
  paste0("Metagenome-level Kruskal-Wallis P: ", kw_test$p.value),
  paste0("Output directory: ", OUT_DIR)
)
writeLines(summary_lines, file.path(OUT_DIR, "04a_run_summary.txt"))
capture.output(sessionInfo(), file = file.path(OUT_DIR, "sessionInfo.txt"))
dir.create(OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)
writeLines(OUT_DIR, file.path(OUTPUT_ROOT, "LATEST_TAXONOMIC_ALPHA.txt"))
msg("Completed successfully: ", OUT_DIR)
