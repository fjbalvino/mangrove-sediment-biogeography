#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# ==============================================================================
# 05C_SUMMARIZE_CORE_ASSOCIATED_PHYLA.R
# Summarize phylum composition within the core-associated fraction from 05A.
#
# Core membership uses >=50 total reads, >=5 REALM || Locality units, occurrence
# in all localities of at least one realm with >=3 localities, and mean relative
# abundance >=0.0001 (0.01%) within that realm.
#
# Locality-balanced composition is calculated by pooling sample counts within
# each locality, calculating each phylum's share of that locality's core reads,
# and averaging those shares with equal weight across REALM || Locality units.
# It does not give equal weight to individual metagenomes or to marine realms.
# Pooled-read composition is exported separately. The plot shows the five phyla
# with the largest equal-locality means and groups the rest as 'Other phyla'.
# ==============================================================================

options(stringsAsFactors = FALSE, scipen = 999, bitmapType = "cairo")
Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
set.seed(42)

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(scales)
})

msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      ..., "\n", sep = "")
}

stopf <- function(...) stop(sprintf(...), call. = FALSE)

extract_phylum <- function(taxon_id) {
  vapply(
    strsplit(as.character(taxon_id), "\\|"),
    function(parts) {
      hit <- parts[grepl("^p__", parts)]
      if (length(hit) == 0L) return("Unclassified")
      value <- trimws(sub("^p__", "", hit[[1L]]))
      if (!nzchar(value)) "Unclassified" else value
    },
    character(1)
  )
}

display_phylum <- function(x) {
  dplyr::recode(
    x,
    "Bacillota" = "Bacillota (Firmicutes)",
    .default = x
  )
}

# ============================== 1) Rutas =====================================

PROJECT_DIR <- Sys.getenv("BIOGEOGRAPHY_BASE", unset = ".")
BASE_DIR <- PROJECT_DIR
RESULTS_DIR <- Sys.getenv("BIOGEOGRAPHY_RESULTS", unset = file.path(PROJECT_DIR, "results"))

SCRIPT_09A_TAG <- "05A_TAXONOMIC_CORE_ASSOCIATED"
LATEST_09A_PTR <- file.path(
  RESULTS_DIR, paste0("LATEST_", SCRIPT_09A_TAG, ".txt")
)

SCRIPT_TAG <- "05C_CORE_ASSOCIATED_PHYLA"
RUN_STAMP  <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUT_DIR    <- file.path(RESULTS_DIR, paste0(SCRIPT_TAG, "_", RUN_STAMP))
TABLE_DIR  <- file.path(OUT_DIR, "tables")
PLOT_DIR   <- file.path(OUT_DIR, "plots")
RDS_DIR    <- file.path(OUT_DIR, "rds")

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_DIR,   recursive = TRUE, showWarnings = FALSE)

LATEST_PTR <- file.path(RESULTS_DIR, paste0("LATEST_", SCRIPT_TAG, ".txt"))

# ============================= 2) Parametros =================================

TOP_N_PHYLA <- 5L
MIN_LABEL_PERCENT <- 4
EXPECTED_MIN_REALM_LOCALITIES <- 3L
BAR_LABEL_COLOUR <- "black"
BAR_LABEL_SIZE   <- 4.8

PHYLA_COLOURS <- c(
  "#5B2A86",
  "#1B9E77",
  "#377EB8",
  "#E6AB02",
  "#D95F02"
)
OTHER_COLOUR <- "#BDBDBD"

RUN_09A_ENV <- trimws(Sys.getenv("RUN_09A", unset = ""))

if (!nzchar(RUN_09A_ENV) && !file.exists(LATEST_09A_PTR)) {
  stopf("No se encontro el puntero 09A: %s", LATEST_09A_PTR)
}

RUN_09A_SOURCE <- if (nzchar(RUN_09A_ENV)) "RUN_09A environment variable" else "LATEST pointer"
RUN_09A <- if (nzchar(RUN_09A_ENV)) {
  RUN_09A_ENV
} else {
  trimws(readLines(LATEST_09A_PTR, warn = FALSE)[1])
}
if (!nzchar(RUN_09A) || !dir.exists(RUN_09A)) {
  stopf("El run 09A resuelto mediante %s no existe: %s", RUN_09A_SOURCE, RUN_09A)
}

INPUT_RDS <- file.path(RUN_09A, "rds", "09A_Figure_3a_analysis_objects.rds")
if (!file.exists(INPUT_RDS)) {
  stopf("No se encontro el objeto de analisis 09A: %s", INPUT_RDS)
}

msg("Starting ", SCRIPT_TAG)
msg("Input 09A run: ", RUN_09A)
msg("Input 09A source: ", RUN_09A_SOURCE)
msg("Input RDS: ", INPUT_RDS)
msg("Output: ", OUT_DIR)

# =============================== 3) Carga ====================================

obj09a <- readRDS(INPUT_RDS)

required_objects <- c(
  "locality_counts",
  "locality_metadata",
  "taxon_classification",
  "taxon_id_harmonization_audit",
  "parameters"
)

missing_objects <- setdiff(required_objects, names(obj09a))
if (length(missing_objects) > 0L) {
  stopf(
    "El objeto 09A no contiene: %s. Ejecute primero el 09A armonizado.",
    paste(missing_objects, collapse = ", ")
  )
}

otu_locality <- as.matrix(obj09a$locality_counts)
locality_meta <- as.data.frame(obj09a$locality_metadata)
taxon_classification <- as.data.frame(obj09a$taxon_classification)
harmonization_audit <- as.data.frame(obj09a$taxon_id_harmonization_audit)

parameter_value <- function(parameter_name) {
  hit <- obj09a$parameters$value[obj09a$parameters$parameter == parameter_name]
  if (length(hit) != 1L) {
    stopf("09A debe contener exactamente un valor para el parametro '%s'.", parameter_name)
  }
  hit[[1L]]
}

observed_min_realm_localities <- as.integer(parameter_value(
  "minimum_realm_localities_to_support_core"
))
if (is.na(observed_min_realm_localities) ||
    observed_min_realm_localities != EXPECTED_MIN_REALM_LOCALITIES) {
  stopf(
    "El run 09A usa minimum_realm_localities=%s; Figure 3 requiere %d.",
    observed_min_realm_localities,
    EXPECTED_MIN_REALM_LOCALITIES
  )
}

if (is.null(rownames(otu_locality)) || is.null(colnames(otu_locality))) {
  stopf("locality_counts debe tener rownames y colnames.")
}
if (any(!is.finite(otu_locality)) || any(otu_locality < 0)) {
  stopf("locality_counts contiene valores no finitos o negativos.")
}
if (anyDuplicated(rownames(otu_locality)) || anyDuplicated(colnames(otu_locality))) {
  stopf("locality_counts contiene localidades o taxon_id duplicados.")
}

required_meta <- c("locality_id", "Locality", "REALM")
if (!all(required_meta %in% names(locality_meta))) {
  stopf(
    "Faltan columnas en locality_metadata: %s",
    paste(setdiff(required_meta, names(locality_meta)), collapse = ", ")
  )
}

locality_meta <- locality_meta[
  match(rownames(otu_locality), locality_meta$locality_id),
  , drop = FALSE
]
if (anyNA(locality_meta$locality_id)) {
  stopf("Fallo la alineacion entre locality_counts y locality_metadata.")
}

required_class <- c(
  "taxon_id", "total_reads", "n_localities", "core_associated", "component"
)
if (!all(required_class %in% names(taxon_classification))) {
  stopf(
    "Faltan columnas en taxon_classification: %s",
    paste(setdiff(required_class, names(taxon_classification)), collapse = ", ")
  )
}
if (anyDuplicated(taxon_classification$taxon_id)) {
  stopf("taxon_classification contiene taxon_id duplicados.")
}

n_ids_merged_upstream <- sum(
  harmonization_audit$mapping_status == "Merged into full lineage",
  na.rm = TRUE
)
if (n_ids_merged_upstream == 0L) {
  stopf("09A no reporta identificadores abreviados armonizados; revise el input.")
}

core_ids <- taxon_classification$taxon_id[
  taxon_classification$core_associated %in% TRUE
]
core_ids <- as.character(core_ids)

if (length(core_ids) == 0L) stopf("No se encontraron taxa core-associated.")

missing_core <- setdiff(core_ids, colnames(otu_locality))
if (length(missing_core) > 0L) {
  stopf(
    "Hay %d taxa core-associated ausentes de locality_counts. Primeros: %s",
    length(missing_core), paste(head(missing_core, 20), collapse = ", ")
  )
}

msg("Locality units: ", nrow(otu_locality))
msg("Detected taxonomic features: ", ncol(otu_locality))
msg("Eligible features: ", nrow(taxon_classification))
msg("Core-associated features: ", length(core_ids))
msg("Taxon IDs merged upstream: ", n_ids_merged_upstream)

# ===================== 4) Composicion core por phylum ========================

core_counts <- otu_locality[, core_ids, drop = FALSE]
core_reads_by_locality <- rowSums(core_counts)

if (any(core_reads_by_locality <= 0)) {
  bad <- rownames(core_counts)[core_reads_by_locality <= 0]
  stopf(
    "Hay %d localidades sin lecturas core-associated: %s",
    length(bad), paste(head(bad, 20), collapse = ", ")
  )
}

phylum_by_taxon <- extract_phylum(core_ids)
names(phylum_by_taxon) <- core_ids

phylum_counts <- t(rowsum(
  t(core_counts),
  group = phylum_by_taxon[colnames(core_counts)],
  reorder = FALSE
))

phylum_relative_locality <- sweep(
  phylum_counts,
  1,
  rowSums(phylum_counts),
  "/"
) * 100

feature_counts_by_phylum <- table(phylum_by_taxon)
pooled_reads_by_phylum <- colSums(phylum_counts)

phylum_summary <- tibble(
  phylum = colnames(phylum_counts),
  display_name = display_phylum(phylum),
  n_core_features = as.integer(feature_counts_by_phylum[phylum]),
  n_localities_detected = as.integer(colSums(phylum_counts > 0)),
  core_reads = as.numeric(pooled_reads_by_phylum),
  pooled_read_percent = 100 * core_reads / sum(core_reads),
  locality_mean_percent = as.numeric(colMeans(phylum_relative_locality)),
  locality_median_percent = as.numeric(apply(
    phylum_relative_locality, 2, median
  )),
  locality_sd_percent = as.numeric(apply(
    phylum_relative_locality, 2, sd
  ))
) %>%
  arrange(desc(locality_mean_percent), phylum)

if (abs(sum(phylum_summary$locality_mean_percent) - 100) > 1e-8) {
  stopf("Las proporciones locality-balanced no suman 100%%.")
}
if (abs(sum(phylum_summary$pooled_read_percent) - 100) > 1e-8) {
  stopf("Las proporciones pooled no suman 100%%.")
}

# Confirmar que los conteos recuperan exactamente las lecturas core de 09A.
reported_core_reads <- sum(
  taxon_classification$total_reads[
    taxon_classification$core_associated %in% TRUE
  ]
)
observed_core_reads <- sum(core_counts)

if (!isTRUE(all.equal(
  as.numeric(reported_core_reads),
  as.numeric(observed_core_reads),
  tolerance = 0
))) {
  stopf(
    "Las lecturas core difieren entre 09A (%s) y locality_counts (%s).",
    reported_core_reads, observed_core_reads
  )
}

locality_phylum_table <- as.data.frame(phylum_relative_locality) %>%
  rownames_to_column("locality_id") %>%
  pivot_longer(
    cols = -locality_id,
    names_to = "phylum",
    values_to = "within_locality_core_read_percent"
  ) %>%
  left_join(
    as_tibble(locality_meta) %>% select(locality_id, Locality, REALM),
    by = "locality_id"
  ) %>%
  select(locality_id, Locality, REALM, phylum,
         within_locality_core_read_percent)

# ========================== 5) Datos del panel ================================

top_phyla <- phylum_summary %>%
  slice_head(n = TOP_N_PHYLA) %>%
  pull(phylum)

plot_summary <- phylum_summary %>%
  mutate(
    plot_group = if_else(phylum %in% top_phyla, phylum, "Other phyla")
  ) %>%
  group_by(plot_group) %>%
  summarise(
    n_core_features = sum(n_core_features),
    core_reads = sum(core_reads),
    pooled_read_percent = sum(pooled_read_percent),
    locality_mean_percent = sum(locality_mean_percent),
    .groups = "drop"
  ) %>%
  mutate(
    display_name = if_else(
      plot_group == "Other phyla",
      "Other phyla",
      display_phylum(plot_group)
    )
  )

display_order <- c(
  display_phylum(top_phyla),
  if (any(plot_summary$plot_group == "Other phyla")) "Other phyla"
)

plot_summary <- plot_summary %>%
  mutate(
    display_name = factor(display_name, levels = display_order),
    label = if_else(
      locality_mean_percent >= MIN_LABEL_PERCENT,
      sprintf("%.1f%%", locality_mean_percent),
      ""
    )
  ) %>%
  arrange(display_name)

if (abs(sum(plot_summary$locality_mean_percent) - 100) > 1e-8) {
  stopf("Los grupos mostrados en el panel no suman 100%%.")
}

fill_values <- c(
  setNames(PHYLA_COLOURS[seq_along(top_phyla)], display_phylum(top_phyla)),
  "Other phyla" = OTHER_COLOUR
)
fill_values <- fill_values[display_order]

write_tsv(
  phylum_summary,
  file.path(TABLE_DIR, "11_Figure_3b_full_phylum_composition.tsv")
)
write_tsv(
  plot_summary %>% mutate(display_name = as.character(display_name)),
  file.path(TABLE_DIR, "11_Figure_3b_plot_data.tsv")
)
write_tsv(
  locality_phylum_table,
  file.path(TABLE_DIR, "11_Figure_3b_locality_phylum_composition.tsv")
)

# =============================== 6) Grafico ==================================

p <- ggplot(
  plot_summary,
  aes(
    x = locality_mean_percent,
    y = "Core-associated reads",
    fill = display_name
  )
) +
  geom_col(
    width = 0.50,
    colour = "white",
    linewidth = 0.35,
    position = position_stack(reverse = TRUE)
  ) +
  geom_text(
    aes(label = label),
    colour = BAR_LABEL_COLOUR,
    fontface = "bold",
    size = BAR_LABEL_SIZE,
    position = position_stack(vjust = 0.5, reverse = TRUE)
  ) +
  scale_fill_manual(
    values = fill_values,
    breaks = display_order,
    name = NULL,
    drop = FALSE
  ) +
  scale_x_continuous(
    name = "Mean percentage of core-associated reads",
    limits = c(0, 100),
    breaks = seq(0, 100, 25),
    labels = label_percent(scale = 1, accuracy = 1),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_discrete(name = NULL) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(colour = "black", size = 10),
    axis.text.y = element_text(colour = "black", size = 10.5),
    axis.title.x = element_text(size = 11, margin = margin(t = 7)),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.justification = "center",
    legend.text = element_text(size = 9.5),
    legend.key.width = grid::unit(0.55, "cm"),
    legend.spacing.x = grid::unit(0.16, "cm"),
    plot.margin = margin(8, 12, 4, 8)
  )

png_file <- file.path(
  PLOT_DIR, "Figure_3b_core_associated_phylum_composition.png"
)
pdf_file <- file.path(
  PLOT_DIR, "Figure_3b_core_associated_phylum_composition.pdf"
)
tiff_file <- file.path(
  PLOT_DIR, "Figure_3b_core_associated_phylum_composition.tiff"
)

ggsave(
  png_file, p, width = 7.2, height = 3.2,
  units = "in", dpi = 600, bg = "white"
)
ggsave(
  pdf_file, p, width = 7.2, height = 3.2,
  units = "in", device = cairo_pdf, bg = "white"
)
ggsave(
  tiff_file, p, width = 7.2, height = 3.2,
  units = "in", dpi = 600, compression = "lzw", bg = "white"
)

# ============================= 7) Provenance =================================

parameters <- tibble(
  parameter = c(
    "input_core_classification",
    "composition_unit",
    "primary_composition_metric",
    "sensitivity_composition_metric",
    "taxonomic_rank",
    "top_n_phyla",
    "minimum_segment_label_percent",
    "expected_minimum_realm_localities",
    "bar_label_colour",
    "bar_label_size"
  ),
  value = c(
    "09A harmonized core-associated classification",
    "REALM || Locality",
    "Mean within-locality percentage of core-associated reads",
    "Pooled percentage of core-associated reads",
    "Phylum parsed from canonical lineage identifier",
    TOP_N_PHYLA,
    MIN_LABEL_PERCENT,
    EXPECTED_MIN_REALM_LOCALITIES,
    BAR_LABEL_COLOUR,
    BAR_LABEL_SIZE
  )
)

top_summary_text <- paste(
  sprintf(
    "%s=%.3f%%",
    phylum_summary$display_name[seq_len(min(5L, nrow(phylum_summary)))],
    phylum_summary$locality_mean_percent[seq_len(min(5L, nrow(phylum_summary)))]
  ),
  collapse = " | "
)

run_summary <- tibble(
  item = c(
    "script",
    "input_09A_run",
    "input_09A_source",
    "input_09A_rds",
    "n_locality_units",
    "n_detected_taxa",
    "n_eligible_taxa",
    "n_core_associated_taxa",
    "n_core_associated_reads",
    "n_phyla_in_core",
    "n_taxon_ids_merged_upstream",
    "top_phyla_locality_balanced",
    "output_directory"
  ),
  value = c(
    SCRIPT_TAG,
    RUN_09A,
    RUN_09A_SOURCE,
    INPUT_RDS,
    nrow(otu_locality),
    ncol(otu_locality),
    nrow(taxon_classification),
    length(core_ids),
    observed_core_reads,
    nrow(phylum_summary),
    n_ids_merged_upstream,
    top_summary_text,
    OUT_DIR
  )
)

write_tsv(parameters, file.path(TABLE_DIR, "11_Figure_3b_parameters.tsv"))
write_tsv(run_summary, file.path(TABLE_DIR, "11_run_summary.tsv"))

saveRDS(
  list(
    input_09A_run = RUN_09A,
    parameters = parameters,
    phylum_by_taxon = phylum_by_taxon,
    phylum_counts = phylum_counts,
    phylum_relative_locality = phylum_relative_locality,
    phylum_summary = phylum_summary,
    plot_summary = plot_summary,
    plot = p
  ),
  file.path(RDS_DIR, "11_Figure_3b_analysis_objects.rds")
)

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
writeLines(OUT_DIR, LATEST_PTR)

msg("Core-associated reads: ", observed_core_reads)
msg("Phyla detected in core-associated fraction: ", nrow(phylum_summary))
msg("Top phyla: ", top_summary_text)
msg("Figure: ", png_file)
msg("=== SCRIPT COMPLETED SUCCESSFULLY ===")
