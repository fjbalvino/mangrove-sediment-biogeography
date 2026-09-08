#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# ==============================================================================
# Script 09B: SUPPLEMENTARY_REALM_RESTRICTED_TAXA_SET401.R
# Proyecto: Biogeografia_2026 / MAGENTA
#
# Objetivo:
#   Recalcular la sensibilidad de taxa estrictamente restringidos a un marine
#   realm usando la tabla NO FILTRADA por taxon y las localidades auditadas por
#   el script 09A. Genera una figura suplementaria que separa los taxa
#   restringidos detectados en una sola localidad de aquellos recuperados en
#   dos o mas localidades del mismo realm.
#
# Definicion:
#   Realm-restricted = taxon detectado en uno y solo un marine realm.
#   Esta clasificacion no requiere abundancia minima ni prevalencia minima y se
#   calcula sobre todos los taxa no cero. El resultado filtrado (>=50 lecturas
#   globales y >=5 localidades) se conserva como comparacion tabular.
#
# Dependencia deliberada:
#   Se reutiliza el objeto 09A para garantizar exactamente la misma alineacion
#   por run_accession, definicion REALM || Locality y universo taxonomico.
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

# ============================== 1) Rutas =====================================

BASE_DIR <- Sys.getenv("BIOGEOGRAPHY_BASE", unset = "/home/fjbalvino/Biogeografia_2026")

SCRIPT_09A_TAG <- "09A_FIG3A_TAXONOMIC_CORE_BARS_SET401"
LATEST_09A_PTR <- file.path(
  BASE_DIR, "resultados", paste0("LATEST_", SCRIPT_09A_TAG, ".txt")
)

SCRIPT_TAG <- "09B_SUPPLEMENTARY_REALM_RESTRICTED_TAXA_SET401"
RUN_STAMP  <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUT_DIR    <- file.path(BASE_DIR, "resultados", paste0(SCRIPT_TAG, "_", RUN_STAMP))
TABLE_DIR  <- file.path(OUT_DIR, "tables")
PLOT_DIR   <- file.path(OUT_DIR, "plots")
RDS_DIR    <- file.path(OUT_DIR, "rds")

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_DIR,   recursive = TRUE, showWarnings = FALSE)

LATEST_PTR <- file.path(BASE_DIR, "resultados", paste0("LATEST_", SCRIPT_TAG, ".txt"))

if (!file.exists(LATEST_09A_PTR)) {
  stopf("No se encontro el puntero del script 09A: %s", LATEST_09A_PTR)
}

RUN_09A <- trimws(readLines(LATEST_09A_PTR, warn = FALSE)[1])
if (!nzchar(RUN_09A) || !dir.exists(RUN_09A)) {
  stopf("El run 09A indicado por LATEST no existe: %s", RUN_09A)
}

INPUT_RDS <- file.path(RUN_09A, "rds", "09A_Figure_3a_analysis_objects.rds")
if (!file.exists(INPUT_RDS)) stopf("No se encontro el objeto 09A: %s", INPUT_RDS)

msg("Starting ", SCRIPT_TAG)
msg("Input 09A run: ", RUN_09A)
msg("Input RDS: ", INPUT_RDS)
msg("Output: ", OUT_DIR)

# ============================= 2) Carga ======================================

obj09a <- readRDS(INPUT_RDS)

required_objects <- c(
  "locality_counts",
  "locality_metadata",
  "taxon_classification",
  "parameters"
)

missing_objects <- setdiff(required_objects, names(obj09a))
if (length(missing_objects) > 0L) {
  stopf(
    "El objeto 09A no contiene: %s",
    paste(missing_objects, collapse = ", ")
  )
}

otu_locality <- as.matrix(obj09a$locality_counts)
locality_meta <- as.data.frame(obj09a$locality_metadata)
taxon_primary <- as.data.frame(obj09a$taxon_classification)

if (is.null(rownames(otu_locality)) || is.null(colnames(otu_locality))) {
  stopf("locality_counts debe tener rownames y colnames.")
}
if (any(!is.finite(otu_locality)) || any(otu_locality < 0)) {
  stopf("locality_counts contiene valores no finitos o negativos.")
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
if (anyNA(locality_meta$REALM) || any(!nzchar(locality_meta$REALM))) {
  stopf("Hay localidades sin REALM.")
}

taxa_nonzero <- colSums(otu_locality) > 0
otu_locality <- otu_locality[, taxa_nonzero, drop = FALSE]
presence_locality <- otu_locality > 0
taxon_ids <- colnames(otu_locality)

if (anyDuplicated(taxon_ids)) stopf("Hay taxon_id duplicados.")

# Orden geografico usado en las figuras principales, de oeste a este.
realm_order_canonical <- c(
  "Temperate Northern Atlantic",
  "Tropical Eastern Pacific",
  "Tropical Atlantic",
  "Indo-Pacific",
  "Western Indo-Pacific",
  "Central Indo-Pacific"
)

realm_levels_observed <- unique(as.character(locality_meta$REALM))
realm_levels <- c(
  realm_order_canonical[realm_order_canonical %in% realm_levels_observed],
  sort(setdiff(realm_levels_observed, realm_order_canonical))
)

# ================== 3) Presencia por marine realm ============================

realm_presence <- matrix(
  FALSE,
  nrow = length(realm_levels),
  ncol = length(taxon_ids),
  dimnames = list(realm_levels, taxon_ids)
)

for (realm_i in realm_levels) {
  idx <- which(locality_meta$REALM == realm_i)
  realm_presence[realm_i, ] <- colSums(
    presence_locality[idx, , drop = FALSE]
  ) > 0
}

n_realms_detected <- colSums(realm_presence)
n_localities_detected <- colSums(presence_locality)
global_reads <- colSums(otu_locality)

is_realm_restricted <- n_realms_detected == 1L

restricted_realm <- rep(NA_character_, length(taxon_ids))
if (any(is_realm_restricted)) {
  restricted_realm[is_realm_restricted] <- apply(
    realm_presence[, is_realm_restricted, drop = FALSE],
    2,
    function(x) rownames(realm_presence)[which(x)]
  )
}

eligible_taxa <- as.character(taxon_primary$taxon_id)

taxon_audit <- tibble(
  taxon_id = taxon_ids,
  total_reads = as.numeric(global_reads),
  n_localities_detected = as.integer(n_localities_detected),
  n_realms_detected = as.integer(n_realms_detected),
  realm_restricted_unfiltered = is_realm_restricted,
  restricted_realm = restricted_realm,
  restricted_locality_class = case_when(
    !is_realm_restricted ~ NA_character_,
    n_localities_detected == 1L ~ "One locality",
    n_localities_detected >= 2L ~ "Multiple localities"
  ),
  passes_primary_taxon_filter = taxon_id %in% eligible_taxa,
  realm_restricted_after_primary_filter =
    is_realm_restricted & taxon_id %in% eligible_taxa
)

if (any(taxon_audit$n_realms_detected < 1L)) {
  stopf("Hay taxa no cero sin presencia en ningun realm.")
}

# ======================== 4) Resumen por realm ================================

realm_summary <- bind_rows(lapply(realm_levels, function(realm_i) {
  loc_idx <- which(locality_meta$REALM == realm_i)
  taxa_in_realm <- realm_presence[realm_i, ]
  restricted_ids <- taxon_audit$taxon_id[
    taxon_audit$realm_restricted_unfiltered &
      taxon_audit$restricted_realm == realm_i
  ]
  restricted_filtered_ids <- taxon_audit$taxon_id[
    taxon_audit$realm_restricted_after_primary_filter &
      taxon_audit$restricted_realm == realm_i
  ]

  n_total_taxa <- sum(taxa_in_realm)
  n_restricted <- length(restricted_ids)
  n_restricted_filtered <- length(restricted_filtered_ids)
  realm_reads <- sum(otu_locality[loc_idx, , drop = FALSE])
  restricted_reads <- if (n_restricted == 0L) {
    0
  } else {
    sum(otu_locality[loc_idx, restricted_ids, drop = FALSE])
  }

  tibble(
    REALM = realm_i,
    n_localities = length(loc_idx),
    qualifies_to_support_core = length(loc_idx) >= 3L,
    n_total_taxa_detected = n_total_taxa,
    n_realm_restricted_unfiltered = n_restricted,
    n_restricted_one_locality = sum(
      taxon_audit$restricted_realm == realm_i &
        taxon_audit$restricted_locality_class == "One locality",
      na.rm = TRUE
    ),
    n_restricted_multiple_localities = sum(
      taxon_audit$restricted_realm == realm_i &
        taxon_audit$restricted_locality_class == "Multiple localities",
      na.rm = TRUE
    ),
    n_realm_restricted_after_primary_filter = n_restricted_filtered,
    realm_restricted_percent_of_realm_richness = ifelse(
      n_total_taxa > 0,
      100 * n_restricted / n_total_taxa,
      NA_real_
    ),
    total_reads_in_realm = realm_reads,
    realm_restricted_reads = restricted_reads,
    realm_restricted_percent_of_realm_reads = ifelse(
      realm_reads > 0,
      100 * restricted_reads / realm_reads,
      NA_real_
    )
  )
}))

if (any(
  realm_summary$n_restricted_one_locality +
    realm_summary$n_restricted_multiple_localities !=
    realm_summary$n_realm_restricted_unfiltered
)) {
  stopf("No reconcilian las categorias de ocupacion local.")
}

if (any(
  realm_summary$n_realm_restricted_after_primary_filter >
    realm_summary$n_realm_restricted_unfiltered
)) {
  stopf("El resumen filtrado excede el resumen no filtrado.")
}

# ============================== 5) Grafico ===================================

plot_df <- realm_summary %>%
  select(
    REALM,
    n_restricted_one_locality,
    n_restricted_multiple_localities
  ) %>%
  pivot_longer(
    cols = starts_with("n_restricted_"),
    names_to = "occupancy_code",
    values_to = "n_taxa"
  ) %>%
  mutate(
    occupancy = recode(
      occupancy_code,
      n_restricted_one_locality = "One locality",
      n_restricted_multiple_localities = "Multiple localities"
    ),
    occupancy = factor(
      occupancy,
      levels = c("Multiple localities", "One locality")
    ),
    REALM = factor(REALM, levels = rev(realm_levels))
  )

label_df <- realm_summary %>%
  mutate(
    REALM = factor(REALM, levels = rev(realm_levels)),
    label = if_else(
      n_realm_restricted_unfiltered == 0L,
      "0",
      sprintf(
        "%s / %s (%.3f%%)",
        comma(n_realm_restricted_unfiltered),
        comma(n_total_taxa_detected),
        realm_restricted_percent_of_realm_richness
      )
    )
  )

p <- ggplot(plot_df, aes(x = n_taxa, y = REALM, fill = occupancy)) +
  geom_col(
    width = 0.68,
    colour = "white",
    linewidth = 0.35,
    position = position_stack(reverse = TRUE)
  ) +
  geom_text(
    data = label_df,
    aes(
      x = n_realm_restricted_unfiltered,
      y = REALM,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = -0.12,
    size = 3.5,
    colour = "grey20"
  ) +
  scale_fill_manual(
    values = c(
      "Multiple localities" = "#542788",
      "One locality" = "#9E9AC8"
    ),
    breaks = c("Multiple localities", "One locality"),
    name = NULL,
    drop = FALSE
  ) +
  scale_x_continuous(
    name = "Realm-restricted taxonomic features",
    breaks = breaks_pretty(n = 5),
    expand = expansion(mult = c(0, 0.32))
  ) +
  scale_y_discrete(name = NULL) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y = element_text(colour = "black", size = 10.5),
    axis.text.x = element_text(colour = "black", size = 10),
    axis.title.x = element_text(size = 11, margin = margin(t = 7)),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.justification = "center",
    legend.text = element_text(size = 10),
    legend.key.width = grid::unit(0.55, "cm"),
    legend.spacing.x = grid::unit(0.20, "cm"),
    plot.margin = margin(8, 46, 5, 8)
  ) +
  coord_cartesian(clip = "off")

png_file <- file.path(PLOT_DIR, "Supplementary_realm_restricted_taxa.png")
pdf_file <- file.path(PLOT_DIR, "Supplementary_realm_restricted_taxa.pdf")
tiff_file <- file.path(PLOT_DIR, "Supplementary_realm_restricted_taxa.tiff")

ggsave(png_file, p, width = 8.4, height = 5.1, units = "in", dpi = 600,
       bg = "white")
ggsave(pdf_file, p, width = 8.4, height = 5.1, units = "in",
       device = cairo_pdf, bg = "white")
ggsave(tiff_file, p, width = 8.4, height = 5.1, units = "in", dpi = 600,
       compression = "lzw", bg = "white")

# ============================= 6) Exportacion ================================

write_tsv(
  taxon_audit,
  file.path(TABLE_DIR, "09B_realm_restricted_taxon_audit_unfiltered.tsv")
)
write_tsv(
  realm_summary,
  file.path(TABLE_DIR, "09B_realm_restricted_summary.tsv")
)

run_summary <- tibble(
  item = c(
    "script",
    "input_09A_run",
    "input_09A_rds",
    "n_localities",
    "n_realms",
    "n_taxa_unfiltered_nonzero",
    "n_taxa_passing_primary_filter",
    "n_realm_restricted_unfiltered",
    "n_realm_restricted_after_primary_filter",
    "output_directory"
  ),
  value = c(
    SCRIPT_TAG,
    RUN_09A,
    INPUT_RDS,
    nrow(otu_locality),
    length(realm_levels),
    ncol(otu_locality),
    length(eligible_taxa),
    sum(taxon_audit$realm_restricted_unfiltered),
    sum(taxon_audit$realm_restricted_after_primary_filter),
    OUT_DIR
  )
)

write_tsv(run_summary, file.path(TABLE_DIR, "09B_run_summary.tsv"))

saveRDS(
  list(
    input_09A_run = RUN_09A,
    locality_metadata = locality_meta,
    realm_presence = realm_presence,
    taxon_audit = taxon_audit,
    realm_summary = realm_summary,
    plot = p
  ),
  file.path(RDS_DIR, "09B_realm_restricted_analysis_objects.rds")
)

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
writeLines(OUT_DIR, LATEST_PTR)

msg("Realm-restricted taxa, unfiltered: ",
    sum(taxon_audit$realm_restricted_unfiltered))
msg("Realm-restricted taxa after primary filter: ",
    sum(taxon_audit$realm_restricted_after_primary_filter))
msg("Figure: ", png_file)
msg("=== SCRIPT COMPLETED SUCCESSFULLY ===")
