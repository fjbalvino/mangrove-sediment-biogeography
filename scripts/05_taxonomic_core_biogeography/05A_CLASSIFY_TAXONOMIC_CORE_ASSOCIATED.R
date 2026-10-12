#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# ==============================================================================
# 05A_CLASSIFY_TAXONOMIC_CORE_ASSOCIATED.R
# Classify realm-supported taxonomic cores and plot their richness/read shares.
#
# Inputs: data/phyloseq_taxonomic.rds (phyloseq count object) and
# data/metadata_taxonomic.csv, joined by run_accession. Metadata must contain
# REALM and core_locality (or Locality when core_locality is absent). Profiles
# without a resolved locality are retained in the source cohort but excluded
# from locality-based core calculations and listed in the sampling audit.
#
# Counts are pooled within each REALM || Locality unit. Eligibility requires
# >=50 total reads and detection in >=5 such units. Core-associated taxa occur
# in every sampled locality of at least one realm with >=3 localities and have
# mean relative abundance >=0.0001 (a proportion equivalent to 0.01%) in that realm.
# The mean gives each locality equal weight, using proportions of all retained
# counts within that locality before the global eligibility filter.
#
# A core-associated taxon need not occur in every realm. The read-share panel
# pools reads across the eligible taxa; it is not an equal-locality average.
# Abbreviated taxon labels are merged with full lineages only when the terminal
# name has an unambiguous one-to-one match; counts are summed before filtering.
# ==============================================================================

options(stringsAsFactors = FALSE, scipen = 999, bitmapType = "cairo")
Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
set.seed(42)

suppressPackageStartupMessages({
  library(phyloseq)
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

collapse_unique <- function(x) {
  x <- trimws(as.character(x))
  x <- sort(unique(x[!is.na(x) & nzchar(x)]))
  if (length(x) == 0L) NA_character_ else paste(x, collapse = " | ")
}

terminal_taxon_label <- function(x) {
  x <- trimws(as.character(x))
  terminal <- sub("^.*\\|", "", x)
  trimws(sub("^[[:alpha:]]__", "", terminal))
}

# ============================== 1) Rutas =====================================

PROJECT_DIR <- Sys.getenv("BIOGEOGRAPHY_BASE", unset = ".")
BASE_DIR <- PROJECT_DIR
RESULTS_DIR <- Sys.getenv("BIOGEOGRAPHY_RESULTS", unset = file.path(PROJECT_DIR, "results"))

# Taxonomic count object for the core-analysis sample manifest.
PS_RDS <- Sys.getenv("MAGENTA_CORE_PHYLOSEQ", unset = file.path(PROJECT_DIR, "data", "phyloseq_taxonomic.rds"))
SAMPLE_MANIFEST <- Sys.getenv("CORE_SAMPLE_MANIFEST", unset = file.path(PROJECT_DIR, "config", "samples", "taxonomic_accessions.txt"))
if (!file.exists(SAMPLE_MANIFEST)) stopf("Core sample manifest not found: %s", SAMPLE_MANIFEST)
EXPECTED_ACCESSIONS <- trimws(readLines(SAMPLE_MANIFEST, warn = FALSE))
EXPECTED_ACCESSIONS <- EXPECTED_ACCESSIONS[nzchar(EXPECTED_ACCESSIONS) & !startsWith(EXPECTED_ACCESSIONS, "#")]
if (!length(EXPECTED_ACCESSIONS) || anyDuplicated(EXPECTED_ACCESSIONS)) stopf("Core sample manifest must contain unique nonempty identifiers.")

# Sample metadata; run_accession is the unique sample key.
METADATA_CSV <- Sys.getenv(
  "MAGENTA_METADATA_CSV",
  unset = file.path(PROJECT_DIR, "data", "metadata_taxonomic.csv")
)

SCRIPT_TAG <- "05A_TAXONOMIC_CORE_ASSOCIATED"
OUT_DIR    <- file.path(RESULTS_DIR, SCRIPT_TAG)
TABLE_DIR  <- file.path(OUT_DIR, "tables")
PLOT_DIR   <- file.path(OUT_DIR, "plots")
RDS_DIR    <- file.path(OUT_DIR, "rds")

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_DIR,   recursive = TRUE, showWarnings = FALSE)

LATEST_PTR <- file.path(RESULTS_DIR, paste0("LATEST_", SCRIPT_TAG, ".txt"))

# ============================ 2) Parametros ==================================

MIN_GLOBAL_READS      <- 50L
MIN_GLOBAL_LOCALITIES <- 5L
MIN_REALM_LOCALITIES  <- 3L
CORE_REALM_PREVALENCE <- 1.0
CORE_MEAN_REL_ABUND   <- 0.0001  # Proportion: 0.01%, not 0.0001%.
EXPECTED_N_SAMPLES    <- length(EXPECTED_ACCESSIONS)

CORE_COLOUR    <- "#B2182B"
NONCORE_COLOUR <- "#737373"
BAR_LABEL_COLOUR <- "black"
BAR_LABEL_SIZE   <- 4.8

msg("Starting ", SCRIPT_TAG)
msg("Phyloseq: ", PS_RDS)
msg("Metadata: ", METADATA_CSV)
msg("Output: ", OUT_DIR)

if (!file.exists(PS_RDS)) stopf("No se encontro el phyloseq: %s", PS_RDS)
if (!file.exists(METADATA_CSV)) stopf("No se encontro la metadata: %s", METADATA_CSV)

# ========================= 3) Carga y alineacion =============================

ps0 <- readRDS(PS_RDS)
if (!inherits(ps0, "phyloseq")) stopf("El objeto no es phyloseq.")

ps_sd <- as(sample_data(ps0), "data.frame") %>%
  rownames_to_column("phyloseq_sample_name")

ps_key_col <- resolve_col(
  ps_sd,
  c("run_accession", "Run_accession", "RUN_ACCESSION"),
  "run_accession dentro del phyloseq",
  required = FALSE
)

if (is.na(ps_key_col)) {
  ps_run <- as.character(sample_names(ps0))
  msg("run_accession no esta en sample_data; se usaran sample_names(ps).")
} else {
  ps_run <- trimws(as.character(ps_sd[[ps_key_col]]))
}

if (anyNA(ps_run) || any(!nzchar(ps_run))) {
  stopf("Hay run_accession vacios en el phyloseq.")
}
if (anyDuplicated(ps_run)) {
  dup <- unique(ps_run[duplicated(ps_run)])
  stopf("run_accession duplicados en phyloseq: %s", paste(head(dup, 20), collapse = ", "))
}

# Canonizar sample_names con la llave primaria.
sample_names(ps0) <- ps_run

meta_raw <- read_csv(METADATA_CSV, show_col_types = FALSE, progress = FALSE)

meta_key_col <- resolve_col(meta_raw, c("run_accession"), "run_accession")
locality_col <- resolve_col(meta_raw, c("core_locality", "Locality", "locality"), "resolved locality")
realm_col    <- resolve_col(meta_raw, c("REALM", "Realm", "realm"), "REALM")
eco_col      <- resolve_col(
  meta_raw,
  c("EcologicalClass", "Ecological_Class", "ecological_class"),
  "EcologicalClass",
  required = FALSE
)

meta_needed <- meta_raw %>%
  transmute(
    run_accession = trimws(as.character(.data[[meta_key_col]])),
    Locality      = trimws(as.character(.data[[locality_col]])),
    REALM         = trimws(as.character(.data[[realm_col]])),
    EcologicalClass = if (is.na(eco_col)) {
      NA_character_
    } else {
      trimws(as.character(.data[[eco_col]]))
    }
  ) %>%
  filter(!is.na(run_accession), nzchar(run_accession))

# Las filas duplicadas de metadata se permiten solo si Locality/REALM coinciden.
meta_audit <- meta_needed %>%
  group_by(run_accession) %>%
  summarise(
    n_rows = n(),
    n_locality = n_distinct(Locality[!is.na(Locality) & nzchar(Locality)]),
    n_realm = n_distinct(REALM[!is.na(REALM) & nzchar(REALM)]),
    n_ecological_classes = n_distinct(EcologicalClass[!is.na(EcologicalClass) & nzchar(EcologicalClass)]),
    localities = collapse_unique(Locality),
    realms = collapse_unique(REALM),
    ecological_classes = collapse_unique(EcologicalClass),
    .groups = "drop"
  )

write_tsv(meta_audit, file.path(TABLE_DIR, "09A_metadata_run_accession_audit.tsv"))

bad_meta <- meta_audit %>% filter(n_locality > 1L | n_realm > 1L | n_ecological_classes > 1L)
if (nrow(bad_meta) > 0L) {
  stopf(
    paste0("Hay %d run_accession con Locality o REALM incompatibles. ",
           "Revise 09A_metadata_run_accession_audit.tsv."),
    nrow(bad_meta)
  )
}

meta_one <- meta_needed %>%
  group_by(run_accession) %>%
  summarise(
    Locality = collapse_unique(Locality),
    REALM = collapse_unique(REALM),
    EcologicalClass = collapse_unique(EcologicalClass),
    .groups = "drop"
  )

idx <- match(sample_names(ps0), meta_one$run_accession)
missing_meta <- sample_names(ps0)[is.na(idx)]
if (length(missing_meta) > 0L) {
  stopf(
    "Faltan %d run_accession del phyloseq en la metadata. Primeros: %s",
    length(missing_meta), paste(head(missing_meta, 20), collapse = ", ")
  )
}

meta_aligned <- meta_one[idx, , drop = FALSE]
if (!identical(meta_aligned$run_accession, sample_names(ps0))) {
  stopf("Fallo la alineacion por run_accession.")
}

if (!setequal(meta_aligned$run_accession, EXPECTED_ACCESSIONS)) {
  stopf("The input samples do not match the core sample manifest: %s", SAMPLE_MANIFEST)
}

if (nrow(meta_aligned) != EXPECTED_N_SAMPLES) {
  stopf(
    "Figure 3 requiere %d metagenomas; el objeto alineado contiene %d.",
    EXPECTED_N_SAMPLES,
    nrow(meta_aligned)
  )
}

# Unresolved geographic regions are not independent sampled localities.
locality_resolved <- !is.na(meta_aligned$Locality) & nzchar(meta_aligned$Locality)
locality_inclusion <- meta_aligned %>%
  mutate(included_in_locality_analysis = locality_resolved,
         reason = if_else(locality_resolved, "Locality identifier available",
                          "Sampling locality unresolved"))
write_tsv(locality_inclusion,
          file.path(TABLE_DIR, "09A_locality_analysis_sample_inclusion.tsv"))
N_INPUT_SAMPLES <- nrow(meta_aligned)
N_LOCALITY_UNRESOLVED <- sum(!locality_resolved)
if (N_LOCALITY_UNRESOLVED > 0L) {
  msg("Profiles excluded from locality-based analysis: ", N_LOCALITY_UNRESOLVED)
  meta_aligned <- meta_aligned[locality_resolved, , drop = FALSE]
  ps0 <- prune_samples(meta_aligned$run_accession, ps0)
  stopifnot(identical(sample_names(ps0), meta_aligned$run_accession))
}
if (!nrow(meta_aligned)) stopf("No profiles have a resolved sampling locality.")
if (anyNA(meta_aligned$REALM) || any(!nzchar(meta_aligned$REALM))) {
  stopf("Hay muestras sin REALM despues del match por run_accession.")
}

# The input object must contain mangrove sediment samples only.
if (!all(is.na(meta_aligned$EcologicalClass))) {
  is_mangrove_sediment <- tolower(meta_aligned$EcologicalClass) == "mangrove sediment"
  if (any(!is_mangrove_sediment, na.rm = TRUE)) {
    bad_ids <- meta_aligned$run_accession[!is_mangrove_sediment &
                                           !is.na(is_mangrove_sediment)]
    stopf(
      paste0("El objeto contiene %d muestras fuera de EcologicalClass == ",
             "'Mangrove sediment'. Primeras: %s"),
      length(bad_ids), paste(head(bad_ids, 20), collapse = ", ")
    )
  }
}

write_tsv(meta_aligned, file.path(TABLE_DIR, "09A_metadata_aligned_by_run_accession.tsv"))

# ======================= 4) Matriz por localidad =============================

otu_sample <- as(otu_table(ps0), "matrix")
if (taxa_are_rows(ps0)) otu_sample <- t(otu_sample)
otu_sample <- otu_sample[sample_names(ps0), , drop = FALSE]

if (any(!is.finite(otu_sample)) || any(otu_sample < 0)) {
  stopf("La tabla taxonomica contiene valores no finitos o negativos.")
}

# Taxon labels may contain full lineages (k__...|s__Species) or terminal names.
# Merge a terminal name with its full lineage only when the mapping is unique.
# Ambiguous mappings stop the analysis to prevent conflation of distinct taxa.
raw_taxon_id <- colnames(otu_sample)
if (is.null(raw_taxon_id) || anyNA(raw_taxon_id) || any(!nzchar(raw_taxon_id))) {
  stopf("La tabla taxonomica contiene taxon_id vacios.")
}
if (anyDuplicated(raw_taxon_id)) {
  stopf("La tabla taxonomica contiene taxon_id exactamente duplicados.")
}

taxon_terminal <- terminal_taxon_label(raw_taxon_id)
is_full_lineage <- grepl("|", raw_taxon_id, fixed = TRUE)
lineage_by_terminal <- split(
  raw_taxon_id[is_full_lineage],
  taxon_terminal[is_full_lineage]
)

canonical_taxon_id <- raw_taxon_id
mapping_status <- ifelse(is_full_lineage, "Full lineage retained", "No lineage match")
ambiguous_terminal <- character(0)

for (i in which(!is_full_lineage)) {
  candidates <- lineage_by_terminal[[taxon_terminal[[i]]]]
  if (is.null(candidates) || length(candidates) == 0L) next
  if (length(candidates) > 1L) {
    ambiguous_terminal <- c(ambiguous_terminal, taxon_terminal[[i]])
    mapping_status[[i]] <- "Ambiguous lineage match"
  } else {
    canonical_taxon_id[[i]] <- candidates[[1L]]
    mapping_status[[i]] <- "Merged into full lineage"
  }
}

taxon_id_harmonization_audit <- tibble(
  original_taxon_id = raw_taxon_id,
  terminal_taxon_label = taxon_terminal,
  original_id_format = if_else(
    is_full_lineage, "Full lineage", "Abbreviated label"
  ),
  canonical_taxon_id = canonical_taxon_id,
  mapping_status = mapping_status,
  reads_before_harmonization = as.numeric(colSums(otu_sample))
)

write_tsv(
  taxon_id_harmonization_audit,
  file.path(TABLE_DIR, "09A_taxon_id_harmonization_audit.tsv")
)

if (length(ambiguous_terminal) > 0L) {
  stopf(
    paste0(
      "Hay %d etiquetas abreviadas con mas de un linaje candidato. ",
      "Revise 09A_taxon_id_harmonization_audit.tsv. Primeras: %s"
    ),
    length(unique(ambiguous_terminal)),
    paste(head(unique(ambiguous_terminal), 20), collapse = ", ")
  )
}

n_ids_merged <- sum(canonical_taxon_id != raw_taxon_id)
reads_in_ids_merged <- sum(
  colSums(otu_sample)[canonical_taxon_id != raw_taxon_id]
)

if (n_ids_merged > 0L) {
  otu_sample <- t(rowsum(
    t(otu_sample),
    group = canonical_taxon_id,
    reorder = FALSE
  ))
}

if (anyDuplicated(colnames(otu_sample))) {
  stopf("La armonizacion produjo taxon_id canonicos duplicados.")
}

msg("Taxon IDs merged into full lineages: ", n_ids_merged)
msg("Reads recovered from abbreviated duplicate IDs: ", reads_in_ids_merged)
msg("Canonical taxonomic features: ", ncol(otu_sample))

# Auditar etiquetas de localidad asociadas con mas de un realm. El analisis
# conserva REALM || Locality como unidad para evitar fusionarlas en silencio,
# pero estas asignaciones deben resolverse en la metadata antes del manuscrito.
locality_realm_audit <- meta_aligned %>%
  group_by(Locality) %>%
  summarise(
    n_samples = n(),
    n_realms = n_distinct(REALM),
    realms = collapse_unique(REALM),
    .groups = "drop"
  ) %>%
  arrange(desc(n_realms), Locality)

write_tsv(
  locality_realm_audit,
  file.path(TABLE_DIR, "09A_locality_realm_assignment_audit.tsv")
)

n_locality_realm_conflicts <- sum(locality_realm_audit$n_realms > 1L)
if (n_locality_realm_conflicts > 0L) {
  msg(
    "WARNING: locality labels assigned to multiple realms: ",
    paste(
      locality_realm_audit$Locality[locality_realm_audit$n_realms > 1L],
      collapse = "; "
    )
  )
}

# Prefijo REALM evita colisionar etiquetas de localidad repetidas entre realms.
locality_id <- paste(meta_aligned$REALM, meta_aligned$Locality, sep = " || ")
otu_locality <- rowsum(otu_sample, group = locality_id, reorder = FALSE)

locality_meta <- tibble(
  locality_id = locality_id,
  Locality = meta_aligned$Locality,
  REALM = meta_aligned$REALM
) %>%
  distinct()

if (anyDuplicated(locality_meta$locality_id)) {
  stopf("locality_id no es unica despues de combinar REALM y Locality.")
}

locality_meta <- locality_meta[match(rownames(otu_locality), locality_meta$locality_id), ]
if (anyNA(locality_meta$locality_id)) stopf("Fallo la alineacion de localidades.")

taxa_nonzero <- colSums(otu_locality) > 0
otu_locality <- otu_locality[, taxa_nonzero, drop = FALSE]

locality_totals <- rowSums(otu_locality)
if (any(locality_totals <= 0)) stopf("Hay localidades con cero lecturas clasificadas.")

# Denominador de abundancia: comunidad completa antes del prefiltro de taxa.
rel_locality_all <- sweep(otu_locality, 1, locality_totals, "/")

global_reads <- colSums(otu_locality)
n_localities <- colSums(otu_locality > 0)

eligible <- global_reads >= MIN_GLOBAL_READS &
  n_localities >= MIN_GLOBAL_LOCALITIES

if (!any(eligible)) stopf("Ningun taxon paso el prefiltro global.")

otu_eligible <- otu_locality[, eligible, drop = FALSE]
rel_eligible <- rel_locality_all[, eligible, drop = FALSE]
eligible_ids <- colnames(otu_eligible)

msg("Samples aligned: ", nrow(meta_aligned))
msg("Realm-locality units: ", nrow(otu_locality))
msg("Detected non-zero taxa: ", ncol(otu_locality))
msg("Eligible taxa (>=50 reads and >=5 localities): ", ncol(otu_eligible))

# ===================== 5) Clasificacion core-associated ======================

realm_levels <- sort(unique(locality_meta$REALM))

realm_counts <- locality_meta %>%
  count(REALM, name = "n_localities") %>%
  arrange(desc(n_localities), REALM) %>%
  mutate(
    qualifies_to_support_core = n_localities >= MIN_REALM_LOCALITIES
  )

qualifying_realms <- realm_counts %>%
  filter(qualifies_to_support_core) %>%
  pull(REALM)

if (length(qualifying_realms) == 0L) {
  stopf(
    "Ningun marine realm tiene >= %d localidades.",
    MIN_REALM_LOCALITIES
  )
}

msg(
  "Realms qualifying to support core (>=",
  MIN_REALM_LOCALITIES,
  " localities): ",
  paste(qualifying_realms, collapse = "; ")
)

realm_prevalence <- matrix(
  NA_real_, nrow = length(realm_levels), ncol = length(eligible_ids),
  dimnames = list(realm_levels, eligible_ids)
)
realm_mean_rel_abund <- realm_prevalence

for (realm_i in realm_levels) {
  loc_idx <- which(locality_meta$REALM == realm_i)
  realm_prevalence[realm_i, ] <- colMeans(
    otu_eligible[loc_idx, , drop = FALSE] > 0
  )
  realm_mean_rel_abund[realm_i, ] <- colMeans(
    rel_eligible[loc_idx, , drop = FALSE]
  )
}

core_support_unrestricted <- realm_prevalence >= CORE_REALM_PREVALENCE &
  realm_mean_rel_abund >= CORE_MEAN_REL_ABUND

# Un realm con menos de MIN_REALM_LOCALITIES no aporta evidencia de
# persistencia regional y, por tanto, no puede sostener por si solo el core.
core_support <- core_support_unrestricted
core_support[!rownames(core_support) %in% qualifying_realms, ] <- FALSE

is_core_unrestricted <- colSums(core_support_unrestricted) > 0L
is_core <- colSums(core_support) > 0L

supporting_realms <- vapply(
  seq_along(eligible_ids),
  function(j) paste(rownames(core_support)[core_support[, j]], collapse = " | "),
  character(1)
)
supporting_realms[!nzchar(supporting_realms)] <- NA_character_

taxon_classification <- tibble(
  taxon_id = eligible_ids,
  total_reads = as.numeric(global_reads[eligible_ids]),
  n_localities = as.integer(n_localities[eligible_ids]),
  core_associated = is_core,
  core_associated_without_minimum_realm_size = is_core_unrestricted,
  reclassified_by_minimum_realm_size = is_core_unrestricted & !is_core,
  component = if_else(is_core, "Core-associated", "Non-core"),
  supporting_realms = supporting_realms,
  max_realm_prevalence = apply(realm_prevalence, 2, max, na.rm = TRUE),
  max_realm_mean_relative_abundance = apply(
    realm_mean_rel_abund, 2, max, na.rm = TRUE
  )
)

if (!any(taxon_classification$core_associated)) {
  stopf("No se identificaron taxa core-associated; revise metadata y umbrales.")
}

# Exportar matrices realm x taxon en formato largo para auditoria.
realm_taxon_audit <- expand_grid(
  REALM = realm_levels,
  taxon_id = eligible_ids
) %>%
  mutate(
    prevalence = as.vector(t(realm_prevalence)),
    mean_relative_abundance = as.vector(t(realm_mean_rel_abund)),
    n_realm_localities = realm_counts$n_localities[
      match(REALM, realm_counts$REALM)
    ],
    realm_qualifies_to_support_core =
      n_realm_localities >= MIN_REALM_LOCALITIES,
    meets_taxon_thresholds = prevalence >= CORE_REALM_PREVALENCE &
      mean_relative_abundance >= CORE_MEAN_REL_ABUND,
    supports_core = realm_qualifies_to_support_core & meets_taxon_thresholds
  )

# Sensibilidad al numero minimo de localidades requerido para que un realm
# pueda sostener la clasificacion. El umbral primario es 3; se conservan 1-5
# para mostrar cuanto depende el resultado de realms escasamente muestreados.
core_minrealm1 <- is_core_unrestricted

minimum_realm_localities_sensitivity <- bind_rows(lapply(
  1:5,
  function(min_realm_localities_i) {
    qualifying_i <- realm_counts %>%
      filter(n_localities >= min_realm_localities_i) %>%
      pull(REALM)

    support_i <- core_support_unrestricted
    support_i[!rownames(support_i) %in% qualifying_i, ] <- FALSE
    core_i <- colSums(support_i) > 0L

    core_reads_i <- sum(global_reads[eligible_ids][core_i])
    total_reads_i <- sum(global_reads[eligible_ids])

    tibble(
      minimum_realm_localities = min_realm_localities_i,
      is_primary = min_realm_localities_i == MIN_REALM_LOCALITIES,
      n_qualifying_realms = length(qualifying_i),
      qualifying_realms = paste(qualifying_i, collapse = " | "),
      n_core_associated = sum(core_i),
      n_non_core = sum(!core_i),
      core_richness_percent = 100 * mean(core_i),
      core_classified_reads = core_reads_i,
      eligible_classified_reads = total_reads_i,
      core_classified_reads_percent = 100 * core_reads_i / total_reads_i,
      taxa_reclassified_vs_minrealm1 = sum(core_minrealm1 & !core_i),
      reads_reclassified_vs_minrealm1 = sum(
        global_reads[eligible_ids][core_minrealm1 & !core_i]
      )
    )
  }
))

write_tsv(
  taxon_classification,
  file.path(TABLE_DIR, "09A_taxonomic_core_classification.tsv")
)
write_tsv(
  realm_taxon_audit,
  file.path(TABLE_DIR, "09A_realm_taxon_core_audit.tsv")
)
write_tsv(
  minimum_realm_localities_sensitivity,
  file.path(TABLE_DIR, "09A_minimum_realm_localities_sensitivity.tsv")
)

# ===================== 6) Resumen para la Figura 3a ==========================

# Analisis primario: denominadores restringidos al universo de taxa elegibles,
# exactamente como establece el metodo de clasificacion.
figure_summary <- taxon_classification %>%
  group_by(component) %>%
  summarise(
    richness = n(),
    classified_reads = sum(total_reads),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(richness, classified_reads),
    names_to = "metric_code",
    values_to = "value"
  ) %>%
  group_by(metric_code) %>%
  mutate(
    denominator = sum(value),
    percentage = 100 * value / denominator
  ) %>%
  ungroup() %>%
  mutate(
    metric = recode(
      metric_code,
      richness = "Taxonomic richness",
      classified_reads = "Reads in eligible taxa"
    ),
    component = factor(component, levels = c("Core-associated", "Non-core")),
    metric = factor(
      metric,
      levels = c("Reads in eligible taxa", "Taxonomic richness")
    )
  )

# Preserve the distinction between eligible non-core features and features
# excluded by the global read/occupancy criteria.
all_detected_classification <- tibble(
  taxon_id = colnames(otu_locality),
  total_reads = as.numeric(global_reads[colnames(otu_locality)]),
  n_localities = as.integer(n_localities[colnames(otu_locality)]),
  passes_global_eligibility = colnames(otu_locality) %in% eligible_ids,
  core_associated = colnames(otu_locality) %in% eligible_ids[is_core],
  classification = case_when(
    core_associated ~ "Eligible core-associated",
    passes_global_eligibility ~ "Eligible non-core",
    TRUE ~ "Excluded by global eligibility criteria"
  )
)
write_tsv(all_detected_classification, file.path(TABLE_DIR, "09A_all_detected_feature_status.tsv"))
write_tsv(tibble(input = c("taxonomic_counts", "metadata", "sample_manifest"),
                 path = c(PS_RDS, METADATA_CSV, SAMPLE_MANIFEST),
                 md5 = unname(tools::md5sum(c(PS_RDS, METADATA_CSV, SAMPLE_MANIFEST)))),
          file.path(TABLE_DIR, "09A_input_checksums.tsv"))

# Sensibilidad: core elegible frente a todo lo restante en la tabla detectada.
all_detected_summary <- tibble(
  taxon_id = colnames(otu_locality),
  total_reads = as.numeric(global_reads[colnames(otu_locality)]),
  component = if_else(
    taxon_id %in% taxon_classification$taxon_id[
      taxon_classification$core_associated
    ],
    "Core-associated eligible taxa",
    "All remaining detected taxa"
  )
) %>%
  group_by(component) %>%
  summarise(
    richness = n(),
    classified_reads = sum(total_reads),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(richness, classified_reads),
    names_to = "metric_code",
    values_to = "value"
  ) %>%
  group_by(metric_code) %>%
  mutate(
    denominator = sum(value),
    percentage = 100 * value / denominator
  ) %>%
  ungroup()

write_tsv(
  figure_summary %>% mutate(across(where(is.factor), as.character)),
  file.path(TABLE_DIR, "09A_Figure_3a_component_summary.tsv")
)
write_tsv(
  all_detected_summary,
  file.path(TABLE_DIR, "09A_all_detected_taxa_denominator_sensitivity.tsv")
)

# ============================== 7) Grafico ===================================

plot_df <- figure_summary %>%
  arrange(metric, component) %>%
  mutate(label = sprintf("%.1f%%", percentage))

p <- ggplot(
  plot_df,
  aes(x = percentage, y = metric, fill = component)
) +
  geom_col(
    width = 0.56,
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
    values = c(
      "Core-associated" = CORE_COLOUR,
      "Non-core" = NONCORE_COLOUR
    ),
    breaks = c("Core-associated", "Non-core"),
    name = NULL,
    drop = FALSE
  ) +
  scale_x_continuous(
    name = "Percentage (%)",
    limits = c(0, 100),
    breaks = seq(0, 100, 25),
    labels = label_percent(scale = 1, accuracy = 1),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_discrete(name = NULL) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y = element_text(colour = "black", size = 11),
    axis.text.x = element_text(colour = "black", size = 10),
    axis.title.x = element_text(size = 11, margin = margin(t = 7)),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.justification = "center",
    legend.text = element_text(size = 10),
    legend.key.width = grid::unit(0.55, "cm"),
    legend.spacing.x = grid::unit(0.18, "cm"),
    plot.margin = margin(8, 12, 4, 8)
  )

png_file  <- file.path(PLOT_DIR, "Figure_3a_taxonomic_core_richness_abundance.png")
pdf_file  <- file.path(PLOT_DIR, "Figure_3a_taxonomic_core_richness_abundance.pdf")
tiff_file <- file.path(PLOT_DIR, "Figure_3a_taxonomic_core_richness_abundance.tiff")

ggsave(png_file,  p, width = 7.2, height = 3.4, units = "in", dpi = 600, bg = "white")
ggsave(pdf_file,  p, width = 7.2, height = 3.4, units = "in", device = cairo_pdf,
       bg = "white")
ggsave(tiff_file, p, width = 7.2, height = 3.4, units = "in", dpi = 600,
       compression = "lzw", bg = "white")

# =========================== 8) Provenance ===================================

parameters <- tibble(
  parameter = c(
    "minimum_global_reads",
    "minimum_global_localities",
    "minimum_realm_localities_to_support_core",
    "core_realm_prevalence",
    "core_mean_relative_abundance",
    "taxon_identifier_harmonization",
    "classification_unit",
    "primary_denominator",
    "sample_key",
    "locality_metadata_column",
    "expected_n_samples",
    "bar_label_colour",
    "bar_label_size"
  ),
  value = c(
    MIN_GLOBAL_READS,
    MIN_GLOBAL_LOCALITIES,
    MIN_REALM_LOCALITIES,
    CORE_REALM_PREVALENCE,
    CORE_MEAN_REL_ABUND,
    "Abbreviated terminal labels merged 1:1 into full lineage IDs; counts summed before filtering",
    "REALM || Locality; replicate counts summed",
    "Eligible taxa after >=50 reads and >=5 localities",
    "run_accession",
    locality_col,
    EXPECTED_N_SAMPLES,
    BAR_LABEL_COLOUR,
    BAR_LABEL_SIZE
  )
)

core_richness_pct <- figure_summary %>%
  filter(metric_code == "richness", component == "Core-associated") %>%
  pull(percentage)

core_reads_pct <- figure_summary %>%
  filter(metric_code == "classified_reads", component == "Core-associated") %>%
  pull(percentage)

run_summary <- tibble(
  item = c(
    "script", "phyloseq", "metadata", "n_input_samples", "n_locality_unresolved", "n_samples", "n_locality_labels",
    "n_realm_locality_units", "n_locality_labels_spanning_multiple_realms",
    "n_realms", "minimum_realm_localities_to_support_core",
    "n_qualifying_realms", "qualifying_realms",
    "n_taxon_ids_merged", "reads_recovered_from_merged_ids",
    "n_taxa_detected", "n_taxa_eligible", "n_core_associated",
    "n_non_core", "core_richness_percent", "core_classified_reads_percent",
    "output_directory"
  ),
  value = c(
    SCRIPT_TAG, PS_RDS, METADATA_CSV, N_INPUT_SAMPLES, N_LOCALITY_UNRESOLVED, nrow(meta_aligned),
    n_distinct(meta_aligned$Locality), nrow(otu_locality),
    n_locality_realm_conflicts,
    length(realm_levels), MIN_REALM_LOCALITIES, length(qualifying_realms),
    paste(qualifying_realms, collapse = " | "),
    n_ids_merged, reads_in_ids_merged,
    ncol(otu_locality), ncol(otu_eligible), sum(is_core),
    sum(!is_core), sprintf("%.6f", core_richness_pct),
    sprintf("%.6f", core_reads_pct), OUT_DIR
  )
)

write_tsv(realm_counts, file.path(TABLE_DIR, "09A_realm_locality_counts.tsv"))
write_tsv(parameters, file.path(TABLE_DIR, "09A_core_definition_parameters.tsv"))
write_tsv(run_summary, file.path(TABLE_DIR, "09A_run_summary.tsv"))

saveRDS(
  list(
    parameters = parameters,
    metadata_aligned = meta_aligned,
    sample_inclusion = locality_inclusion,
    taxon_id_harmonization_audit = taxon_id_harmonization_audit,
    locality_realm_assignment_audit = locality_realm_audit,
    locality_metadata = locality_meta,
    locality_counts = otu_locality,
    realm_counts = realm_counts,
    qualifying_realms = qualifying_realms,
    taxon_classification = taxon_classification,
    realm_prevalence = realm_prevalence,
    realm_mean_relative_abundance = realm_mean_rel_abund,
    minimum_realm_localities_sensitivity =
      minimum_realm_localities_sensitivity,
    figure_summary = figure_summary,
    plot = p
  ),
  file.path(RDS_DIR, "09A_Figure_3a_analysis_objects.rds")
)

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
writeLines(OUT_DIR, LATEST_PTR)

msg(sprintf("Core-associated richness: %.1f%%", core_richness_pct))
msg(sprintf("Core-associated classified reads: %.1f%%", core_reads_pct))
msg("Figure: ", png_file)
msg("=== SCRIPT COMPLETED SUCCESSFULLY ===")
