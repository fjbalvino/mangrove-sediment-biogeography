#!/usr/bin/env Rscript

# ============================================================
# 153_FIG4A_FUNCTIONAL_KO_CLASSES_RECONCILED.R
#
# Objetivo
#   Recalcular y graficar las clases funcionales de KOs usando
#   una definición reconciliada con el core taxonómico:
#
#   1) absolute_global_core:
#      detectado en el 100% de los metagenomas funcionales.
#   2) near_global:
#      detectado en >=95% pero <100% de los metagenomas.
#   3) realm_core_only:
#      no pertenece a 1 o 2 y está presente en todas las
#      localidades de al menos un marine realm con >=3
#      localidades, con abundancia relativa media >=0.0001.
#   4) functional_periphery:
#      multi_realm_non_core + realm_specific.
#
# Reconciliación metodológica clave
#   - Llave primaria: run_accession.
#   - Las réplicas se suman dentro de REALM x Locality.
#   - La abundancia relativa se calcula después de agregar.
#   - Cada localidad recibe el mismo peso dentro de su realm.
#   - El universo de muestras lo define exclusivamente la matriz
#     funcional usada por los scripts 131-152: metagenomas que
#     superaron el control de calidad del ensamblaje.
#   - No se aplica balanceo geográfico ni submuestreo de Guadalupe.
#     Todas las muestras de Guadalupe presentes en la matriz KO se
#     conservan.
#   - No se copian los prefiltros taxonómicos de 50 lecturas y
#     cinco localidades porque las abundancias KO no son lecturas
#     taxonómicas. Se conservan todos los KOs no-cero de la matriz.
#
# Salidas principales
#   figures/Figure_4a_functional_KO_classes_reconciled.png
#   figures/Figure_4a_functional_KO_classes_reconciled.pdf
#   figures/Figure_4a_functional_KO_classes_reconciled.tiff
#   tables/153_KO_classification_primary_minloc3.tsv
#   tables/153_Figure_4a_plot_summary.tsv
#   tables/153_class_counts_by_min_localities.tsv
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

# -----------------------------
# Argumentos
# -----------------------------
get_arg <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- which(args == flag)
  if (length(hit) == 0L || hit[length(hit)] == length(args)) {
    return(default)
  }
  args[hit[length(hit)] + 1L]
}

parse_int_vector <- function(x) {
  out <- suppressWarnings(as.integer(strsplit(x, ",", fixed = TRUE)[[1]]))
  out <- sort(unique(out[is.finite(out) & out >= 1L]))
  if (length(out) == 0L) stop("minloc_sensitivity no contiene enteros validos.")
  out
}

BASE <- "/home/fjbalvino/Biogeografia_2026"
RESULTS_DIR <- file.path(BASE, "resultados")

METADATA_FILE <- get_arg(
  "--metadata",
  "/data/Magenta/data/metadata/Metadata_2026_copia.csv"
)

MATRIX_DIR <- get_arg(
  "--matrix_dir",
  "/data/Magenta/results/035.metagenome_assembly_annotation/final_results"
)

KO_MATRIX <- get_arg("--ko_matrix", NA_character_)
RUN_COL_ARG <- get_arg("--run_col", "run_accession")
REALM_COL_ARG <- get_arg("--realm_col", NA_character_)
LOCALITY_COL_ARG <- get_arg("--locality_col", NA_character_)

NEAR_GLOBAL_THRESHOLD <- as.numeric(
  get_arg("--near_global_threshold", "0.95")
)

MIN_MEAN_REL_REALM <- as.numeric(
  get_arg("--min_mean_rel_realm", "0.0001")
)

PRESENCE_THRESHOLD <- as.numeric(
  get_arg("--presence_threshold", "0")
)

PRIMARY_MIN_REALM_LOCALITIES <- as.integer(
  get_arg("--min_realm_localities", "3")
)

MINLOC_SENSITIVITY <- parse_int_vector(
  get_arg("--minloc_sensitivity", "1,2,3,4")
)

if (!PRIMARY_MIN_REALM_LOCALITIES %in% MINLOC_SENSITIVITY) {
  MINLOC_SENSITIVITY <- sort(unique(c(
    MINLOC_SENSITIVITY,
    PRIMARY_MIN_REALM_LOCALITIES
  )))
}

EXPECTED_SAMPLES <- as.integer(get_arg("--expected_samples", "194"))
EXPECTED_KOS <- as.integer(get_arg("--expected_kos", "17508"))
STRICT_EXPECTED <- as.integer(get_arg("--strict_expected", "1")) == 1L
LEGACY_ANNOTATION <- get_arg("--legacy_annotation", NA_character_)
DPI <- as.integer(get_arg("--dpi", "600"))
PANEL_WIDTH_MM <- as.numeric(get_arg("--panel_width_mm", "87"))
PANEL_HEIGHT_MM <- as.numeric(get_arg("--panel_height_mm", "100"))
BASE_FONT_SIZE <- as.numeric(get_arg("--base_font_size", "8"))
TEXT_SIZE_MM <- BASE_FONT_SIZE / 2.845276
SEED <- as.integer(get_arg("--seed", "123"))

if (!is.finite(NEAR_GLOBAL_THRESHOLD) ||
    NEAR_GLOBAL_THRESHOLD <= 0 ||
    NEAR_GLOBAL_THRESHOLD > 1) {
  stop("near_global_threshold debe estar en (0, 1].")
}

if (!is.finite(MIN_MEAN_REL_REALM) || MIN_MEAN_REL_REALM < 0) {
  stop("min_mean_rel_realm debe ser >= 0.")
}

if (!is.finite(PRIMARY_MIN_REALM_LOCALITIES) ||
    PRIMARY_MIN_REALM_LOCALITIES < 1L) {
  stop("min_realm_localities debe ser >= 1.")
}

if (!is.finite(BASE_FONT_SIZE) ||
    BASE_FONT_SIZE < 8 ||
    BASE_FONT_SIZE > 11) {
  stop("base_font_size debe estar entre 8 y 11 pt.")
}

set.seed(SEED)

# -----------------------------
# Directorios y log
# -----------------------------
SCRIPT_ID <- "153_FIG4A_FUNCTIONAL_KO_CLASSES_RECONCILED"
STAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUTDIR <- file.path(RESULTS_DIR, paste0(SCRIPT_ID, "_", STAMP))
FIG_DIR <- file.path(OUTDIR, "figures")
TAB_DIR <- file.path(OUTDIR, "tables")
RDS_DIR <- file.path(OUTDIR, "rds")
LOG_DIR <- file.path(OUTDIR, "logs")

for (d in c(OUTDIR, FIG_DIR, TAB_DIR, RDS_DIR, LOG_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

LATEST_FILE <- file.path(RESULTS_DIR, paste0("LATEST_", SCRIPT_ID, ".txt"))
writeLines(OUTDIR, LATEST_FILE)

LOG_FILE <- file.path(LOG_DIR, "153_run.log")
log_con <- file(LOG_FILE, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")

on.exit({
  while (sink.number(type = "message") > 0L) sink(type = "message")
  while (sink.number(type = "output") > 0L) sink(type = "output")
  close(log_con)
}, add = TRUE)

log_msg <- function(...) {
  cat(
    sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0(..., collapse = ""),
    "\n",
    sep = ""
  )
}

write_tsv <- function(x, path) {
  data.table::fwrite(
    as.data.table(x),
    file = path,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )
}

write_status <- function(status, extra = character()) {
  writeLines(
    c(
      status,
      paste0("timestamp=", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      paste0("output_directory=", OUTDIR),
      extra
    ),
    file.path(OUTDIR, "STATUS.txt")
  )
}

# -----------------------------
# Utilidades de columnas e IDs
# -----------------------------
is_missing_like <- function(x) {
  y <- trimws(as.character(x))
  is.na(x) | y == "" | tolower(y) %in% c("na", "n/a", "nan", "null", "none")
}

clean_chr <- function(x) {
  y <- trimws(as.character(x))
  y[is_missing_like(y)] <- NA_character_
  y
}

norm_name <- function(x) {
  y <- tolower(trimws(as.character(x)))
  y <- gsub("[^a-z0-9]+", "_", y)
  gsub("^_+|_+$", "", y)
}

first_nonmissing <- function(x) {
  y <- clean_chr(x)
  y <- y[!is.na(y)]
  if (length(y) == 0L) NA_character_ else y[1L]
}

resolve_col <- function(df, candidates, label, required = TRUE) {
  nms <- names(df)
  idx <- match(norm_name(candidates), norm_name(nms))
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0L) {
    if (required) {
      stop(
        "No se encontro ", label, ". Candidatas: ",
        paste(candidates, collapse = ", ")
      )
    }
    return(NA_character_)
  }
  nms[idx[1L]]
}

match_ids_strict <- function(input_ids, reference_ids) {
  input_ids <- clean_chr(input_ids)
  reference_ids <- clean_chr(reference_ids)

  if (anyDuplicated(reference_ids[!is.na(reference_ids)])) {
    stop("La referencia contiene run_accession duplicados.")
  }

  ref_norm <- norm_name(reference_ids)
  if (anyDuplicated(ref_norm[!is.na(ref_norm)])) {
    stop("run_accession colisionan despues de normalizar IDs.")
  }

  out <- data.table(
    input_id = input_ids,
    matched_id = NA_character_,
    match_mode = NA_character_
  )

  idx_exact <- match(input_ids, reference_ids)
  hit <- !is.na(idx_exact)
  out$matched_id[hit] <- reference_ids[idx_exact[hit]]
  out$match_mode[hit] <- "exact"

  remaining <- which(is.na(out$matched_id) & !is.na(input_ids))
  if (length(remaining) > 0L) {
    idx_ci <- match(tolower(input_ids[remaining]), tolower(reference_ids))
    hit_ci <- !is.na(idx_ci)
    out$matched_id[remaining[hit_ci]] <- reference_ids[idx_ci[hit_ci]]
    out$match_mode[remaining[hit_ci]] <- "case_insensitive"
  }

  remaining <- which(is.na(out$matched_id) & !is.na(input_ids))
  if (length(remaining) > 0L) {
    idx_norm <- match(norm_name(input_ids[remaining]), ref_norm)
    hit_norm <- !is.na(idx_norm)
    out$matched_id[remaining[hit_norm]] <- reference_ids[idx_norm[hit_norm]]
    out$match_mode[remaining[hit_norm]] <- "normalized"
  }

  out
}

# -----------------------------
# Lectura y detección de matriz KO
# -----------------------------
read_table_auto <- function(path) {
  if (!file.exists(path)) stop("No existe: ", path)
  as.data.frame(
    data.table::fread(path, data.table = FALSE, check.names = FALSE),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

guess_ko_matrix_file <- function(matrix_dir) {
  if (!dir.exists(matrix_dir)) stop("No existe matrix_dir: ", matrix_dir)

  files <- list.files(
    matrix_dir,
    pattern = "\\.(tsv|txt|csv)(\\.gz)?$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )

  if (length(files) == 0L) {
    stop("No se encontraron matrices tabulares en: ", matrix_dir)
  }

  bn <- basename(files)
  score <- rep(0, length(files))
  score <- score + ifelse(grepl("kegg", bn, ignore.case = TRUE), 4, 0)
  score <- score + ifelse(grepl("(^|[^a-z])ko([^a-z]|$)|ortholog", bn, ignore.case = TRUE), 4, 0)
  score <- score + ifelse(grepl("matrix|abundance|count", bn, ignore.case = TRUE), 2, 0)
  score <- score - ifelse(grepl("module|pathway|brite", bn, ignore.case = TRUE), 12, 0)
  score <- score - ifelse(grepl("summary|metadata|readme|log", bn, ignore.case = TRUE), 8, 0)

  audit <- data.table(file = files, basename = bn, score = score)
  setorder(audit, -score, basename)
  write_tsv(audit, file.path(TAB_DIR, "153_KO_matrix_candidates.tsv"))

  if (audit$score[1L] <= 0) {
    stop(
      "No fue posible identificar la matriz KO con confianza. ",
      "Pasa --ko_matrix explicitamente."
    )
  }

  audit$file[1L]
}

numeric_matrix_from_df <- function(df) {
  m <- as.matrix(data.frame(
    lapply(df, function(x) suppressWarnings(as.numeric(as.character(x)))),
    check.names = FALSE
  ))
  m[is.na(m)] <- 0
  storage.mode(m) <- "numeric"
  m
}

collapse_duplicate_kos <- function(mat, ko_ids) {
  if (!anyDuplicated(ko_ids)) {
    rownames(mat) <- ko_ids
    return(mat)
  }
  rowsum(mat, group = ko_ids, reorder = FALSE)
}

load_ko_matrix <- function(path, metadata_ids) {
  raw <- read_table_auto(path)
  ko_by_col <- vapply(
    raw,
    function(x) sum(grepl("^K[0-9]{5}$", clean_chr(x))),
    integer(1)
  )

  if (max(ko_by_col) >= 10L) {
    ko_col <- names(ko_by_col)[which.max(ko_by_col)]
    ko_ids <- clean_chr(raw[[ko_col]])
    keep_ko <- grepl("^K[0-9]{5}$", ko_ids)
    raw <- raw[keep_ko, , drop = FALSE]
    ko_ids <- ko_ids[keep_ko]

    candidate_cols <- setdiff(names(raw), ko_col)
    matching <- match_ids_strict(candidate_cols, metadata_ids)
    matching <- matching[!is.na(matched_id)]

    if (nrow(matching) == 0L) {
      stop("No se encontraron columnas funcionales asociables por run_accession.")
    }
    if (anyDuplicated(matching$matched_id)) {
      stop("Varias columnas de la matriz KO se asociaron al mismo run_accession.")
    }

    mat <- numeric_matrix_from_df(raw[, matching$input_id, drop = FALSE])
    colnames(mat) <- matching$matched_id
    mat <- collapse_duplicate_kos(mat, ko_ids)

    return(list(
      matrix = mat,
      matching = matching,
      orientation = "KOs_by_rows",
      ko_id_column = ko_col
    ))
  }

  ko_cols <- names(raw)[grepl("^K[0-9]{5}$", names(raw))]
  if (length(ko_cols) < 10L) {
    stop("No se pudo detectar la orientacion de la matriz KO.")
  }

  sample_candidate_cols <- setdiff(names(raw), ko_cols)
  if (length(sample_candidate_cols) == 0L) {
    stop("La matriz no contiene una columna candidata de run_accession.")
  }

  sample_match_counts <- vapply(raw[, sample_candidate_cols, drop = FALSE], function(x) {
    sum(!is.na(match_ids_strict(clean_chr(x), metadata_ids)$matched_id))
  }, integer(1))

  sample_col <- names(sample_match_counts)[which.max(sample_match_counts)]
  if (max(sample_match_counts) < 5L) {
    stop("No se encontro una columna de run_accession en la matriz KO.")
  }

  matching <- match_ids_strict(clean_chr(raw[[sample_col]]), metadata_ids)
  keep_sample <- !is.na(matching$matched_id)
  matching <- matching[keep_sample]

  if (anyDuplicated(matching$matched_id)) {
    stop("La matriz KO contiene run_accession duplicados.")
  }

  m_samples <- numeric_matrix_from_df(raw[keep_sample, ko_cols, drop = FALSE])
  rownames(m_samples) <- matching$matched_id
  mat <- t(m_samples)
  rownames(mat) <- ko_cols

  list(
    matrix = mat,
    matching = matching,
    orientation = "samples_by_rows",
    ko_id_column = NA_character_,
    sample_id_column = sample_col
  )
}

# -----------------------------
# Clasificación funcional
# -----------------------------
class_order_five <- c(
  "absolute_global_core",
  "near_global",
  "realm_core_only",
  "multi_realm_non_core",
  "realm_specific"
)

plot_order_four <- c(
  "absolute_global_core",
  "near_global",
  "realm_core_only",
  "functional_periphery"
)

plot_labels <- c(
  "absolute_global_core" = "Absolute global core",
  "near_global" = "Near-global",
  "realm_core_only" = "Realm-core only",
  "functional_periphery" = "Functional periphery"
)

plot_axis_labels <- c(
  "absolute_global_core" = "Absolute core",
  "near_global" = "Near-global",
  "realm_core_only" = "Realm-core",
  "functional_periphery" = "Periphery"
)

plot_palette <- c(
  "absolute_global_core" = "#7F0000",
  "near_global" = "#C62828",
  "realm_core_only" = "#EF5350",
  "functional_periphery" = "#7F7F7F"
)

classify_kos <- function(
    min_realm_localities,
    mat,
    locality_counts,
    locality_rel,
    locality_meta,
    global_prevalence_n,
    global_prevalence_prop,
    total_abundance,
    detected_realm_matrix) {

  kos <- rownames(mat)
  realm_counts <- locality_meta[, .(n_localities = .N), by = realm]
  qualifying_realms <- realm_counts[
    n_localities >= min_realm_localities,
    realm
  ]

  support_matrix <- matrix(
    FALSE,
    nrow = length(kos),
    ncol = length(qualifying_realms),
    dimnames = list(kos, qualifying_realms)
  )

  details <- vector("list", length(unique(locality_meta$realm)))
  all_realms <- sort(unique(locality_meta$realm))

  for (i in seq_along(all_realms)) {
    realm_i <- all_realms[i]
    idx <- which(locality_meta$realm == realm_i)
    nloc <- length(idx)
    n_present_localities <- colSums(
      locality_counts[idx, , drop = FALSE] > PRESENCE_THRESHOLD
    )
    present_all <- n_present_localities == nloc
    mean_rel <- colMeans(locality_rel[idx, , drop = FALSE])
    qualifies <- nloc >= min_realm_localities
    supports <- qualifies & present_all & mean_rel >= MIN_MEAN_REL_REALM

    if (qualifies) support_matrix[, realm_i] <- supports

    details[[i]] <- data.table(
      KO = kos,
      realm = realm_i,
      n_localities_realm = nloc,
      min_realm_localities = min_realm_localities,
      realm_qualifies = qualifies,
      n_present_localities = n_present_localities,
      present_all_localities = present_all,
      mean_relative_abundance_realm = mean_rel,
      supports_realm_core = supports
    )
  }

  realm_core_any <- if (ncol(support_matrix) == 0L) {
    rep(FALSE, length(kos))
  } else {
    rowSums(support_matrix) > 0L
  }

  realm_core_realms <- if (ncol(support_matrix) == 0L) {
    rep("", length(kos))
  } else {
    apply(support_matrix, 1L, function(z) {
      paste(colnames(support_matrix)[z], collapse = ";")
    })
  }

  n_realms_detected <- rowSums(detected_realm_matrix)
  realms_detected <- apply(detected_realm_matrix, 1L, function(z) {
    paste(colnames(detected_realm_matrix)[z], collapse = ";")
  })

  is_absolute <- global_prevalence_n == ncol(mat)
  is_near <- !is_absolute & global_prevalence_prop >= NEAR_GLOBAL_THRESHOLD
  is_realm_core <- !is_absolute & !is_near & realm_core_any
  is_multi <- !is_absolute & !is_near & !is_realm_core & n_realms_detected > 1L
  is_specific <- !is_absolute & !is_near & !is_realm_core & n_realms_detected == 1L

  cls <- rep(NA_character_, length(kos))
  cls[is_absolute] <- "absolute_global_core"
  cls[is_near] <- "near_global"
  cls[is_realm_core] <- "realm_core_only"
  cls[is_multi] <- "multi_realm_non_core"
  cls[is_specific] <- "realm_specific"

  if (anyNA(cls)) {
    stop(
      "Hay ", sum(is.na(cls)),
      " KOs sin clasificar con minloc=", min_realm_localities, "."
    )
  }

  annotation <- data.table(
    KO = kos,
    min_realm_localities = min_realm_localities,
    functional_core_class = cls,
    global_prevalence_n = global_prevalence_n,
    global_prevalence_prop = global_prevalence_prop,
    n_functional_samples = ncol(mat),
    total_abundance = total_abundance,
    n_realms_detected = n_realms_detected,
    realms_detected = realms_detected,
    realm_core_realms = realm_core_realms,
    is_absolute_global_core = is_absolute,
    is_near_global = is_near,
    is_realm_core_only = is_realm_core,
    is_multi_realm_non_core = is_multi,
    is_realm_specific = is_specific,
    near_global_threshold = NEAR_GLOBAL_THRESHOLD,
    min_mean_relative_abundance_realm_core = MIN_MEAN_REL_REALM,
    presence_threshold = PRESENCE_THRESHOLD,
    abundance_weighting = "equal_weight_per_realm_locality"
  )

  list(
    annotation = annotation,
    details = rbindlist(details),
    qualifying_realms = qualifying_realms
  )
}

# -----------------------------
# Tema gráfico
# -----------------------------
theme_panel <- function(base_size = BASE_FONT_SIZE) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      axis.line = element_line(linewidth = 0.45, color = "black"),
      axis.ticks = element_line(linewidth = 0.4, color = "black"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = base_size),
      plot.margin = margin(5, 12, 5, 8)
    )
}

# -----------------------------
# Análisis principal
# -----------------------------
main <- function() {
  log_msg("Starting ", SCRIPT_ID)
  log_msg("Metadata: ", METADATA_FILE)
  log_msg("Matrix directory: ", MATRIX_DIR)
  log_msg("Primary min realm localities: ", PRIMARY_MIN_REALM_LOCALITIES)
  log_msg("Mean relative abundance threshold: ", MIN_MEAN_REL_REALM)
  log_msg("Near-global threshold: ", NEAR_GLOBAL_THRESHOLD)
  log_msg("Output: ", OUTDIR)

  if (!file.exists(METADATA_FILE)) stop("No existe metadata: ", METADATA_FILE)

  # 1. Metadata y auditoría de run_accession
  metadata_raw <- as.data.table(read_table_auto(METADATA_FILE))
  run_col <- resolve_col(metadata_raw, RUN_COL_ARG, "run_accession")
  realm_col <- if (!is.na(REALM_COL_ARG)) {
    resolve_col(metadata_raw, REALM_COL_ARG, "marine realm")
  } else {
    resolve_col(
      metadata_raw,
      c("REALM", "MarineRealm", "marine_realm", "realm"),
      "marine realm"
    )
  }
  locality_col <- if (!is.na(LOCALITY_COL_ARG)) {
    resolve_col(metadata_raw, LOCALITY_COL_ARG, "locality")
  } else {
    resolve_col(
      metadata_raw,
      c("Locality", "locality", "site", "location"),
      "locality"
    )
  }

  country_col <- resolve_col(
    metadata_raw,
    c("Country", "country", "COUNTRY"),
    "country",
    required = FALSE
  )

  meta <- metadata_raw[, .(
    run_accession = clean_chr(get(run_col)),
    realm = clean_chr(get(realm_col)),
    locality = clean_chr(get(locality_col))
  )]

  if (!is.na(country_col)) {
    meta[, country := clean_chr(metadata_raw[[country_col]])]
  } else {
    meta[, country := NA_character_]
  }

  meta <- meta[!is.na(run_accession)]

  meta_audit <- meta[, .(
    n_rows = .N,
    n_realms = uniqueN(realm[!is.na(realm)]),
    n_localities = uniqueN(locality[!is.na(locality)]),
    n_countries = uniqueN(country[!is.na(country)]),
    realms = paste(sort(unique(realm[!is.na(realm)])), collapse = ";"),
    localities = paste(sort(unique(locality[!is.na(locality)])), collapse = ";"),
    countries = paste(sort(unique(country[!is.na(country)])), collapse = ";")
  ), by = run_accession]
  write_tsv(meta_audit, file.path(TAB_DIR, "153_metadata_run_accession_audit.tsv"))

  bad_meta <- meta_audit[
    n_realms > 1L | n_localities > 1L | n_countries > 1L
  ]
  if (nrow(bad_meta) > 0L) {
    stop(
      "Hay ", nrow(bad_meta),
      " run_accession con realm/locality/country incompatibles."
    )
  }

  meta_unique <- meta[, .(
    realm = first_nonmissing(realm),
    locality = first_nonmissing(locality),
    country = first_nonmissing(country)
  ), by = run_accession]

  if (anyDuplicated(meta_unique$run_accession)) {
    stop("La metadata deduplicada conserva run_accession duplicados.")
  }

  # 2. Matriz KO emparejada por run_accession
  if (is.na(KO_MATRIX) || !nzchar(KO_MATRIX)) {
    KO_MATRIX <<- guess_ko_matrix_file(MATRIX_DIR)
  }
  log_msg("KO matrix selected: ", KO_MATRIX)

  loaded <- load_ko_matrix(KO_MATRIX, meta_unique$run_accession)
  mat <- loaded$matrix
  write_tsv(loaded$matching, file.path(TAB_DIR, "153_functional_sample_matching.tsv"))

  nonzero_ko <- rowSums(mat > PRESENCE_THRESHOLD) > 0L
  mat <- mat[nonzero_ko, , drop = FALSE]

  if (anyDuplicated(colnames(mat))) stop("run_accession duplicados en matriz KO.")
  if (anyDuplicated(rownames(mat))) stop("KO duplicados despues de colapsar IDs.")

  idx_meta <- match(colnames(mat), meta_unique$run_accession)
  if (anyNA(idx_meta)) stop("Fallo de alineacion final matriz-metadata.")
  meta_aligned <- copy(meta_unique[idx_meta])

  if (!identical(meta_aligned$run_accession, colnames(mat))) {
    stop("El orden de run_accession no coincide entre matriz y metadata.")
  }

  meta_aligned[, guadeloupe_flag := grepl(
    "guadeloupe",
    tolower(paste(country, locality)),
    fixed = TRUE
  )]

  meta_aligned[, cohort_source :=
    "functional_KO_matrix_assembly_QC_pass"]
  meta_aligned[, geographic_balancing_applied := FALSE]
  meta_aligned[, guadeloupe_subsampling_applied := FALSE]

  write_tsv(
    meta_aligned,
    file.path(TAB_DIR, "153_functional_matrix_cohort_audit.tsv")
  )

  log_msg(
    "Guadeloupe samples retained from functional matrix: ",
    meta_aligned[guadeloupe_flag == TRUE, .N]
  )

  missing_geo <- meta_aligned[is.na(realm) | is.na(locality)]
  write_tsv(missing_geo, file.path(TAB_DIR, "153_samples_missing_realm_or_locality.tsv"))
  if (nrow(missing_geo) > 0L) {
    stop(
      "Hay ", nrow(missing_geo),
      " metagenomas funcionales sin realm o locality."
    )
  }

  n_samples <- ncol(mat)
  n_kos <- nrow(mat)
  log_msg("Functional samples matched: ", n_samples)
  log_msg("Non-zero KOs: ", n_kos)

  if (STRICT_EXPECTED && n_samples != EXPECTED_SAMPLES) {
    stop("Se esperaban ", EXPECTED_SAMPLES, " muestras; se encontraron ", n_samples, ".")
  }
  if (STRICT_EXPECTED && n_kos != EXPECTED_KOS) {
    stop("Se esperaban ", EXPECTED_KOS, " KOs; se encontraron ", n_kos, ".")
  }

  write_tsv(meta_aligned, file.path(TAB_DIR, "153_metadata_aligned_by_run_accession.tsv"))

  # 3. Prevalencia global por metagenoma
  present_sample <- mat > PRESENCE_THRESHOLD
  global_prevalence_n <- rowSums(present_sample)
  global_prevalence_prop <- global_prevalence_n / n_samples
  total_abundance <- rowSums(mat)

  # 4. Agregación por REALM x Locality
  meta_aligned[, realm_locality_id := paste(realm, locality, sep = " || ")]
  locality_counts <- rowsum(
    t(mat),
    group = meta_aligned$realm_locality_id,
    reorder = FALSE
  )

  locality_meta <- unique(meta_aligned[, .(
    realm_locality_id,
    realm,
    locality
  )])

  if (anyDuplicated(locality_meta$realm_locality_id)) {
    stop("realm_locality_id no es unica despues de agregar replicas.")
  }

  locality_meta <- locality_meta[
    match(rownames(locality_counts), realm_locality_id)
  ]
  if (anyNA(locality_meta$realm_locality_id)) {
    stop("No se pudo alinear locality_counts con locality_meta.")
  }

  locality_totals <- rowSums(locality_counts)
  if (any(locality_totals <= 0)) {
    stop("Hay unidades realm-locality con abundancia KO total cero.")
  }
  locality_rel <- sweep(locality_counts, 1L, locality_totals, "/")

  realm_locality_counts <- locality_meta[, .(
    n_localities = .N
  ), by = realm][order(-n_localities, realm)]
  realm_locality_counts[, qualifies_primary :=
    n_localities >= PRIMARY_MIN_REALM_LOCALITIES]
  write_tsv(
    realm_locality_counts,
    file.path(TAB_DIR, "153_realm_locality_counts.tsv")
  )

  log_msg("Realm-locality units: ", nrow(locality_meta))
  log_msg(
    "Primary qualifying realms: ",
    paste(realm_locality_counts[qualifies_primary == TRUE, realm], collapse = " | ")
  )

  all_realms <- sort(unique(locality_meta$realm))
  detected_realm_matrix <- sapply(all_realms, function(realm_i) {
    idx <- which(locality_meta$realm == realm_i)
    colSums(locality_counts[idx, , drop = FALSE] > PRESENCE_THRESHOLD) > 0L
  })
  detected_realm_matrix <- as.matrix(detected_realm_matrix)
  rownames(detected_realm_matrix) <- rownames(mat)
  colnames(detected_realm_matrix) <- all_realms

  # 5. Sensibilidad 1-4 localidades y clasificación primaria minloc3
  sensitivity_results <- vector("list", length(MINLOC_SENSITIVITY))
  names(sensitivity_results) <- as.character(MINLOC_SENSITIVITY)
  count_list <- list()
  primary_result <- NULL

  for (minloc_i in MINLOC_SENSITIVITY) {
    log_msg("Classifying KOs with minloc=", minloc_i)
    res_i <- classify_kos(
      min_realm_localities = minloc_i,
      mat = mat,
      locality_counts = locality_counts,
      locality_rel = locality_rel,
      locality_meta = locality_meta,
      global_prevalence_n = global_prevalence_n,
      global_prevalence_prop = global_prevalence_prop,
      total_abundance = total_abundance,
      detected_realm_matrix = detected_realm_matrix
    )
    sensitivity_results[[as.character(minloc_i)]] <- res_i$annotation

    count_list[[as.character(minloc_i)]] <- res_i$annotation[, .(
      n_KOs = .N,
      total_abundance = sum(total_abundance)
    ), by = functional_core_class][
      , min_realm_localities := minloc_i
    ]

    if (minloc_i == PRIMARY_MIN_REALM_LOCALITIES) {
      primary_result <- res_i
    }
  }

  if (is.null(primary_result)) stop("No se genero la clasificacion primaria.")

  class_counts_all <- rbindlist(count_list, fill = TRUE)
  class_counts_all[, functional_core_class := factor(
    functional_core_class,
    levels = class_order_five
  )]
  setorder(class_counts_all, min_realm_localities, functional_core_class)
  class_counts_all[, functional_core_class := as.character(functional_core_class)]
  class_counts_all[, abundance_fraction := total_abundance / sum(total_abundance),
                   by = min_realm_localities]

  write_tsv(
    class_counts_all,
    file.path(TAB_DIR, "153_class_counts_by_min_localities.tsv")
  )
  write_tsv(
    primary_result$annotation,
    file.path(TAB_DIR, "153_KO_classification_primary_minloc3.tsv")
  )
  write_tsv(
    primary_result$details,
    file.path(TAB_DIR, "153_KO_realm_support_audit_primary.tsv")
  )

  # 6. Transiciones minloc1 -> minloc3
  if (all(c("1", as.character(PRIMARY_MIN_REALM_LOCALITIES)) %in%
          names(sensitivity_results))) {
    ann_1 <- sensitivity_results[["1"]][, .(
      KO,
      class_minloc1 = functional_core_class
    )]
    ann_p <- primary_result$annotation[, .(
      KO,
      class_primary = functional_core_class
    )]
    transitions <- merge(ann_1, ann_p, by = "KO", all = TRUE)
    transition_summary <- transitions[, .N, by = .(
      class_minloc1,
      class_primary
    )][order(class_minloc1, class_primary)]
    write_tsv(
      transitions,
      file.path(TAB_DIR, "153_KO_transitions_minloc1_to_primary.tsv")
    )
    write_tsv(
      transition_summary,
      file.path(TAB_DIR, "153_transition_summary_minloc1_to_primary.tsv")
    )
  }

  # 7. Comparación opcional contra clasificación histórica
  legacy_path <- LEGACY_ANNOTATION
  if (is.na(legacy_path) || !nzchar(legacy_path)) {
    latest_132 <- file.path(
      RESULTS_DIR,
      "LATEST_132_KO_functional_core_class_sensitivity.txt"
    )
    if (file.exists(latest_132)) {
      run132 <- trimws(readLines(latest_132, n = 1L, warn = FALSE))
      candidate <- file.path(
        run132,
        "tables",
        "132_KO_core_class_recommended_minloc3.tsv"
      )
      if (file.exists(candidate)) legacy_path <- candidate
    }
  }

  if (!is.na(legacy_path) && nzchar(legacy_path) && file.exists(legacy_path)) {
    legacy <- as.data.table(read_table_auto(legacy_path))

    write_tsv(
      data.table(
        column = names(legacy),
        normalized_column = norm_name(names(legacy))
      ),
      file.path(TAB_DIR, "153_legacy_columns_detected.tsv")
    )

    ko_col_legacy <- resolve_col(
      legacy,
      c("KO", "ko", "KO_id", "feature_id"),
      "legacy KO",
      required = FALSE
    )

    if (is.na(ko_col_legacy)) {
      ko_scores <- vapply(
        legacy,
        function(x) sum(grepl("^K[0-9]{5}$", clean_chr(x))),
        integer(1)
      )
      if (length(ko_scores) > 0L && max(ko_scores) > 0L) {
        ko_col_legacy <- names(ko_scores)[which.max(ko_scores)]
      }
    }

    class_col_legacy <- resolve_col(
      legacy,
      c(
        "functional_core_class",
        "recommended_functional_core_class",
        "recommended_class",
        "functional_core_class_minloc3",
        "core_class_minloc3",
        "core_class",
        "class"
      ),
      "legacy class",
      required = FALSE
    )

    if (is.na(class_col_legacy)) {
      class_candidates <- names(legacy)[
        grepl("class", norm_name(names(legacy)), fixed = TRUE)
      ]

      if (length(class_candidates) > 0L) {
        expected_classes <- c(
          class_order_five,
          "functional_periphery"
        )
        class_scores <- vapply(
          class_candidates,
          function(nm) {
            sum(clean_chr(legacy[[nm]]) %in% expected_classes)
          },
          integer(1)
        )
        if (max(class_scores) > 0L) {
          class_col_legacy <- class_candidates[which.max(class_scores)]
        }
      }
    }

    if (!is.na(ko_col_legacy) && !is.na(class_col_legacy)) {
      legacy_cmp <- legacy[, .(
        KO = clean_chr(get(ko_col_legacy)),
        legacy_class = clean_chr(get(class_col_legacy))
      )]
      legacy_cmp <- unique(legacy_cmp[!is.na(KO)])

      current_cmp <- primary_result$annotation[, .(
        KO,
        reconciled_class = functional_core_class
      )]
      comparison <- merge(legacy_cmp, current_cmp, by = "KO", all = TRUE)
      comparison[, changed := legacy_class != reconciled_class]
      comparison_summary <- comparison[, .N, by = .(
        legacy_class,
        reconciled_class
      )][order(legacy_class, reconciled_class)]

      write_tsv(
        comparison,
        file.path(TAB_DIR, "153_comparison_with_legacy_minloc3.tsv")
      )
      write_tsv(
        comparison_summary,
        file.path(TAB_DIR, "153_comparison_with_legacy_summary.tsv")
      )
      log_msg("Legacy comparison performed: ", legacy_path)
      log_msg("Legacy KO column: ", ko_col_legacy)
      log_msg("Legacy class column: ", class_col_legacy)
      log_msg("KOs changing class: ", comparison[changed == TRUE, .N])
    } else {
      log_msg(
        "Legacy comparison skipped: KO/class columns could not be resolved. ",
        "This optional audit does not affect the reconciled analysis."
      )
    }
  } else {
    log_msg("Legacy annotation not supplied/found; comparison skipped.")
  }

  # 8. Datos para el panel principal de cuatro clases
  panel_data <- copy(primary_result$annotation)
  panel_data[, plot_class := fifelse(
    functional_core_class %in% c("multi_realm_non_core", "realm_specific"),
    "functional_periphery",
    functional_core_class
  )]

  panel_summary <- panel_data[, .(
    n_KOs = .N,
    total_abundance = sum(total_abundance)
  ), by = plot_class]
  panel_summary[, richness_fraction := n_KOs / sum(n_KOs)]
  panel_summary[, abundance_fraction := total_abundance / sum(total_abundance)]
  panel_summary[, plot_class := factor(plot_class, levels = plot_order_four)]
  setorder(panel_summary, plot_class)
  panel_summary[, plot_class := as.character(plot_class)]
  panel_summary[, label := unname(plot_labels[plot_class])]

  shared_summary <- data.table(
    component = c("shared_core_associated", "functional_periphery"),
    n_KOs = c(
      panel_summary[plot_class != "functional_periphery", sum(n_KOs)],
      panel_summary[plot_class == "functional_periphery", sum(n_KOs)]
    ),
    total_abundance = c(
      panel_summary[plot_class != "functional_periphery", sum(total_abundance)],
      panel_summary[plot_class == "functional_periphery", sum(total_abundance)]
    )
  )
  shared_summary[, richness_fraction := n_KOs / sum(n_KOs)]
  shared_summary[, abundance_fraction := total_abundance / sum(total_abundance)]

  write_tsv(
    panel_summary,
    file.path(TAB_DIR, "153_Figure_4a_plot_summary.tsv")
  )
  write_tsv(
    shared_summary,
    file.path(TAB_DIR, "153_shared_vs_periphery_summary.tsv")
  )

  # 9. Gráfico principal sin título ni letra de panel
  plot_counts <- copy(panel_summary)
  plot_counts[, y_class := factor(
    plot_class,
    levels = rev(plot_order_four)
  )]

  p_counts <- ggplot(
    plot_counts,
    aes(x = n_KOs, y = y_class, color = plot_class)
  ) +
    geom_segment(
      aes(x = 1, xend = n_KOs, yend = y_class),
      linewidth = 1.15,
      alpha = 0.68
    ) +
    geom_point(size = 3.2) +
    geom_text(
      aes(label = comma(n_KOs)),
      hjust = -0.16,
      color = "black",
      size = TEXT_SIZE_MM
    ) +
    scale_x_log10(
      limits = c(1, 100000),
      breaks = c(10, 100, 1000, 10000),
      labels = c("10", "100", "1k", "10k"),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_discrete(labels = function(x) {
      unname(plot_axis_labels[as.character(x)])
    }) +
    scale_color_manual(values = plot_palette, guide = "none") +
    coord_cartesian(clip = "off") +
    labs(x = "Number of KOs (log scale)", y = NULL) +
    theme_panel(base_size = BASE_FONT_SIZE) +
    theme(
      axis.text.x = element_text(size = BASE_FONT_SIZE),
      axis.text.y = element_text(
        size = BASE_FONT_SIZE,
        lineheight = 1
      ),
      axis.title.x = element_text(
        size = BASE_FONT_SIZE,
        margin = margin(t = 5)
      ),
      plot.margin = margin(4, 7, 1, 4)
    )

  plot_abundance <- copy(panel_summary)
  plot_abundance[, plot_class := factor(
    plot_class,
    levels = plot_order_four
  )]
  plot_abundance[, pct_label := sprintf(
    "%.1f%%",
    100 * abundance_fraction
  )]
  plot_abundance[, order_index := match(
    as.character(plot_class),
    plot_order_four
  )]
  setorder(plot_abundance, order_index)
  plot_abundance[, xmax := cumsum(abundance_fraction)]
  plot_abundance[, xmin := shift(xmax, fill = 0)]
  plot_abundance[, xmid := (xmin + xmax) / 2]

  p_abundance <- ggplot(plot_abundance) +
    geom_rect(
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = 0.22,
        ymax = 0.78,
        fill = plot_class
      ),
      color = "white",
      linewidth = 0.5
    ) +
    geom_text(
      data = plot_abundance[abundance_fraction >= 0.06],
      aes(
        x = xmid,
        y = 0.50,
        label = pct_label
      ),
      inherit.aes = FALSE,
      color = "black",
      fontface = "bold",
      size = TEXT_SIZE_MM
    ) +
    geom_text(
      data = plot_abundance[abundance_fraction < 0.06],
      aes(
        x = xmid,
        y = 0.50,
        label = pct_label
      ),
      inherit.aes = FALSE,
      color = "black",
      fontface = "bold",
      angle = 90,
      size = TEXT_SIZE_MM
    ) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.25),
      labels = percent_format(accuracy = 1),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = 0.5,
      labels = "KO abundance",
      expand = expansion(mult = c(0, 0))
    ) +
    scale_fill_manual(
      values = plot_palette,
      guide = "none"
    ) +
    labs(x = "Annotated KO abundance (%)", y = NULL) +
    theme_panel(base_size = BASE_FONT_SIZE) +
    theme(
      axis.text.x = element_text(size = BASE_FONT_SIZE),
      axis.text.y = element_text(
        size = BASE_FONT_SIZE,
        lineheight = 1
      ),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      axis.title.x = element_text(
        size = BASE_FONT_SIZE,
        margin = margin(t = 5)
      ),
      legend.position = "none",
      plot.margin = margin(1, 5, 3, 4)
    )

  combined <- p_counts / p_abundance +
    plot_layout(heights = c(2.55, 0.95))

  png_file <- file.path(FIG_DIR, "Figure_4a_functional_KO_classes_reconciled.png")
  pdf_file <- file.path(FIG_DIR, "Figure_4a_functional_KO_classes_reconciled.pdf")
  tiff_file <- file.path(FIG_DIR, "Figure_4a_functional_KO_classes_reconciled.tiff")

  ggsave(
    png_file,
    combined,
    width = PANEL_WIDTH_MM,
    height = PANEL_HEIGHT_MM,
    units = "mm",
    dpi = DPI,
    bg = "white"
  )
  ggsave(
    pdf_file,
    combined,
    width = PANEL_WIDTH_MM,
    height = PANEL_HEIGHT_MM,
    units = "mm",
    device = grDevices::cairo_pdf,
    bg = "white"
  )
  ggsave(
    tiff_file,
    combined,
    width = PANEL_WIDTH_MM,
    height = PANEL_HEIGHT_MM,
    units = "mm",
    dpi = DPI,
    compression = "lzw",
    bg = "white"
  )

  # 10. Sensibilidad minloc para suplemento/QC
  sens_plot_data <- class_counts_all[
    functional_core_class %in% c(
      "absolute_global_core",
      "near_global",
      "realm_core_only"
    )
  ]
  sens_plot_data[, functional_core_class := factor(
    functional_core_class,
    levels = c("absolute_global_core", "near_global", "realm_core_only")
  )]

  p_sens <- ggplot(
    sens_plot_data,
    aes(
      x = min_realm_localities,
      y = n_KOs,
      color = functional_core_class,
      group = functional_core_class
    )
  ) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.6) +
    scale_y_log10(labels = comma) +
    scale_x_continuous(breaks = MINLOC_SENSITIVITY) +
    scale_color_manual(
      values = plot_palette[c(
        "absolute_global_core",
        "near_global",
        "realm_core_only"
      )],
      labels = unname(plot_labels[c(
        "absolute_global_core",
        "near_global",
        "realm_core_only"
      )])
    ) +
    labs(
      x = "Minimum localities per realm",
      y = "Number of KOs (log scale)"
    ) +
    theme_panel(base_size = 12)

  ggsave(
    file.path(FIG_DIR, "Supplementary_minloc_sensitivity_reconciled.png"),
    p_sens,
    width = 150,
    height = 105,
    units = "mm",
    dpi = DPI,
    bg = "white"
  )

  # 11. Parámetros, objetos y provenance
  parameters <- data.table(
    parameter = c(
      "metadata_file",
      "ko_matrix",
      "sample_key",
      "classification_unit",
      "replicate_handling",
      "realm_abundance_weighting",
      "functional_cohort_source",
      "assembly_quality_filter",
      "geographic_balancing_applied",
      "guadeloupe_subsampling_applied",
      "panel_width_mm",
      "panel_height_mm",
      "base_font_size",
      "near_global_threshold",
      "primary_minimum_realm_localities",
      "realm_core_prevalence",
      "realm_core_mean_relative_abundance",
      "presence_threshold",
      "n_functional_samples",
      "n_realm_locality_units",
      "n_KOs",
      "taxonomic_global_read_filter_applied",
      "taxonomic_global_locality_filter_applied",
      "seed"
    ),
    value = c(
      METADATA_FILE,
      KO_MATRIX,
      "run_accession",
      "REALM x Locality",
      "KO counts summed within locality",
      "equal weight per locality",
      "KO matrix used by scripts 131-152",
      "inherited from functional matrix construction",
      "NO",
      "NO",
      PANEL_WIDTH_MM,
      PANEL_HEIGHT_MM,
      BASE_FONT_SIZE,
      NEAR_GLOBAL_THRESHOLD,
      PRIMARY_MIN_REALM_LOCALITIES,
      1,
      MIN_MEAN_REL_REALM,
      PRESENCE_THRESHOLD,
      n_samples,
      nrow(locality_meta),
      n_kos,
      "NO",
      "NO",
      SEED
    )
  )
  write_tsv(parameters, file.path(TAB_DIR, "153_parameters.tsv"))

  saveRDS(
    list(
      primary_annotation = primary_result$annotation,
      panel_summary = panel_summary,
      shared_summary = shared_summary,
      class_counts_by_minloc = class_counts_all,
      metadata_aligned = meta_aligned,
      realm_locality_metadata = locality_meta,
      realm_locality_counts = locality_counts,
      realm_locality_relative_abundance = locality_rel,
      parameters = parameters
    ),
    file.path(RDS_DIR, "153_Figure_4a_analysis_objects.rds")
  )

  capture.output(
    sessionInfo(),
    file = file.path(LOG_DIR, "153_sessionInfo.txt")
  )

  primary_counts <- class_counts_all[
    min_realm_localities == PRIMARY_MIN_REALM_LOCALITIES
  ]
  log_msg("Primary class counts:")
  print(primary_counts)
  log_msg("Four-class Figure 4a summary:")
  print(panel_summary)
  log_msg("Figure: ", png_file)
  log_msg("=== SCRIPT COMPLETED SUCCESSFULLY ===")

  write_status(
    "PASS",
    c(
      paste0("n_functional_samples=", n_samples),
      paste0("n_realm_locality_units=", nrow(locality_meta)),
      paste0("n_KOs=", n_kos),
      paste0("figure=", png_file)
    )
  )
}

tryCatch(
  main(),
  error = function(e) {
    write_status("FAIL", paste0("error=", conditionMessage(e)))
    log_msg("ERROR: ", conditionMessage(e))
    stop(e)
  }
)
