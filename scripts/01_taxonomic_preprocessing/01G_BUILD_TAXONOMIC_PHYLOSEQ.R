#!/usr/bin/env Rscript

# Build the canonical read-based MAGENTA phyloseq from direct Kraken2 counts.
# run_accession is the only sample key. Taxonomic features are represented by
# NCBI lineage strings and restricted to Bacteria and Archaea.

options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(data.table)
  library(phyloseq)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  exact <- which(args == flag)
  prefixed <- grep(paste0("^", flag, "="), args)
  if (length(exact) > 0L) {
    index <- exact[[length(exact)]]
    if (index == length(args)) stop(flag, " requiere un valor.", call. = FALSE)
    return(args[[index + 1L]])
  }
  if (length(prefixed) > 0L) {
    return(sub(paste0("^", flag, "="), "", args[[prefixed[[length(prefixed)]]]]))
  }
  default
}

as_integer_arg <- function(flag, default) {
  value <- suppressWarnings(as.integer(get_arg(flag, default)))
  if (is.na(value)) stop(flag, " debe ser entero.", call. = FALSE)
  value
}

COUNT_TABLE <- get_arg("--kraken-count-table")
METADATA <- get_arg("--metadata")
BASE <- Sys.getenv("BIOGEOGRAPHY_BASE", unset = ".")
DATA_DIR <- Sys.getenv("MAGENTA_DATA_DIR", unset = file.path(BASE, "data"))
TAXDUMP_DIR <- get_arg("--taxdump-dir", Sys.getenv("NCBI_TAXDUMP_DIR", unset = file.path(DATA_DIR, "references", "taxdump")))
OUTPUT_RDS <- get_arg(
  "--out-rds",
  Sys.getenv("PHYLOSEQ_405", unset = file.path(DATA_DIR, "phyloseq_source.rds"))
)
AUDIT_DIR <- get_arg(
  "--audit-dir",
  file.path(
    Sys.getenv("BIOGEOGRAPHY_RESULTS", unset = file.path(BASE, "results")),
    paste0("01G_BUILD_TAXONOMIC_PHYLOSEQ_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  )
)
RAREFACTION_DEPTH <- as_integer_arg("--rarefaction-depth", "18145")
SEED <- as_integer_arg("--seed", "42")
EXPECTED_SAMPLES <- as_integer_arg("--expected-samples", "405")
EXPECTED_TAXA <- as_integer_arg("--expected-taxa", "24045")
OVERWRITE <- identical(tolower(get_arg("--overwrite", "false")), "true")

if (is.null(COUNT_TABLE) || is.null(METADATA)) {
  stop("Se requieren --kraken-count-table y --metadata.", call. = FALSE)
}

dir.create(AUDIT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(OUTPUT_RDS), recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(AUDIT_DIR, "01G_build_phyloseq.log")
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message", append = TRUE)
on.exit({
  while (sink.number(type = "message") > 0L) sink(type = "message")
  while (sink.number(type = "output") > 0L) sink(type = "output")
  close(log_con)
}, add = TRUE)

msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = "")))
}

clean_text <- function(x) {
  value <- trimws(as.character(x))
  value[tolower(value) %in% c("", "na", "n/a", "nan", "none", "null")] <- NA_character_
  value
}

find_taxdump <- function(root, filename) {
  candidates <- c(
    file.path(root, filename),
    file.path(root, "taxonomy", filename),
    file.path(root, "taxdump", filename)
  )
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates) != 1L) {
    stop("No se encontró de forma inequívoca ", filename, " bajo ", root, call. = FALSE)
  }
  normalizePath(candidates[[1L]], mustWork = TRUE)
}

collapse_metadata <- function(metadata, id_col = "run_accession") {
  if (!id_col %in% names(metadata)) stop("Metadata sin run_accession.", call. = FALSE)
  metadata[[id_col]] <- clean_text(metadata[[id_col]])
  if (anyNA(metadata[[id_col]])) stop("Metadata con run_accession vacío.", call. = FALSE)
  groups <- split(seq_len(nrow(metadata)), metadata[[id_col]])
  rows <- lapply(groups, function(index) {
    group <- metadata[index, , drop = FALSE]
    output <- group[1L, , drop = FALSE]
    for (column in names(group)) {
      values <- unique(na.omit(clean_text(group[[column]])))
      if (length(values) > 1L) {
        stop(
          "Metadata contradictoria para ", group[[id_col]][[1L]],
          " en columna ", column, ": ", paste(head(values, 5L), collapse = " | "),
          call. = FALSE
        )
      }
      output[[column]] <- if (length(values) == 1L) values[[1L]] else NA_character_
    }
    output
  })
  output <- rbindlist(rows, fill = TRUE)
  setorderv(output, id_col)
  if (anyDuplicated(output[[id_col]])) stop("No se resolvieron duplicados metadata.", call. = FALSE)
  output
}

msg("Kraken count table: ", COUNT_TABLE)
msg("Metadata: ", METADATA)
msg("Taxdump:  ", TAXDUMP_DIR)
msg("Output:   ", OUTPUT_RDS)

required_files <- c(COUNT_TABLE, METADATA)
if (any(!file.exists(required_files))) {
  stop("Faltan inputs: ", paste(required_files[!file.exists(required_files)], collapse = ", "), call. = FALSE)
}
if (file.exists(OUTPUT_RDS) && !OVERWRITE) {
  stop("El output ya existe; use --overwrite true: ", OUTPUT_RDS, call. = FALSE)
}

metadata_raw <- fread(METADATA, na.strings = c("", "NA", "N/A"), check.names = FALSE)
metadata <- collapse_metadata(metadata_raw)
if (!"EcologicalClass" %in% names(metadata)) {
  stop("Metadata sin EcologicalClass.", call. = FALSE)
}
metadata[, EcologicalClass := clean_text(EcologicalClass)]
metadata_eligible <- metadata[
  EcologicalClass == "Mangrove sediment" | run_accession == "TX1_3"
]
setorder(metadata_eligible, run_accession)

msg("Reading direct Kraken2 count table...")
long_counts <- fread(COUNT_TABLE, na.strings = c("", "NA", "N/A"))
required_counts <- c("run_accession", "taxonomy_id", "count")
if (!all(required_counts %in% names(long_counts))) {
  stop("Tabla Kraken sin columnas: ", paste(setdiff(required_counts, names(long_counts)), collapse = ", "), call. = FALSE)
}
long_counts[, run_accession := clean_text(run_accession)]
long_counts[, taxonomy_id := clean_text(taxonomy_id)]
long_counts[, count := suppressWarnings(as.numeric(count))]
if (anyNA(long_counts$run_accession) || anyNA(long_counts$taxonomy_id) ||
    anyNA(long_counts$count) || any(long_counts$count < 0) ||
    any(abs(long_counts$count - round(long_counts$count)) > 1e-8)) {
  stop("Tabla Kraken con run_accession, taxonomy_id o conteos inválidos.", call. = FALSE)
}
long_counts[, count := as.numeric(round(count))]
long_counts <- long_counts[, .(count = sum(count)), by = .(run_accession, taxonomy_id)]
long_counts <- long_counts[count > 0]

missing_counts <- setdiff(metadata_eligible$run_accession, unique(long_counts$run_accession))
unexpected_counts <- setdiff(unique(long_counts$run_accession), metadata_eligible$run_accession)
if (length(missing_counts) > 0L) {
  fwrite(data.table(run_accession = missing_counts), file.path(AUDIT_DIR, "01G_metadata_samples_missing_kraken_counts.tsv"), sep = "\t")
  stop(length(missing_counts), " muestras elegibles carecen de conteos Kraken2.", call. = FALSE)
}
if (length(unexpected_counts) > 0L) {
  fwrite(data.table(run_accession = unexpected_counts), file.path(AUDIT_DIR, "01G_kraken_samples_not_eligible.tsv"), sep = "\t")
}
long_counts <- long_counts[run_accession %in% metadata_eligible$run_accession]

sample_totals_before <- long_counts[, .(classified_reads_before_rarefaction = sum(count)), by = run_accession]
low_depth <- sample_totals_before[classified_reads_before_rarefaction < RAREFACTION_DEPTH]
fwrite(sample_totals_before, file.path(AUDIT_DIR, "01G_sample_totals_before_rarefaction.tsv"), sep = "\t")
fwrite(low_depth, file.path(AUDIT_DIR, "01G_samples_below_rarefaction_depth.tsv"), sep = "\t")
if (nrow(low_depth) > 0L) {
  stop(nrow(low_depth), " muestras por debajo de ", RAREFACTION_DEPTH, " lecturas.", call. = FALSE)
}

taxids <- sort(unique(long_counts$taxonomy_id))
samples <- metadata_eligible$run_accession
tax_index <- match(long_counts$taxonomy_id, taxids)
sample_index <- match(long_counts$run_accession, samples)
count_matrix <- matrix(0, nrow = length(taxids), ncol = length(samples), dimnames = list(taxids, samples))
count_matrix[cbind(tax_index, sample_index)] <- long_counts$count

names_file <- find_taxdump(TAXDUMP_DIR, "names.dmp")
nodes_file <- find_taxdump(TAXDUMP_DIR, "nodes.dmp")
msg("Reading NCBI taxdump...")
nodes <- fread(nodes_file, sep = "|", header = FALSE, select = 1:3, strip.white = TRUE)
setnames(nodes, c("taxid", "parent", "rank"))
nodes[, `:=`(taxid = as.character(taxid), parent = as.character(parent), rank = trimws(rank))]
names_dt <- fread(names_file, sep = "|", header = FALSE, select = c(1, 2, 4), strip.white = TRUE)
setnames(names_dt, c("taxid", "name", "name_class"))
names_dt[, `:=`(taxid = as.character(taxid), name = trimws(name), name_class = trimws(name_class))]
names_dt <- names_dt[name_class == "scientific name"]
if (anyDuplicated(nodes$taxid) || anyDuplicated(names_dt$taxid)) stop("Taxdump con taxid duplicados.", call. = FALSE)

parent_map <- setNames(nodes$parent, nodes$taxid)
rank_map <- setNames(nodes$rank, nodes$taxid)
name_map <- setNames(names_dt$name, names_dt$taxid)
target_ranks <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")

lineage_for_taxid <- function(taxid) {
  output <- setNames(rep(NA_character_, length(target_ranks)), target_ranks)
  current <- taxid
  visited <- character(0)
  while (!is.na(current) && nzchar(current) && !current %in% visited) {
    visited <- c(visited, current)
    current_rank <- unname(rank_map[[current]])
    current_name <- unname(name_map[[current]])
    if (!is.null(current_rank) && current_rank %in% target_ranks && !is.null(current_name)) {
      output[[current_rank]] <- current_name
    }
    if (identical(current, "1")) break
    next_id <- unname(parent_map[[current]])
    if (is.null(next_id) || identical(next_id, current)) break
    current <- next_id
  }
  output
}

msg("Resolving lineage for ", length(taxids), " taxa...")
taxonomy_matrix <- t(vapply(taxids, lineage_for_taxid, FUN.VALUE = setNames(rep(NA_character_, length(target_ranks)), target_ranks)))
colnames(taxonomy_matrix) <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
taxonomy <- data.frame(taxid = taxids, taxonomy_matrix, check.names = FALSE, stringsAsFactors = FALSE)

keep_domain <- taxonomy$kingdom %in% c("Bacteria", "Archaea")
fwrite(as.data.table(taxonomy[!keep_domain, , drop = FALSE]), file.path(AUDIT_DIR, "01G_non_bacteria_archaea_taxa_removed.tsv"), sep = "\t")
taxonomy <- taxonomy[keep_domain, , drop = FALSE]
count_matrix <- count_matrix[taxonomy$taxid, , drop = FALSE]

prefixes <- c(kingdom = "k__", phylum = "p__", class = "c__", order = "o__", family = "f__", genus = "g__", species = "s__")
lineage_id <- apply(taxonomy[names(prefixes)], 1L, function(row) {
  values <- ifelse(is.na(row) | !nzchar(row), "unclassified", row)
  paste0(prefixes, values, collapse = "|")
})
duplicates <- duplicated(lineage_id) | duplicated(lineage_id, fromLast = TRUE)
lineage_id[duplicates] <- paste0(lineage_id[duplicates], " [taxid=", taxonomy$taxid[duplicates], "]")
if (anyDuplicated(lineage_id)) stop("No fue posible generar lineage_id únicos.", call. = FALSE)
rownames(count_matrix) <- lineage_id
rownames(taxonomy) <- lineage_id

metadata_df <- as.data.frame(metadata_eligible)
rownames(metadata_df) <- metadata_df$run_accession
metadata_df <- metadata_df[colnames(count_matrix), , drop = FALSE]

ps <- phyloseq(
  otu_table(count_matrix, taxa_are_rows = TRUE),
  tax_table(as.matrix(taxonomy[, c("taxid", names(prefixes)), drop = FALSE])),
  sample_data(metadata_df)
)

set.seed(SEED)
ps <- suppressWarnings(
  rarefy_even_depth(
    ps,
    sample.size = RAREFACTION_DEPTH,
    rngseed = SEED,
    replace = FALSE,
    trimOTUs = TRUE,
    verbose = FALSE
  )
)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

if (nsamples(ps) != EXPECTED_SAMPLES) {
  stop("Muestras observadas=", nsamples(ps), "; esperadas=", EXPECTED_SAMPLES, call. = FALSE)
}
if (EXPECTED_TAXA > 0L && ntaxa(ps) != EXPECTED_TAXA) {
  stop("Taxa observados=", ntaxa(ps), "; esperados=", EXPECTED_TAXA, call. = FALSE)
}
if (!identical(sample_names(ps), metadata_eligible$run_accession)) {
  stop("El orden final no coincide con run_accession de metadata.", call. = FALSE)
}
if (any(sample_sums(ps) != RAREFACTION_DEPTH)) {
  stop("La rarefacción final no produjo bibliotecas uniformes.", call. = FALSE)
}

saveRDS(ps, OUTPUT_RDS, compress = "xz")
ps_check <- readRDS(OUTPUT_RDS)
if (nsamples(ps_check) != nsamples(ps) || ntaxa(ps_check) != ntaxa(ps)) {
  stop("Validación post-save falló.", call. = FALSE)
}

fwrite(as.data.table(taxonomy, keep.rownames = "taxon_id"), file.path(AUDIT_DIR, "01G_taxonomy_used.tsv"), sep = "\t")
fwrite(as.data.table(metadata_df, keep.rownames = "phyloseq_sample_name"), file.path(AUDIT_DIR, "01G_metadata_used.tsv"), sep = "\t")

summary <- data.table(
  metric = c(
    "output_rds", "n_samples", "n_taxa", "rarefaction_depth", "total_counts",
    "n_bacteria_taxa", "n_archaea_taxa", "run_accession_unique", "output_md5"
  ),
  value = c(
    OUTPUT_RDS, nsamples(ps), ntaxa(ps), RAREFACTION_DEPTH, sum(sample_sums(ps)),
    sum(as(tax_table(ps), "matrix")[, "kingdom"] == "Bacteria"),
    sum(as(tax_table(ps), "matrix")[, "kingdom"] == "Archaea"),
    !anyDuplicated(sample_data(ps)$run_accession),
    unname(tools::md5sum(OUTPUT_RDS))
  )
)
fwrite(summary, file.path(AUDIT_DIR, "01G_run_summary.tsv"), sep = "\t")
capture.output(sessionInfo(), file = file.path(AUDIT_DIR, "01G_sessionInfo.txt"))

msg("PASS: ", nsamples(ps), " samples x ", ntaxa(ps), " taxa")
msg("Counts: ", format(sum(sample_sums(ps)), scientific = FALSE))
msg("Output MD5: ", unname(tools::md5sum(OUTPUT_RDS)))
