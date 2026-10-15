#!/usr/bin/env Rscript
# Independent base-R verification of Python process scores and map estimates.
args <- commandArgs(trailingOnly=TRUE)
getarg <- function(key) {i <- match(key,args); if(is.na(i)) stop(key); args[i+1L]}
counts_file <- getarg('--counts'); meta_file <- getarg('--metadata')
tables <- getarg('--tables'); catalog_file <- getarg('--catalog')
readtab <- function(n) read.delim(file.path(tables,n),check.names=FALSE,stringsAsFactors=FALSE)
meta <- read.csv(meta_file,check.names=FALSE,stringsAsFactors=FALSE)
raw <- read.delim(gzfile(counts_file),check.names=FALSE)
rownames(raw) <- raw[[1]]
mat <- as.matrix(raw[,meta$run_accession,drop=FALSE])
master <- readtab('180_KO_master_classification.tsv')
shell <- master$KO[master$prevalence_class=='Shell']
clr <- t(log(mat[shell,,drop=FALSE]+1)); clr <- clr-rowMeans(clr)
catalog <- read.delim(catalog_file,check.names=FALSE,stringsAsFactors=FALSE)
catalog <- catalog[catalog$KO %in% shell,]
processes <- unique(catalog$process_id)
scores <- sapply(processes,function(p) rowMeans(clr[,catalog$KO[catalog$process_id==p],drop=FALSE]))
z <- scale(scores)
reported <- readtab('181_sample_process_scores_z.tsv')
reported <- as.matrix(reported[match(meta$run_accession,reported$run_accession),processes])
score_error <- max(abs(z-reported))
stopifnot(score_error<1e-12)
draws <- readtab('181_process_sampling_manifest.tsv')
out <- readtab('181_process_balanced_effects.tsv')
log_total <- log10(colSums(mat))
errors <- numeric(nrow(draws))
for(i in seq_len(nrow(draws))) {
  ids <- strsplit(draws$selected_run_accessions[i],';',fixed=TRUE)[[1]]
  ix <- match(ids,meta$run_accession)
  y <- z[ix,,drop=FALSE]; groups <- factor(meta$REALM[ix]); total <- log_total[ix]
  sst <- colSums(sweep(y,2,colMeans(y),'-')^2)
  fit <- lm.fit(model.matrix(~groups),y)
  raw_r2 <- 1-colSums(fit$residuals^2)/sst
  reduced <- lm.fit(model.matrix(~total),y)
  full <- lm.fit(model.matrix(~total+groups),y)
  conditional <- (colSums(reduced$residuals^2)-colSums(full$residuals^2))/sst
  expected <- out[out$design==draws$design[i] & out$iteration==draws$iteration[i],]
  expected <- expected[match(processes,expected$process_id),]
  errors[i] <- max(abs(raw_r2-expected$R2),abs(conditional-expected$adjusted_unique_R2))
}
stopifnot(max(errors)<1e-12)
units <- paste(meta$REALM,meta$Country,meta$core_locality,sep=' | ')
centroids <- rowsum(clr,units)/as.numeric(table(units)[rownames(rowsum(clr,units))])
center <- colMeans(centroids)
point <- sqrt(rowSums(sweep(centroids,2,center,'-')^2)/ncol(clr))
expected_map <- readtab('182_shell_locality_turnover_summary.tsv')
map_error <- max(abs(point[expected_map$locality_id]-expected_map$turnover_point_estimate))
stopifnot(map_error<1e-12)
check <- data.frame(check=c('standardized_process_scores','all_balanced_process_R2_and_conditional_R2','all_locality_RMS_Aitchison_distances'),
  n_values=c(length(z),nrow(out)*2,nrow(expected_map)),maximum_absolute_error=c(score_error,max(errors),map_error),pass=TRUE)
write.table(check,file.path(tables,'181_182_independent_R_verification.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
print(check);cat('PROCESS_MAP_INDEPENDENT_VERIFICATION_PASS\n')
