#!/usr/bin/env Rscript
# Select reference accessions and reconcile metadata without changing counts.
options(stringsAsFactors=FALSE)
suppressPackageStartupMessages(library(phyloseq))
args <- commandArgs(trailingOnly=TRUE)
arg <- function(flag,default='') {
  i <- match(flag,args)
  if(is.na(i)) return(default)
  if(i==length(args)) stop(flag,' requires a value.')
  args[i+1L]
}
base <- arg('--base',Sys.getenv('BIOGEOGRAPHY_BASE','.'))
input <- arg('--input-rds',Sys.getenv('MAGENTA_SOURCE_PHYLOSEQ',''))
output <- arg('--output-rds',file.path(base,'data/phyloseq_taxonomic.rds'))
manifest <- arg('--manifest',file.path(base,'config/samples/taxonomic_accessions.txt'))
metadata <- arg('--metadata',Sys.getenv('MAGENTA_METADATA_CSV',file.path(base,'data/metadata_taxonomic.csv')))
stopifnot(file.exists(input),file.exists(manifest),file.exists(metadata))
ids <- trimws(readLines(manifest,warn=FALSE)); ids <- ids[nzchar(ids)&!startsWith(ids,'#')]
stopifnot(length(ids)>0L,!anyDuplicated(ids))
ps <- readRDS(input)
stopifnot(inherits(ps,'phyloseq'))
old <- as(sample_data(ps),'data.frame')
stopifnot('run_accession'%in%names(old),!anyNA(old$run_accession),!anyDuplicated(old$run_accession),all(ids%in%old$run_accession))
selected <- prune_samples(old$run_accession%in%ids,ps)
md <- read.csv(metadata,check.names=FALSE,na.strings=c('','NA'))
selected_md <- as(sample_data(selected),'data.frame')
stopifnot(!anyDuplicated(md$run_accession),all(ids%in%md$run_accession))
md <- md[match(selected_md$run_accession,md$run_accession),,drop=FALSE]
rownames(md) <- rownames(selected_md)
counts <- otu_table(selected); taxonomy <- tax_table(selected)
sample_data(selected) <- sample_data(md)
stopifnot(identical(counts,otu_table(selected)),identical(taxonomy,tax_table(selected)),setequal(md$run_accession,ids))
dir.create(dirname(output),recursive=TRUE,showWarnings=FALSE)
tmp <- tempfile(tmpdir=dirname(output),fileext='.rds')
saveRDS(selected,tmp,compress='gzip')
stopifnot(identical(selected,readRDS(tmp)),file.rename(tmp,output))
cat('Reference profiles:',nsamples(selected),'; features:',ntaxa(selected),'\n')
