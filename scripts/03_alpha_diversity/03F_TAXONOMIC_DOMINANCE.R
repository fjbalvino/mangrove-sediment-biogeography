#!/usr/bin/env Rscript
# Descriptive concentration of abundance in the three leading annotation groups.
# Uses the supplied standardized count matrix; no additional rarefaction is applied.
options(stringsAsFactors = FALSE, scipen = 8)
suppressPackageStartupMessages(library(phyloseq))
root <- Sys.getenv("BIOGEOGRAPHY_BASE", unset = ".")
results <- Sys.getenv("BIOGEOGRAPHY_RESULTS", unset = file.path(root, "results"))
alpha_dir <- Sys.getenv("ALPHA_INPUT_DIR", unset = file.path(results, "alpha"))
out <- Sys.getenv("DOMINANCE_OUTPUT_DIR", unset = file.path(alpha_dir, "dominance"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
ps_path <- Sys.getenv("MAGENTA_PHYLOSEQ", unset = Sys.getenv("MAGENTA_PHYLOSEQ_405", unset = file.path(root, "data", "phyloseq_taxonomic.rds")))
ps <- readRDS(ps_path)
meta <- read.csv(file.path(alpha_dir, "04a_alpha_diversity_filtered_samplewise.csv"), check.names = FALSE)
stopifnot("sampling_site_id" %in% names(meta))
meta <- meta[is.finite(meta$lat) & is.finite(meta$lon) & abs(meta$lat) <= 90 & abs(meta$lon) <= 180 & !is.na(meta$REALM) & !is.na(meta$sampling_site_id) & nzchar(meta$sampling_site_id), ]
meta$site_id <- meta$sampling_site_id
pm <- as(sample_data(ps), "data.frame")
stopifnot("run_accession" %in% names(pm), !anyNA(pm$run_accession), !anyDuplicated(pm$run_accession), !anyDuplicated(meta$run_accession))
ix <- match(meta$run_accession, pm$run_accession)
stopifnot(!anyNA(ix))
meta$.ps_sample_id <- rownames(pm)[ix]
write.csv(meta, file.path(out, "dominance_sample_manifest.csv"), row.names = FALSE)
counts<-as(otu_table(ps),'matrix');if(taxa_are_rows(ps))counts<-t(counts)
counts<-counts[meta$.ps_sample_id,,drop=FALSE];counts<-counts[,colSums(counts)>0,drop=FALSE]
tt<-as(tax_table(ps),'matrix')[colnames(counts),,drop=FALSE]
write.csv(data.frame(run_accession=meta$run_accession, supplied_reads=rowSums(counts)), file.path(out, 'dominance_input_library_sizes.csv'), row.names=FALSE)
stopifnot(nrow(counts)>0L,all(rowSums(counts)>0),all(is.finite(counts)),all(counts>=0))
tax_clean<-tt;tax_clean[is.na(tax_clean)]<-'<missing>'
key<-apply(tax_clean,1,paste,collapse='|')
annotation_counts<-t(rowsum(t(counts),group=key,reorder=FALSE))
first<-match(colnames(annotation_counts),key)
taxonomy<-data.frame(annotation_key=colnames(annotation_counts),representative_feature_id=colnames(counts)[first],tt[first,,drop=FALSE],n_source_features=as.numeric(table(factor(key,levels=colnames(annotation_counts)))),check.names=FALSE,row.names=NULL)
stopifnot(max(abs(rowSums(annotation_counts)-rowSums(counts)))==0)
codes<-c('Central Indo-Pacific'='CIP','Indo-Pacific'='IP','Temperate Northern Atlantic'='TNA','Tropical Atlantic'='TA','Tropical Eastern Pacific'='TEP','Western Indo-Pacific'='WIP')
ord <- read.csv(file.path(alpha_dir,'04a_realm_longitude_order_west_to_east.csv'))
realm_order <- ord$REALM[order(ord$west_to_east_position)]
gr<-split(seq_len(nrow(meta)),meta$site_id)
stopifnot(all(vapply(gr,function(i)length(unique(meta$REALM[i]))==1L,logical(1))))
site_meta<-do.call(rbind,lapply(gr,function(i)data.frame(site_id=meta$site_id[i[1]],REALM=meta$REALM[i[1]],n_metagenomes=length(i))))
row.names(site_meta)<-NULL
all_sample<-list();all_site<-list();all_summary<-list();tops<-list()
for(method in c('retained_taxonomic_features','identical_taxonomic_annotation_groups')) {
  x<-if(method=='retained_taxonomic_features')counts else annotation_counts
  rel<-x/rowSums(x)
  c3<-apply(rel,1,function(v)sum(sort(v,decreasing=TRUE)[1:3]))
  ss<-data.frame(run_accession=meta$run_accession,site_id=meta$site_id,REALM=meta$REALM,method=method,C3=c3)
  sm<-site_meta;sm$method<-method;sm$median_sample_C3<-vapply(gr,function(i)median(c3[i]),numeric(1))
  all_sample[[method]]<-ss;all_site[[method]]<-sm
  all_summary[[method]]<-do.call(rbind,lapply(c('ALL_COORDINATE_SITES',realm_order),function(r){dd<-if(r=='ALL_COORDINATE_SITES')sm else sm[sm$REALM==r,];data.frame(method=method,REALM=r,realm_code=if(r=='ALL_COORDINATE_SITES')'ALL' else unname(codes[r]),n_sites=nrow(dd),n_metagenomes=sum(dd$n_metagenomes),realm_at_least_5_sites=nrow(dd)>=5,median_site_C3=median(dd$median_sample_C3),lower_quartile=unname(quantile(dd$median_sample_C3,.25)),upper_quartile=unname(quantile(dd$median_sample_C3,.75)),minimum=min(dd$median_sample_C3),maximum=max(dd$median_sample_C3))}))
  if(method=='identical_taxonomic_annotation_groups') {
    site_rel<-do.call(rbind,lapply(gr,function(i)colMeans(rel[i,,drop=FALSE])))
    for(r in c('ALL_COORDINATE_SITES',realm_order)) {
      keep<-if(r=='ALL_COORDINATE_SITES')rep(TRUE,nrow(site_meta)) else site_meta$REALM==r
      avg<-colMeans(site_rel[keep,,drop=FALSE]);oi<-head(order(avg,decreasing=TRUE),3)
      top<-cbind(data.frame(REALM=r,realm_code=if(r=='ALL_COORDINATE_SITES')'ALL' else unname(codes[r]),n_sites=sum(keep),rank=1:3,site_balanced_mean_relative_abundance=avg[oi],top3_realm_mean_sum=sum(avg[oi]),row.names=NULL),taxonomy[oi,,drop=FALSE])
      tops[[r]]<-top
    }
  }
}
write_out<-function(x,n)write.csv(x,file.path(out,n),row.names=FALSE,na='')
write_out(do.call(rbind,all_sample),'dominance_C3_samplewise.csv')
write_out(do.call(rbind,all_site),'dominance_C3_site_medians.csv')
summary<-do.call(rbind,all_summary);top<-do.call(rbind,tops)
write_out(summary,'dominance_C3_descriptive_by_realm.csv')
write_out(top,'dominance_top3_annotation_groups_site_balanced.csv')
compact<-summary[summary$method=='identical_taxonomic_annotation_groups'&summary$REALM!='ALL_COORDINATE_SITES',c('REALM','realm_code','n_sites','n_metagenomes','realm_at_least_5_sites','median_site_C3','lower_quartile','upper_quartile')]
compact$median_C3_percent<-100*compact$median_site_C3
compact$lower_quartile_percent<-100*compact$lower_quartile
compact$upper_quartile_percent<-100*compact$upper_quartile
compact$C3_median_IQR_percent<-sprintf('%.2f (%.2f–%.2f)',compact$median_C3_percent,compact$lower_quartile_percent,compact$upper_quartile_percent)
write_out(compact,'dominance_C3_compact_realms.csv')
write_out(taxonomy,'dominance_taxonomy_annotation_group_map.csv')
mapping<-data.frame(feature_id=colnames(counts),annotation_key=key,tt,check.names=FALSE)
write_out(mapping,'dominance_retained_feature_taxonomy_map.csv')
audit<-data.frame(n_metagenomes=nrow(counts),n_sites=nrow(site_meta),n_retained_nonzero_features=ncol(counts),n_exact_annotation_groups=ncol(annotation_counts),n_annotation_groups_multiple_features=sum(taxonomy$n_source_features>1),max_features_per_annotation_group=max(taxonomy$n_source_features),max_count_conservation_error=max(abs(rowSums(annotation_counts)-rowSums(counts))),max_C3_change_on_aggregation=max(abs(all_sample[['retained_taxonomic_features']]$C3-all_sample[['identical_taxonomic_annotation_groups']]$C3)),n_sites_median_C3_above_half=sum(all_site[['identical_taxonomic_annotation_groups']]$median_sample_C3>.5),new_rarefaction=FALSE)
write_out(audit,'dominance_taxonomic_resolution_audit.csv')
print(audit);print(summary);print(top[,c('realm_code','rank','Species','Genus','n_source_features','site_balanced_mean_relative_abundance','top3_realm_mean_sum')],row.names=FALSE)

capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
writeLines(c(
  "C3 is the sum of the three largest within-profile relative abundances.",
  "Features with identical complete taxonomic annotations are combined; unresolved labels are retained.",
  "Each resolved sampling site contributes the median C3 of its metagenomes to realm summaries.",
  "Realm-leading identities are ranked separately by the mean relative abundance across equally weighted site means.",
  "The sum of realm-leading mean abundances and the median per-profile C3 are distinct summaries.",
  "These analyses use the supplied count matrix; library sizes are exported, and no new rarefaction is applied."
), file.path(out, "analysis_notes.txt"))
