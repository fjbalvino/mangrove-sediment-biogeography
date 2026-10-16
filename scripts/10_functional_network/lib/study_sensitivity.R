suppressPackageStartupMessages({library(data.table);library(igraph)})
a <- commandArgs(trailingOnly=TRUE); mode <- a[1]; root <- normalizePath(a[2])
script <- sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1L])
source(file.path(dirname(normalizePath(script)), 'network_partition.R'))
if(mode=='design') {
  meta <- fread(file.path(root,'inputs/metadata.tsv'))
  clean <- function(x) {x<-trimws(as.character(x));x[is.na(x)|x=='']<-'Unassigned';gsub('[^A-Za-z0-9_.-]+','_',x)}
  meta[,BALANCE_STRATUM:=paste(clean(REALM),clean(Country),clean(Locality),sep='__')]
  stopifnot(!anyNA(meta$study_id),nrow(meta)==191)
  units <- meta[,.(n_removed=.N,realms_removed=paste(sort(unique(REALM)),collapse=';'),
    sample_ids=paste(run_accession,collapse=';')),by=study_id]
  setorder(units,-n_removed,study_id);units[,case:=sprintf('omit_%02d',seq_len(.N))]
  units[,unit_type:=ifelse(grepl('^PRJ',study_id),'BioProject','unresolved project label')]
  fwrite(units,file.path(root,'study_units.tsv'),sep='\t')
  dir.create(file.path(root,'designs'),showWarnings=FALSE)
  for(case in c('baseline',units$case)) {
    m<-if(case=='baseline')copy(meta) else meta[study_id!=units$study_id[match(case,units$case)]]
    set.seed(123);groups<-split(seq_len(nrow(m)),m$BALANCE_STRATUM)
    design<-rbindlist(lapply(seq_len(500),function(b) {
      ix<-unlist(lapply(groups,function(ix) if(length(ix)<=4)ix else sample(ix,4,replace=FALSE)),use.names=FALSE)
      m[unique(ix),.(bootstrap_id=b,run_accession,BALANCE_STRATUM)]
    }))
    stopifnot(min(design[,.N,by=bootstrap_id]$N)>=5)
    fwrite(design,file.path(root,'designs',paste0(case,'.tsv.gz')),sep='\t')
  }
}

if(mode=='evaluate') {
  run173 <- normalizePath(a[3]);run171<-normalizePath(a[4])
  ref <- fread(file.path(run173,'tables/173_network_nodes_annotated.tsv'))
  vertices <- fread(file.path(root,'inputs/candidates.tsv'))$KO
  edges_ref <- fread(file.path(run171,'tables/171_consensus_edges.tsv.gz'))
  key <- function(e) paste(pmin(e$KO1,e$KO2),pmax(e$KO1,e$KO2),sep='__')
  ref_keys<-key(edges_ref);units<-fread(file.path(root,'study_network_summary.tsv'))
  modules<-sort(unique(ref$module_label)); summaries<-list();scores<-list(); members<-list()
  one_case <- function(i) {
    ed<-fread(file.path(root,'networks',paste0(units$case[i],'_edges.tsv')))
    g<-graph_from_data_frame(ed[,.(KO1,KO2,weight=mean_abs_weight)],directed=FALSE,vertices=data.frame(name=vertices))
    fit<-fit_network_partition(g);m<-fit$membership;connected<-names(which(degree(g)>0))
    edgekeys<-key(ed)
    summaries[[i]]<-cbind(units[i],data.table(edge_Jaccard=length(intersect(ref_keys,edgekeys))/length(union(ref_keys,edgekeys)),
      primary_edge_recall=mean(ref_keys%in%edgekeys),major_node_retention=mean(ref$KO%in%connected),
      ARI_reference_major=compare(as.integer(factor(ref$module_label)),m[ref$KO],method='adjusted.rand'),
      n_major=sum(table(m[connected])>=10L),modularity=fit$Q,best_seed=fit$best_seed))
    scores[[i]]<-rbindlist(lapply(modules,function(mod) {
      kos<-ref[module_label==mod]$KO;mm<-m[kos];n<-length(mm)
      joint<-sum(table(mm)*(table(mm)-1))/(n*(n-1))
      groups<-unique(mm);jac<-vapply(groups,function(z) {k2<-names(m)[m==z];length(intersect(kos,k2))/length(union(kos,k2))},numeric(1))
      data.table(case=units$case[i],study_id=units$study_id[i],module_label=mod,n_KOs=n,
        within_module_coassignment=joint,best_overlap_Jaccard=max(jac),connected_KO_fraction=mean(kos%in%connected))
    }))
    members[[i]]<-data.table(case=units$case[i],KO=names(m),membership=as.integer(m),connected=names(m)%in%connected)
    cat('PARTITION',i,'/',nrow(units),units$study_id[i],'ARI',summaries[[i]]$ARI_reference_major,'\n');flush.console()
    list(summary=summaries[[i]],scores=scores[[i]],members=members[[i]])
  }
  ans<-parallel::mclapply(seq_len(nrow(units)),one_case,mc.cores=4,mc.set.seed=FALSE)
  if(any(vapply(ans,inherits,logical(1),'try-error')))stop('Study partition failed')
  summaries<-lapply(ans,`[[`,'summary');scores<-lapply(ans,`[[`,'scores');members<-lapply(ans,`[[`,'members')
  fwrite(rbindlist(summaries),file.path(root,'study_exclusion_results.tsv'),sep='\t')
  fwrite(rbindlist(scores),file.path(root,'study_exclusion_module_scores.tsv'),sep='\t')
  fwrite(rbindlist(members),file.path(root,'study_exclusion_memberships.tsv.gz'),sep='\t')
  writeLines(capture.output(sessionInfo()),file.path(root,'study_partition_sessionInfo.txt'))
}
