suppressPackageStartupMessages(library(data.table))
a<-commandArgs(trailingOnly=TRUE)
source_script<-a[1];run156<-a[2];root<-a[3];case<-a[4]
for(e in parse(source_script)) if(is.call(e) && identical(e[[1]],as.name('<-')) && identical(e[[2]],as.name('compute_proportionality_edges')))eval(e)
vertices<-fread(file.path(root,'inputs/candidates.tsv'))$KO
M<-fread(file.path(run156,'tables/156_KO_CLR_network_filtered.tsv.gz'))
Y<-t(as.matrix(M[match(vertices,M$KO),-1]));colnames(Y)<-vertices
D<-fread(file.path(root,'designs',paste0(case,'.tsv.gz')))
upper<-which(upper.tri(matrix(FALSE,length(vertices),length(vertices))),arr.ind=TRUE)
rows<-lapply(seq(1,500,by=10),function(b) {
  ids<-D[bootstrap_id==b]$run_accession
  out<-compute_proportionality_edges(Y[ids,,drop=FALSE],upper,2000,0)
  out[,bootstrap_id:=b];out
})
fwrite(rbindlist(rows),file.path(root,'influential_omission_R_check.tsv.gz'),sep='\t')
cat('Independent R edge selections completed for 50 draws\n')
