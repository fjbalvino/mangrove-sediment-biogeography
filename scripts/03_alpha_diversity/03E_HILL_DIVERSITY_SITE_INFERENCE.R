#!/usr/bin/env Rscript
# Hill q0/q2 figures and primary realm comparisons using medians per resolved sampling site.
# Realms with fewer than five resolved sampling sites are retained in plots only.
options(stringsAsFactors = FALSE, scipen = 8, bitmapType = "cairo")
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
root <- Sys.getenv("BIOGEOGRAPHY_BASE", unset = ".")
results <- Sys.getenv("BIOGEOGRAPHY_RESULTS", unset = file.path(root, "results"))
input_dir <- Sys.getenv("ALPHA_INPUT_DIR", unset = file.path(results, "alpha"))
out <- Sys.getenv("ALPHA_SUPPLEMENT_OUTPUT_DIR", unset = file.path(input_dir, "supplementary"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out_dir <- out
dat <- read.csv(file.path(input_dir, "04a_alpha_diversity_filtered_samplewise.csv"), check.names = FALSE)
stopifnot(nrow(dat) > 0L, "sampling_site_id" %in% names(dat), !anyDuplicated(dat$run_accession), all(is.finite(as.matrix(dat[c("q0", "q1", "q2")]))))
codes<-c('Central Indo-Pacific'='CIP','Indo-Pacific'='IP','Temperate Northern Atlantic'='TNA','Tropical Atlantic'='TA','Tropical Eastern Pacific'='TEP','Western Indo-Pacific'='WIP')
ord <- read.csv(file.path(input_dir, '04a_realm_longitude_order_west_to_east.csv'))
realm_order <- ord$REALM[order(ord$west_to_east_position)]
coord_ok<-is.finite(dat$lat)&is.finite(dat$lon)&abs(dat$lat)<=90&abs(dat$lon)<=180&!is.na(dat$REALM)&!is.na(dat$sampling_site_id)&nzchar(dat$sampling_site_id)
dat$coordinate_valid<-coord_ok
dat$coordinate_site_id<-NA_character_
dat$coordinate_site_id[coord_ok]<-dat$sampling_site_id[coord_ok]
valid<-dat[coord_ok,]
stopifnot(nrow(valid)>0L, !anyNA(valid$coordinate_site_id))
gr<-split(seq_len(nrow(valid)),valid$coordinate_site_id)
stopifnot(all(vapply(gr,function(i)length(unique(valid$REALM[i]))==1L,logical(1))))
site<-do.call(rbind,lapply(gr,function(i)data.frame(coordinate_site_id=valid$coordinate_site_id[i[1]],REALM=valid$REALM[i[1]],realm_code=unname(codes[valid$REALM[i[1]]]),lat=median(valid$lat[i]),lon=median(valid$lon[i]),n_samples=length(i),q0=median(valid$q0[i]),q1=median(valid$q1[i]),q2=median(valid$q2[i]))))
row.names(site)<-NULL
site$REALM<-factor(site$REALM,levels=realm_order)
counts<-do.call(rbind,lapply(realm_order,function(r)data.frame(REALM=r,realm_code=unname(codes[r]),n_selected_metagenomes=sum(dat$REALM==r),n_missing_or_invalid_coordinates=sum(dat$REALM==r&!dat$coordinate_valid),n_coordinate_metagenomes=sum(valid$REALM==r),n_coordinate_sites=sum(site$REALM==r))))
counts$included_primary_inference<-counts$n_coordinate_sites>=5L
eligible<-counts$REALM[counts$included_primary_inference]
site$included_primary_inference<-as.character(site$REALM)%in%eligible
primary<-site[site$included_primary_inference,]
primary$REALM<-droplevels(primary$REALM)
dat$inference_status<-ifelse(!dat$coordinate_valid,'excluded_missing_or_invalid_coordinates',ifelse(dat$REALM%in%eligible,'included_primary_coordinate_site_inference','excluded_realm_fewer_than_5_coordinate_sites'))
stopifnot(nrow(primary)>nlevels(primary$REALM), nlevels(primary$REALM)>=2L)

dunn_holm<-function(values,groups) {
  groups<-droplevels(factor(groups));n<-length(values);r<-rank(values,ties.method='average');ties<-as.numeric(table(values));v<-n*(n+1)/12-sum(ties^3-ties)/(12*(n-1));nn<-table(groups);mr<-tapply(r,groups,mean)
  pp<-combn(levels(groups),2,simplify=FALSE)
  ans<-do.call(rbind,lapply(pp,function(p){z<-(mr[p[1]]-mr[p[2]])/sqrt(v*(1/nn[p[1]]+1/nn[p[2]]));data.frame(group_1=p[1],group_2=p[2],code_1=unname(codes[p[1]]),code_2=unname(codes[p[2]]),n_sites_1=unname(nn[p[1]]),n_sites_2=unname(nn[p[2]]),z=unname(z),p_raw=2*pnorm(-abs(z)))}))
  ans$p_holm<-p.adjust(ans$p_raw,'holm');ans$significant_005<-ans$p_holm<=.05;ans
}
omnibus<-list();pairwise<-list();desc<-list()
for(q in c('q0','q1','q2')) {
  kw<-kruskal.test(primary[[q]],primary$REALM)
  omnibus[[q]]<-data.frame(diversity=q,unit='median_per_resolved_sampling_site',n_sites=nrow(primary),n_realms=nlevels(primary$REALM),n_metagenomes_represented=sum(primary$n_samples),H=unname(kw$statistic),df=unname(kw$parameter),p_value=kw$p.value,epsilon_squared_H=(unname(kw$statistic)-nlevels(primary$REALM)+1)/(nrow(primary)-nlevels(primary$REALM)))
  pw<-dunn_holm(primary[[q]],primary$REALM);pw$diversity<-q;pairwise[[q]]<-pw
  desc[[q]]<-do.call(rbind,lapply(levels(primary$REALM),function(r){v<-primary[[q]][primary$REALM==r];data.frame(diversity=q,REALM=r,realm_code=unname(codes[r]),n_sites=length(v),n_metagenomes=sum(primary$n_samples[primary$REALM==r]),median=median(v),lower_quartile=unname(quantile(v,.25)),upper_quartile=unname(quantile(v,.75)),IQR=IQR(v))}))
}
omnibus<-do.call(rbind,omnibus);row.names(omnibus)<-NULL;omnibus$p_holm_across_three_q<-p.adjust(omnibus$p_value,'holm')
pairwise<-do.call(rbind,pairwise);row.names(pairwise)<-NULL
desc<-do.call(rbind,desc);row.names(desc)<-NULL
save_csv<-function(x,n)write.csv(x,file.path(out,n),row.names=FALSE,na='')
save_csv(omnibus,'alpha_coordinate_site_inference.csv')
save_csv(pairwise,'alpha_coordinate_site_Dunn_Holm.csv')
save_csv(desc,'alpha_coordinate_site_descriptive.csv')
save_csv(site,'alpha_coordinate_site_medians.csv')
save_csv(counts,'alpha_coordinate_site_cohort_counts.csv')
save_csv(dat,'alpha_coordinate_site_inference_sample_manifest.csv')
flow<-data.frame(step=c('selected_alpha_cohort','excluded_unresolved_sampling_site','resolved_site_cohort','excluded_realms_fewer_than_5_sites','primary_realms_at_least_5_sites'),n_metagenomes=c(nrow(dat),sum(!coord_ok),nrow(valid),sum(site$n_samples[!site$included_primary_inference]),sum(primary$n_samples)),n_coordinate_sites=c(NA,NA,nrow(site),sum(!site$included_primary_inference),nrow(primary)))
save_csv(flow,'alpha_coordinate_site_cohort_flow.csv')

# Metagenome-level summaries describe the displayed distributions. The tests
# above use coordinate-site medians and provide the primary realm inference.
ord <- read.csv(file.path(input_dir, '04a_realm_longitude_order_west_to_east.csv'))
ord <- ord[order(ord$west_to_east_position), ]
realm_levels <- ord$REALM
dat$REALM <- factor(dat$REALM, levels = realm_levels)
sites <- site
cols <- c('Central Indo-Pacific'='#66C2A5','Indo-Pacific'='#FC8D62','Temperate Northern Atlantic'='#8DA0CB','Tropical Atlantic'='#E78AC3','Tropical Eastern Pacific'='#A6D854','Western Indo-Pacific'='#FFD92F')
write_table <- function(d,n) write.csv(d,file.path(out_dir,n),row.names=FALSE,na='')
dunn_holm <- function(values,groups) {
  groups <- droplevels(factor(groups)); n<-length(values); r<-rank(values,ties.method='average')
  ties<-as.numeric(table(values)); v<-n*(n+1)/12-sum(ties^3-ties)/(12*(n-1))
  nn<-table(groups); mr<-tapply(r,groups,mean)
  pp<-combn(levels(groups),2,simplify=FALSE)
  ans<-do.call(rbind,lapply(pp,function(p) {z<-(mr[p[1]]-mr[p[2]])/sqrt(v*(1/nn[p[1]]+1/nn[p[2]])); data.frame(group_1=p[1],group_2=p[2],code_1=unname(codes[p[1]]),code_2=unname(codes[p[2]]),n_1=unname(nn[p[1]]),n_2=unname(nn[p[2]]),z=unname(z),p_raw=2*pnorm(-abs(z)))}))
  ans$p_holm<-p.adjust(ans$p_raw,'holm');ans$significant_005<-ans$p_holm<=.05;ans
}
compact_letters <- function(pw,group_order) {
  k<-length(group_order); pm<-matrix(1,k,k,dimnames=list(group_order,group_order))
  for(i in seq_len(nrow(pw))) {pm[pw$group_1[i],pw$group_2[i]]<-pw$p_holm[i];pm[pw$group_2[i],pw$group_1[i]]<-pw$p_holm[i]}
  nonsig<-pm>.05;diag(nonsig)<-TRUE
  cliques<-list()
  for(mask in seq_len(2^k-1L)) {members<-group_order[as.logical(intToBits(mask)[seq_len(k)])];if(all(nonsig[members,members,drop=FALSE]))cliques[[length(cliques)+1L]]<-members}
  maximal<-Filter(function(x)!any(vapply(cliques,function(y)length(y)>length(x)&&all(x%in%y),logical(1))),cliques)
  maximal<-maximal[order(vapply(maximal,function(x)min(match(x,group_order)),numeric(1)))]
  ans<-setNames(rep('',k),group_order)
  for(i in seq_along(maximal)) for(g in maximal[[i]]) ans[g]<-paste0(ans[g],letters[i])
  data.frame(REALM=names(ans),letters=unname(ans))
}
kw_list<-list(); pw_list<-list(); ds_list<-list(); letters_list<-list()
for(unit in 'metagenome') {
  dd<-dat
  for(q in c('q0','q1','q2')) {
    kw<-kruskal.test(dd[[q]],dd$REALM)
    kw_list[[paste(unit,q)]]<-data.frame(unit=unit,diversity=q,n=nrow(dd),n_realms=nlevels(dd$REALM),H=unname(kw$statistic),df=unname(kw$parameter),p_value=kw$p.value,epsilon_squared_H=(unname(kw$statistic)-nlevels(dd$REALM)+1)/(nrow(dd)-nlevels(dd$REALM)))
    pw<-dunn_holm(dd[[q]],dd$REALM); pw$unit<-unit;pw$diversity<-q
    pw_list[[paste(unit,q)]]<-pw
    ds<-do.call(rbind,lapply(realm_levels,function(r) {v<-dd[[q]][dd$REALM==r];data.frame(unit=unit,diversity=q,REALM=r,realm_code=unname(codes[r]),n_units=length(v),n_samples=sum(dat$REALM==r),n_sites=sum(sites$REALM==r),median=median(v),lower_quartile=unname(quantile(v,.25)),upper_quartile=unname(quantile(v,.75)),IQR=IQR(v),minimum=min(v),maximum=max(v))}))
    cl<-compact_letters(pw,ds$REALM[order(-ds$median)]);cl$unit<-unit;cl$diversity<-q
    letters_list[[paste(unit,q)]]<-cl
    ds$letters<-cl$letters[match(ds$REALM,cl$REALM)];ds_list[[paste(unit,q)]]<-ds
  }
}
kw_all<-do.call(rbind,kw_list);row.names(kw_all)<-NULL
kw_all$p_holm_across_q_within_unit<-ave(kw_all$p_value,kw_all$unit,FUN=function(p)p.adjust(p,'holm'))
pw_all<-do.call(rbind,pw_list);row.names(pw_all)<-NULL
ds_all<-do.call(rbind,ds_list);row.names(ds_all)<-NULL
cl_all<-do.call(rbind,letters_list);row.names(cl_all)<-NULL
write_table(kw_all,'alpha_Kruskal_Wallis_metagenome.csv')
write_table(pw_all,'alpha_Dunn_Holm_metagenome.csv')
write_table(ds_all,'alpha_descriptive_metagenome.csv')
write_table(cl_all,'alpha_compact_letters_metagenome.csv')
GeomHalfViolin<-ggproto('GeomHalfViolin',GeomViolin,draw_group=function(self,data,...,draw_quantiles=NULL){data<-transform(data,xminv=x-violinwidth*(x-xmin),xmaxv=x);p<-rbind(transform(data,x=xminv)[order(data$y),],transform(data,x=xmaxv)[order(data$y,decreasing=TRUE),]);p<-rbind(p,p[1,]);GeomPolygon$draw_panel(p,...)})
half_violin<-function(...) layer(stat='ydensity',geom=GeomHalfViolin,position='identity',params=list(...))
set.seed(3030) # Reproducible point jitter, shared with the q1 figure
dat$realm_index<-as.numeric(dat$REALM)
dat$point_x<-dat$realm_index+runif(nrow(dat),.10,.28)
make_panel<-function(q,tag) {
  dd<-dat;dd$value<-dd[[q]]
  ds<-ds_all[ds_all$unit=='metagenome'&ds_all$diversity==q,];ds$x<-match(ds$REALM,realm_levels)
  ymax<-max(dd$value);ds$letter_y<-ds$maximum+.040*ymax
  k<-kw_all[kw_all$unit=='metagenome'&kw_all$diversity==q,]
  ylab<-if(q=='q0')'Hill diversity (q = 0)' else 'Hill diversity (q = 2)'
  ggplot(dd,aes(x=realm_index,y=value,group=REALM,fill=REALM))+
    half_violin(width=.72,trim=TRUE,scale='width',alpha=.55,colour=NA)+
    geom_boxplot(width=.15,outlier.shape=NA,alpha=.82,colour='#333333',linewidth=.38)+
    geom_point(aes(x=point_x,y=value,fill=REALM),inherit.aes=FALSE,shape=21,size=1.60,stroke=.22,colour='white',alpha=.86)+
    geom_text(data=ds,aes(x=x,y=letter_y,label=letters),inherit.aes=FALSE,fontface='bold',size=3.8,colour='#222222')+
    annotate('text',x=length(realm_levels)+.32,y=1.155*ymax,label=sprintf('Kruskal–Wallis: H(%d) = %.2f, %s',k$df,k$H,if(k$p_value<.001)'P < 0.001' else sprintf('P = %.3f',k$p_value)),hjust=1,vjust=1,size=3.4,colour='#222222')+
    scale_fill_manual(values=cols,drop=FALSE)+
    scale_x_continuous(breaks=seq_along(realm_levels),labels=paste0(unname(codes[realm_levels]),'\nn = ',as.numeric(table(dat$REALM))),limits=c(.48,length(realm_levels)+.48),expand=expansion(mult=0))+
    scale_y_continuous(limits=c(0,1.21*ymax),expand=expansion(mult=0),labels=scales::label_comma())+
    labs(tag=tag,x='Marine realms ordered west to east',y=ylab)+coord_cartesian(clip='off')+
    theme_minimal(base_size=11.5,base_family='sans')+
    theme(legend.position='none',panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),panel.grid.major.y=element_line(colour='#E8E8E8',linewidth=.35),axis.line=element_line(colour='#333333',linewidth=.42),axis.ticks=element_line(colour='#333333',linewidth=.35),axis.text.x=element_text(colour='#222222',size=10.5,lineheight=1.05,margin=margin(t=6)),axis.text.y=element_text(colour='#333333',size=10.5),axis.title.x=element_text(size=11,margin=margin(t=8)),axis.title.y=element_text(size=11.5,margin=margin(r=8)),plot.tag=element_text(face='bold',size=15),plot.tag.position=c(.006,.99),plot.margin=margin(t=8,r=10,b=6,l=8),plot.background=element_rect(fill='white',colour=NA))
}
figure<-make_panel('q0','A')/make_panel('q2','B')
ggsave(file.path(out_dir,'Figure_S_Hill_q0_q2_by_realm.png'),figure,width=7.2,height=8.6,units='in',dpi=450,bg='white')
ggsave(file.path(out_dir,'Figure_S_Hill_q0_q2_by_realm.pdf'),figure,width=7.2,height=8.6,units='in',device=cairo_pdf,bg='white')
capture.output(sessionInfo(),file=file.path(out_dir,'sessionInfo.txt'))
writeLines(c(
  "Primary inference uses one median per resolved sampling site, as recorded in sampling_site_id.",
  "Coordinates define sites unless source records identify repeated sampling of one locality.",
  sprintf("%d realms with at least five sites contribute %d sites to the tests.",nlevels(primary$REALM),nrow(primary)),
  sprintf("All %d observed realms are retained in metagenome-level descriptive plots.",length(realm_levels)),
  "Omnibus P values are Holm-adjusted across q0, q1 and q2; pairwise Dunn P values are Holm-adjusted within each diversity order.",
  "No additional rarefaction is applied in this stage."
), file.path(out, "analysis_notes.txt"))
message("Saved Hill diversity supplementary analyses in ", out)
