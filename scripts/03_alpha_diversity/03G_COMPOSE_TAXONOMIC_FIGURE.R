#!/usr/bin/env Rscript
# Compose Figure 2 from completed analyses and a separately supplied iCAMP summary.
# No analysis is estimated here. iCAMP contributions refer to their own source cohort.
options(stringsAsFactors = FALSE, bitmapType = 'cairo')
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(grid)})
root <- Sys.getenv('BIOGEOGRAPHY_BASE', unset='.')
results <- Sys.getenv('BIOGEOGRAPHY_RESULTS', unset=file.path(root,'results'))
out <- Sys.getenv('FIGURE2_OUTPUT_DIR',unset=results)
alpha_dir <- Sys.getenv('ALPHA_INPUT_DIR',unset=file.path(results,'alpha'))
gam_dir <- Sys.getenv('GAM_INPUT_DIR',unset=file.path(results,'gam'))
decay_dir <- Sys.getenv('DECAY_INPUT_DIR',unset=file.path(results,'decay'))
nmds_dir <- Sys.getenv('NMDS_INPUT_DIR',unset=file.path(results,'nmds'))
icamp_file <- Sys.getenv('ICAMP_SUMMARY_TSV',unset=file.path(results,'icamp','icamp_process_contributions.tsv'))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
inputs <- c(file.path(alpha_dir,'Figure_2a_ggplot.rds'),file.path(gam_dir,'tables','17A_all_environmental_partial_effect_curves.csv'),file.path(gam_dir,'tables','17A_full_model_input.csv'),file.path(gam_dir,'tables','17A_selected_top4_effects.csv'),file.path(decay_dir,'rds','06_Figure_2c_ggplot_object.rds'),file.path(nmds_dir,'rds','05_Figure_2d_ggplot_object.rds'),icamp_file)
if(any(!file.exists(inputs))) stop('Missing figure inputs: ',paste(inputs[!file.exists(inputs)],collapse=', '),call.=FALSE)
cols <- c('Central Indo-Pacific'='#66C2A5','Temperate Northern Atlantic'='#8DA0CB','Tropical Atlantic'='#E78AC3','Tropical Eastern Pacific'='#A6D854','Western Indo-Pacific'='#FFD92F')
codes <- c('Central Indo-Pacific'='CIP','Temperate Northern Atlantic'='TNA','Tropical Atlantic'='TA','Tropical Eastern Pacific'='TEP','Western Indo-Pacific'='WIP')
common <- theme(axis.title=element_text(size=10),axis.text=element_text(size=9,colour='black'),plot.title=element_blank(),plot.tag=element_blank(),plot.margin=margin(7,7,7,7),plot.background=element_rect(fill='white',colour=NA))
pa <- readRDS(inputs[1]) + labs(tag=NULL,title=NULL,x='Marine realm',y='Hill diversity (q = 1)') + common
# Main-panel annotations are kept in supplementary tables; letters remain on a.
pa$layers <- pa$layers[-5L]
curves <- read.csv(inputs[2]); model_input <- read.csv(inputs[3]); top <- read.csv(inputs[4]);top<-top[order(top$effect_rank),]
xlabels <- c(bio15='Precipitation seasonality\n(BIO15)',soc_mean='Soil organic carbon\n(dg/kg)',ph_mean='Soil pH × 10',bio13='Precipitation of wettest\nmonth (BIO13)',bio2='Mean diurnal temperature\nrange (BIO2)',bio3='Isothermality (BIO3)')
lim<-range(curves$lower[curves$variable%in%top$variable],curves$upper[curves$variable%in%top$variable]);lim<-lim+c(-.06,.08)*diff(lim)
b_panels<-lapply(seq_len(nrow(top)),function(i){v<-top$variable[i];d<-curves[curves$variable==v,];rug<-data.frame(x=model_input[[v]]);rug<-rug[rug$x>=min(d$x)&rug$x<=max(d$x),,drop=FALSE];label<-if(v%in%names(xlabels))xlabels[[v]] else top$panel_label[i]
 ggplot(d,aes(x=x,y=effect))+geom_hline(yintercept=0,linetype='dashed',linewidth=.25,colour='#888888')+geom_ribbon(aes(ymin=lower,ymax=upper),fill='#56B4E9',alpha=.28)+geom_line(colour='#0072B2',linewidth=.55)+geom_rug(data=rug,aes(x=x),inherit.aes=FALSE,linewidth=.14,alpha=.35)+scale_x_continuous(breaks=scales::breaks_pretty(n=3))+scale_y_continuous(breaks=c(-.4,0,.4))+coord_cartesian(ylim=lim)+labs(x=label,y=if(i%%2==1)'Partial effect' else NULL)+theme_classic(base_size=9)+theme(axis.title=element_text(size=9),axis.text=element_text(size=8,colour='black'),axis.line=element_line(linewidth=.3),plot.margin=margin(4,6,3,4))
})
pb <- wrap_plots(b_panels,ncol=2)
pc <- readRDS(inputs[5]);pc$layers<-pc$layers[1:2]
pc <- pc + labs(tag=NULL) + scale_y_continuous(name='Community similarity',limits=c(0,1),breaks=c(0,.25,.5,.75,1),expand=expansion(mult=c(0,.015))) + scale_x_continuous(name='Geographic distance (km)',breaks=log10(c(1,10,100,1000,10000)),labels=scales::label_comma()(c(1,10,100,1000,10000)),expand=expansion(mult=c(.025,.025))) + guides(colour=guide_legend(title=NULL,override.aes=list(alpha=1,size=2))) + common + theme(legend.position='bottom',legend.text=element_text(size=9),legend.margin=margin(0,0,0,0),legend.key.width=unit(.45,'cm'))
pd <- readRDS(inputs[6]) + labs(tag=NULL) + scale_fill_manual(values=cols,labels=codes,drop=TRUE,name=NULL) + guides(fill=guide_legend(nrow=1,override.aes=list(size=2.2,alpha=1))) + common + theme(legend.position='bottom',legend.text=element_text(size=8.5),legend.key.width=unit(.25,'cm'),legend.key.height=unit(.25,'cm'),legend.spacing.x=unit(.1,'cm'),legend.margin=margin(0,0,0,0))
processes <- c('heterogeneous_selection','homogeneous_selection','dispersal_limitation','homogenizing_dispersal','drift_weak_selection_or_undominated')
ic <- read.delim(icamp_file,check.names=FALSE);ic<-ic[match(processes,ic$process),]
stopifnot(nrow(ic)==5,!anyNA(ic$process),abs(sum(ic$mean_fraction)-1)<1e-7)
ic$process<-factor(ic$process,levels=processes);ic$xmax<-cumsum(ic$mean_percent);ic$xmin<-c(0,head(ic$xmax,-1));ic$mid<-(ic$xmin+ic$xmax)/2
process_labels<-c('Heterogeneous selection','Homogeneous selection','Dispersal limitation','Homogenizing dispersal','Drift, weak selection or\nother undominated processes')
process_cols<-setNames(c('#CC79A7','#E69F00','#56B4E9','#009E73','#B3B3B3'),processes)
pe<-ggplot(ic)+geom_rect(aes(xmin=xmin,xmax=xmax,ymin=0,ymax=1,fill=process),colour='white',linewidth=.4)+geom_text(aes(x=mid,y=.5,label=sprintf('%.1f%%',mean_percent)),size=3.0)+scale_fill_manual(values=process_cols,labels=process_labels,drop=FALSE)+scale_x_continuous(limits=c(0,100),expand=expansion(mult=0))+scale_y_continuous(limits=c(0,1),expand=expansion(mult=0))+guides(fill=guide_legend(nrow=2,byrow=TRUE,title=NULL))+theme_void(base_size=9.5)+theme(legend.position='bottom',legend.text=element_text(size=9),legend.key.size=unit(.23,'cm'),legend.spacing.x=unit(.15,'cm'),legend.margin=margin(3,0,0,0),plot.margin=margin(3,7,3,7))
# Explicit viewports give the four GAM curves room without changing data geometry.
draw <- function(){grid.newpage();print(pa,vp=viewport(x=.195,y=.79,width=.38,height=.40));print(pb,vp=viewport(x=.69,y=.79,width=.60,height=.40));print(pc,vp=viewport(x=.255,y=.365,width=.49,height=.40));print(pd,vp=viewport(x=.745,y=.365,width=.49,height=.40));print(pe,vp=viewport(x=.505,y=.080,width=.96,height=.145));for(z in list(c('a',.013,.982),c('b',.402,.982),c('c',.013,.565),c('d',.505,.565),c('e',.013,.162)))grid.text(z[1],x=as.numeric(z[2]),y=as.numeric(z[3]),just=c('left','top'),gp=gpar(fontfamily='sans',fontface='bold',fontsize=13))}
png(file.path(out,'Figure_2.png'),width=9,height=7.5,units='in',res=600,bg='white');draw();dev.off()
cairo_pdf(file.path(out,'Figure_2.pdf'),width=9,height=7.5,bg='white',onefile=TRUE);draw();dev.off()
icamp_out <- file.path(results, 'icamp')
dir.create(icamp_out, recursive=TRUE, showWarnings=FALSE)
ggsave(file.path(icamp_out,'Figure_2e_icamp.png'),pe,width=8,height=1.7,dpi=600,bg='white')
ggsave(file.path(icamp_out,'Figure_2e_icamp.pdf'),pe,width=8,height=1.7,device=cairo_pdf,bg='white')
write.csv(data.frame(panel=c('a','b','b','b','c','d','e'),input=inputs),file.path(out,'Figure_2_source_inputs.csv'),row.names=FALSE)
writeLines(c('Figure 2 | Geographic patterns, environmental associations and community assembly of mangrove-associated microbiomes.', '(a) Hill diversity (q = 1) across five marine realms (341 profiles). Half-violins, boxes and points show distributions, medians and interquartile ranges, and individual profiles. Letters summarize descriptive metagenome-level Dunn–Holm comparisons; primary tests use independent site medians.', '(b) Four environmental smooths with the largest fitted amplitudes in the Gaussian GAM with 14 environmental smooths and a locality random effect (274 profiles, 33 localities; 73.2% deviance explained). Lines and bands show centered partial effects on log1p(Shannon entropy) and pointwise 95% intervals; rugs indicate observations. Curves span the central 95% of predictor values; strong predictor overlap limits individual attribution.', '(c) Bray–Curtis and Sørensen similarity versus geographic distance across 323 profiles from 122 sampling sites. Lines show descriptive linear fits.', '(d) Aitchison NMDS of 323 profiles, colored by marine realm; contours show a smoothed longitude surface (stress = 0.054).', '(e) Mean abundance-weighted iCAMP contributions from the separate 400-metagenome analysis (1,892 genera; 24,692 within-realm pairs; 999 randomizations). These retained estimates were not recalculated for the 341-profile reference dataset.'),file.path(out,'Figure_2_caption.txt'))
cat('Figure 2 saved in ',normalizePath(out),'\n',sep='')
