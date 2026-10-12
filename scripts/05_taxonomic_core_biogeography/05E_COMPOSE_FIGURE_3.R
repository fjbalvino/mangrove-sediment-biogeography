#!/usr/bin/env Rscript
# Assemble Figure 3 from the final analysis plot objects; no statistics are refitted.
suppressPackageStartupMessages({library(ggplot2);library(patchwork)})
results <- Sys.getenv('BIOGEOGRAPHY_RESULTS', unset='results/taxonomic_core')
read_run <- function(tag) {
  pointer <- file.path(results, paste0('LATEST_', tag, '.txt'))
  stopifnot(file.exists(pointer))
  path <- readLines(pointer, warn=FALSE)[[1]]
  stopifnot(dir.exists(path))
  path
}
a_run <- read_run('05A_TAXONOMIC_CORE_ASSOCIATED')
b_run <- read_run('05C_CORE_ASSOCIATED_PHYLA')
decay_run <- read_run('05D_CORE_NONCORE_DISTANCE_DECAY')
a <- readRDS(file.path(a_run,'rds/09A_Figure_3a_analysis_objects.rds'))$plot
b <- readRDS(file.path(b_run,'rds/11_Figure_3b_analysis_objects.rds'))$plot
decay <- readRDS(file.path(decay_run,'rds/12_Figure_3cd_analysis_objects.rds'))
stopifnot('Reads in eligible taxa' %in% as.character(a$data$metric))

prepare_decay <- function(plot,result,x_label) {
  plot$layers <- Filter(function(layer) !inherits(layer$geom,'GeomText'),plot$layers)
  # Clip the displayed response range after fitting, retaining every observation.
  y_scale <- plot$scales$get_scales('y')
  y_scale$limits <- NULL
  y_scale$name <- 'Community similarity\n(1 − dissimilarity)'
  plot <- plot + coord_cartesian(ylim=c(0,1))
  s <- result$statistics
  core <- s[s$fraction=='Core-associated',]
  non <- s[s$fraction=='Non-core',]
  contrast <- result$contrast
  subtitle <- paste(
    sprintf('Core: r = %.3f; slope = %.3f; P = %.4f',core$mantel_r,core$similarity_slope,core$permutation_p),
    sprintf('Non-core: r = %.3f; slope = %.3f; P = %.4f',non$mantel_r,non$similarity_slope,non$permutation_p),
    sprintf('Δ slope = %.3f; P = %.4f',contrast$delta_slope,contrast$permutation_p),
    sep='\n')
  plot + labs(x=x_label,y='Community similarity\n(1 − dissimilarity)',subtitle=subtitle) +
    theme(plot.subtitle=element_text(size=14,lineheight=1.1,margin=margin(b=10)),
          axis.title=element_text(size=15),axis.text=element_text(size=14),
          legend.text=element_text(size=14),plot.margin=margin(8,24,8,10))
}
c <- prepare_decay(decay$geographic$plot,decay$geographic$result,'Geographic distance (km, log10 scale)')
d <- prepare_decay(decay$environmental$plot,decay$environmental$result,'Climatic distance')
a <- a + coord_cartesian(clip='off') + theme(axis.text=element_text(size=14,colour='black'),
               axis.title=element_text(size=15),
               legend.text=element_text(size=14),plot.margin=margin(8,24,8,10))
b <- b + coord_cartesian(clip='off') + labs(x='Mean percentage of core-associated reads') +
  guides(fill=guide_legend(nrow=3,ncol=2,byrow=TRUE)) +
  theme(axis.text=element_text(size=14,colour='black'),
        axis.title=element_text(size=15),
        legend.text=element_text(size=14),plot.margin=margin(8,24,8,10))
# Explicit child elements override sizes inherited from the source plots.
a <- a + theme(axis.text.x=element_text(size=14),axis.text.y=element_text(size=14),axis.title.x=element_text(size=15))
b <- b + theme(axis.text.x=element_text(size=14),axis.text.y=element_text(size=14),axis.title.x=element_text(size=15))
top <- a | b
bottom <- (c | d) + plot_layout(guides='collect') & theme(legend.position='bottom')
p <- (top / bottom) + plot_layout(heights=c(0.68,1)) +
  plot_annotation(tag_levels='a')
p <- p & theme(plot.tag=element_text(face='bold',size=18))
ggsave(file.path(results,'Figure_3.png'),p,width=13.2,height=9.4,units='in',dpi=450,bg='white')
ggsave(file.path(results,'Figure_3.pdf'),p,width=13.2,height=9.4,units='in',device=cairo_pdf,bg='white')
saveRDS(p,file.path(results,'Figure_3_plot.rds'))
writeLines(c(a_run,b_run,decay_run),file.path(results,'Figure_3_sources.txt'))
cat('FIGURE_3_COMPOSITE_PASS\n')
