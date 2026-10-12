#!/usr/bin/env Rscript
# Plot and tabulate an existing GAM; this stage does not refit the model.
options(stringsAsFactors = FALSE, bitmapType = "cairo")
suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(mgcv)})
root <- Sys.getenv("BIOGEOGRAPHY_BASE", unset = ".")
results <- Sys.getenv("BIOGEOGRAPHY_RESULTS", unset = file.path(root, "results"))
src <- Sys.getenv("GAM_INPUT_DIR", unset = "")
if (!nzchar(src)) {
  pointer <- file.path(results, "LATEST_SHANNON_GAM.txt")
  if (!file.exists(pointer)) stop("Run 03B first or supply GAM_INPUT_DIR.", call. = FALSE)
  src <- trimws(readLines(pointer, n = 1L, warn = FALSE))
}
out <- Sys.getenv("GAM_FIGURE_OUTPUT_DIR", unset = file.path(src, "manuscript"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
read_tab <- function(n) read.csv(file.path(src, "tables", n), check.names = FALSE)
curves <- read_tab("17A_all_environmental_partial_effect_curves.csv")
ranking <- read_tab("17A_environmental_effect_ranking.csv")
input <- read_tab("17A_full_model_input.csv")
stats <- read_tab("17A_all_smooth_terms.csv")
metrics <- read_tab("17A_model_metrics.csv")
kcheck <- read_tab("17A_k_check.csv")
conc <- read_tab("17A_concurvity.csv")
stopifnot(nrow(input)==metrics$n_after_4sd_rule, nrow(ranking)==14, metrics$convergence_accepted)
stats$variable <- sub("\\)$", "", sub("^s\\(", "", stats$term))
stats$BH_q_14_environmental <- NA_real_
env <- stats$variable != "Locality"
stats$BH_q_14_environmental[env] <- p.adjust(stats$`p-value`[env], method = "BH")
stats$Predictor <- ranking$panel_label[match(stats$variable, ranking$variable)]
stats$Predictor[!env] <- "Sampling locality (random effect)"
stats$effect_amplitude_central95 <- ranking$central_95_effect_amplitude[match(stats$variable,ranking$variable)]
stats <- stats[order(match(stats$variable,c(ranking$variable,"Locality"))),]
write.csv(stats,file.path(out,"GAM_smooth_statistics_full_precision.csv"),row.names=FALSE)
fmt_n <- function(x) ifelse(is.na(x),"—",ifelse(abs(x)<.001,"<0.001",sprintf("%.3f",x)))
fmt_p <- function(x) ifelse(is.na(x),"—",ifelse(x<.001,"<0.001",sprintf("%.3f",x)))
display <- data.frame(Predictor=stats$Predictor, edf=fmt_n(stats$edf),
                     `Ref. df`=fmt_n(stats$Ref.df), F=fmt_n(stats$F),
                     P=fmt_p(stats$`p-value`), `BH q`=fmt_p(stats$BH_q_14_environmental),
                     check.names=FALSE)
write.csv(display,file.path(out,"GAM_smooth_statistics_display.csv"),row.names=FALSE)
diagnostics <- data.frame(Predictor=stats$Predictor,
  `k prime`=kcheck$`k'`[match(stats$term,kcheck$term)],
  edf=stats$edf,
  `k index`=kcheck$`k-index`[match(stats$term,kcheck$term)],
  `k P`=kcheck$`p-value`[match(stats$term,kcheck$term)],
  `Concurvity worst`=as.numeric(conc[conc$metric=="worst",stats$term]),
  `Concurvity observed`=as.numeric(conc[conc$metric=="observed",stats$term]),
  `Concurvity estimate`=as.numeric(conc[conc$metric=="estimate",stats$term]),check.names=FALSE)
write.csv(diagnostics,file.path(out,"GAM_diagnostics_full_precision.csv"),row.names=FALSE)
dd <- diagnostics
dd$edf <- fmt_n(dd$edf)
for(n in c("k index","Concurvity worst","Concurvity observed","Concurvity estimate")) dd[[n]] <- ifelse(is.na(dd[[n]]),"—",sprintf("%.3f",dd[[n]]))
dd$`k P` <- ifelse(is.na(dd$`k P`),"—",ifelse(dd$`k P`==0,"<0.0025",format(dd$`k P`,digits=3,trim=TRUE)))
write.csv(dd,file.path(out,"GAM_diagnostics_display.csv"),row.names=FALSE)
write.csv(metrics,file.path(out,"GAM_model_metrics.csv"),row.names=FALSE)

# Audit the saved curves against the saved object (read/predict only).
model <- readRDS(file.path(src,"models","17A_full_environmental_GAM_Shannon.rds"))
vars <- ranking$variable
validation <- lapply(vars,function(v){
  cv <- curves[curves$variable==v,]
  sel <- c(1,ceiling(nrow(cv)/2),nrow(cv))
  nd <- input[rep(1,3),]
  for(w in ranking$variable) nd[[w]] <- median(input[[w]])
  nd$Locality <- factor(levels(model$model$Locality)[1],levels=levels(model$model$Locality))
  nd[[v]] <- cv$x[sel]
  pred <- predict(model,newdata=nd,type="terms",se.fit=TRUE,exclude="s(Locality)")
  tm <- paste0("s(",v,")")
  data.frame(variable=v,max_abs_effect_difference=max(abs(pred$fit[,tm]-cv$effect[sel])),
    max_abs_se_difference=max(abs(pred$se.fit[,tm]-cv$standard_error[sel])),
    max_abs_interval_error=max(abs(cv$lower-cv$effect+1.96*cv$standard_error),abs(cv$upper-cv$effect-1.96*cv$standard_error)),
    min_x=cv$x[1],max_x=cv$x[nrow(cv)])
})
validation <- do.call(rbind,validation)
stopifnot(max(validation$max_abs_effect_difference)<1e-10,max(validation$max_abs_se_difference)<1e-10,max(validation$max_abs_interval_error)<1e-10)
write.csv(validation,file.path(out,"GAM_curve_validation.csv"),row.names=FALSE)

xlabels <- c(silt_mean="Soil silt content",ph_mean="Soil pH × 10",
 gbif_richness_species_2km="GBIF species richness\n(2 km)",
 bio15="Precipitation seasonality\n(BIO15)",bio3="Isothermality\n(BIO3)",
 bio13="Precipitation of wettest month\n(BIO13)",bio2="Mean diurnal temperature range\n(BIO2)")
fallback_labels <- setNames(ranking$panel_label, ranking$variable)
xlabels <- c(xlabels, fallback_labels[!names(fallback_labels) %in% names(xlabels)])
make_figure <- function(variables,ncol,height,filebase) {
  cvs <- curves[curves$variable %in% variables,]
  limits <- range(c(cvs$lower,cvs$upper))
  span <- diff(limits)
  limits <- limits+c(-.08,.24)*span
  ps <- lapply(seq_along(variables),function(i){
    v <- variables[i]; cv <- curves[curves$variable==v,]
    st <- stats[stats$variable==v,]
    rug <- data.frame(x=input[[v]])
    rug <- rug[rug$x>=min(cv$x)&rug$x<=max(cv$x),,drop=FALSE]
    ann <- paste0("F = ",sprintf("%.2f",st$F),"\n",if(st$`p-value`<.001) "P < 0.001" else paste0("P = ",sprintf("%.3f",st$`p-value`)))
    ggplot(cv,aes(x=x,y=effect))+
      geom_hline(yintercept=0,linewidth=.32,linetype="dashed",colour="#777777")+
      geom_ribbon(aes(ymin=lower,ymax=upper),fill="#56B4E9",alpha=.27)+
      geom_line(colour="#0072B2",linewidth=.65)+
      geom_rug(data=rug,aes(x=x),inherit.aes=FALSE,sides="b",alpha=.3,linewidth=.2)+
      annotate("text",x=Inf,y=Inf,label=ann,hjust=1.06,vjust=1.1,size=2.8,lineheight=1.08)+
      coord_cartesian(ylim=limits,expand=TRUE)+
      scale_y_continuous(breaks=scales::breaks_pretty(n=4))+
      scale_x_continuous(breaks=scales::breaks_pretty(n=3))+
      labs(x=unname(xlabels[v]),y=if ((i-1)%%ncol==0) "Partial effect on\nlog1p(Shannon entropy)" else NULL,tag=letters[i])+
      theme_classic(base_size=9.5,base_family="sans")+
      theme(axis.text=element_text(size=8.6,colour="black"),axis.title=element_text(size=9.3),
            axis.title.y=element_text(margin=margin(r=5)),axis.title.x=element_text(margin=margin(t=6)),
            plot.tag=element_text(face="bold",size=11),plot.tag.position=c(.02,.97),
            plot.margin=margin(t=11,r=7,b=5,l=7),axis.line=element_line(linewidth=.35),
            axis.ticks=element_line(linewidth=.3))
  })
  fig <- wrap_plots(ps,ncol=ncol)
  ggsave(file.path(out,paste0(filebase,".png")),fig,width=6.8,height=height,dpi=400,units="in",bg="white")
  ggsave(file.path(out,paste0(filebase,".pdf")),fig,width=6.8,height=height,units="in",device=cairo_pdf,bg="white")
}
top4 <- ranking$variable[order(ranking$effect_rank)][seq_len(4L)]
additional <- stats$variable[!is.na(stats$BH_q_14_environmental) & stats$BH_q_14_environmental < .05 & !stats$variable %in% top4]
if (GAM_PLOT_SCOPE == "top4") make_figure(top4, 2, 5.4, "Figure_2b_GAM_top4")
if (GAM_PLOT_SCOPE == "additional") {
  if (length(additional)) make_figure(additional, min(3L, length(additional)), 3.15 * ceiling(length(additional)/3), "Figure_S_GAM_additional_effects")
  else message("No additional environmental smooths meet BH-adjusted P < 0.05.")
}
writeLines(c("Partial-effect curves are conditional on all model terms.",
             "Four main-figure predictors are ranked by fitted amplitude across their central 95% range.",
             "Additional terms are selected by BH-adjusted P < 0.05 across 14 environmental smooths.",
             "Numerical convergence does not establish model adequacy; review concurvity and k-check diagnostics."),
           file.path(out, "GAM_interpretation.txt"))
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
cat("Saved figures, statistics and curve audit in",out,"\n")
