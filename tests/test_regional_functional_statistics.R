#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(vegan))
args <- commandArgs(trailingOnly=TRUE)
script <- if(length(args)) args[1] else "scripts/09_functional_regional_differentiation/09A_FUNCTIONAL_PREVALENCE_REGIONAL_DIFFERENTIATION.R"
needed <- c("numeric_matrix_from_df", "permanova_stats_fast")
for(expr in parse(script)) {
  if(is.call(expr) && identical(expr[[1]], as.name("<-")) && is.symbol(expr[[2]]) && as.character(expr[[2]]) %in% needed) eval(expr)
}
set.seed(431)
y <- matrix(rnorm(18*30),18,30)
g <- factor(rep(c("A","B","C"),each=6))
y[g=="B",1:4] <- y[g=="B",1:4]+1
st <- permanova_stats_fast(y,g)
ad <- adonis2(dist(y)~g,permutations=0)
stopifnot(abs(st["R2"]-ad$R2[1])<1e-12,abs(st["F"]-ad$F[1])<1e-12)
z <- seq_len(18)+rnorm(18)
d2 <- as.matrix(dist(y))^2
G <- -0.5*(d2-rowMeans(d2)-rep(colMeans(d2),each=18)+mean(d2))
q0 <- qr.Q(qr(cbind(1,z)));q1 <- qr.Q(qr(model.matrix(~z+g)))
unique_R2 <- (sum(diag(crossprod(q1,G%*%q1)))-sum(diag(crossprod(q0,G%*%q0))))/sum(diag(G))
adj <- adonis2(dist(y)~z+g,by="margin",permutations=19)
stopifnot(abs(unique_R2-adj["g","R2"])<1e-12)
for(bad in list(data.frame(x=c(1,NA)),data.frame(x=c(1,-1)),data.frame(x=c(1,"error")))) {
  stopifnot(inherits(try(numeric_matrix_from_df(bad),silent=TRUE),"try-error"))
}
cat("PASS: Euclidean R2/F, adjusted term SS, and missing/negative/non-numeric count rejection.\n")
