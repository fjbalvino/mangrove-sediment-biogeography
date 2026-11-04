# Conditional-design calibration; no real KO association tests.
args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(args[1]); mode <- if(length(args)>1) args[2] else 'final'
out <- file.path(root,'calibration_final')
library(clubSandwich)
g <- read.delim(file.path(out,'design.tsv'),check.names=FALSE)
u <- read.delim(file.path(out,'run_weights.tsv'),check.names=FALSE)
alpha <- as.matrix(read.delim(file.path(out,'detection_intercepts.tsv'),row.names=1,check.names=FALSE))
mu <- as.matrix(read.delim(file.path(out,'positive_intercepts.tsv'),row.names=1,check.names=FALSE))
K <- ncol(alpha); N <- nrow(g); U <- nrow(u)
stopifnot(K>0,N>0,U>0,identical(g$analysis_group,rownames(alpha)),identical(rownames(alpha),rownames(mu)))
g$study <- factor(g$study_id)
X <- as.matrix(g[c('x1','x2','x3')]); run.group <- match(u$analysis_group,g$analysis_group)
W <- matrix(0,N,U);W[cbind(run.group,seq_len(U))] <- u$descriptive_run_weight_within_group
stopifnot(max(abs(rowSums(W)-1))<1e-12)
markers <- readLines(file.path(out,'signal_genes.txt'))
dominant <- readLines(file.path(out,'dominant_study.txt'))
stopifnot(all(markers %in% colnames(alpha)))
is.signal <- colnames(alpha) %in% markers
spatial.chol <- function(lat,lon) {
 a <- lat*pi/180;b <- lon*pi/180
 h <- outer(a,a,function(x,y)sin((x-y)/2)^2) + outer(cos(a),cos(a))*outer(b,b,function(x,y)sin((x-y)/2)^2)
 dist <- 6371*2*asin(sqrt(pmin(1,pmax(0,h))))
 C <- .5*diag(length(a))+.5*exp(-dist/50)
 t(chol(C))
}
Lgroup <- spatial.chol(g$lat,g$lon)
Lrun <- spatial.chol(g$lat[run.group],g$lon[run.group])
studies <- levels(g$study); si <- as.integer(g$study)
normal.genes <- function(n) sqrt(.4)*matrix(rnorm(n),n,K)+sqrt(.6)*matrix(rnorm(n*K),n,K)
idx <- expand.grid(predictor=c('x1','x2','x3'),method=c('conventional','CR2'),component=c('detection','positive'),gene=colnames(alpha),stringsAsFactors=FALSE)
J <- nrow(idx);stopifnot(J==K*12)
fit.component <- function(y,component) {
 ans <- matrix(NA_real_,6,4,dimnames=list(NULL,c('beta','se','df','p')))
 ok <- is.finite(y)
 counts <- table(g$study[ok]); active <- names(counts)[counts>=2]
 if(component=='detection') {
  bounds <- tapply(y,g$study,function(v)all(v==0)|all(v==1))
  active <- setdiff(active,names(bounds)[bounds])
 }
 dd <- droplevels(g[ok & g$study %in% active,]);dd$y <- y[ok & g$study %in% active]
 p <- nlevels(dd$study)+if(component=='detection')4 else 3
 if(nrow(dd)-p<=1 || nlevels(dd$study)<2)return(ans)
 fit <- tryCatch(suppressWarnings(if(component=='detection')
   glm(y~0+study+x1+x2+x3+recovery,data=dd,family=quasibinomial(),control=glm.control(maxit=60,epsilon=1e-9))
   else lm(y~0+study+x1+x2+x3,data=dd)),error=function(e)NULL)
 if(is.null(fit)||fit$rank!=p)return(ans)
 if(component=='detection' && (!fit$converged || any(fitted(fit)<1e-8|fitted(fit)>1-1e-8)))return(ans)
 nm <- c('x1','x2','x3'); b <- coef(fit)[nm]
 conventional <- tryCatch({
  se <- sqrt(diag(vcov(fit)))[nm];df <- rep(df.residual(fit),3)
  cbind(b,se,df,2*pt(-abs(b/se),df))
 },error=function(e)matrix(NA_real_,3,4))
 robust <- tryCatch(suppressWarnings({
  v <- vcovCR(fit,cluster=dd$study,type='CR2')
  z <- coef_test(fit,vcov=v,test='Satterthwaite',coefs=nm)
  cbind(z$beta,z$SE,z$df_Satt,z$p_Satt)
 }),error=function(e)matrix(NA_real_,3,4))
 ans[1:3,] <- conventional;ans[4:6,] <- robust;ans
}
one <- function(scenario,signal,replicate,seed) {
 set.seed(seed)
 latent <- normal.genes(U)
 if(scenario=='study_dependence')latent <- sqrt(.4)*normal.genes(length(studies))[si[run.group],]+sqrt(.3)*normal.genes(N)[run.group,]+sqrt(.3)*latent
 if(scenario=='spatial_50km')latent <- Lrun%*%latent
 bd <- matrix(0,3,K);bp <- bd
 if(signal=='mixed'){bd[,is.signal]<-log(2);bp[,is.signal]<-.5}
 probs <- plogis(alpha+.75*g$recovery+X%*%bd)
 fraction <- W%*%(pnorm(latent)<probs[run.group,])
 err <- normal.genes(N)
 if(scenario=='study_dependence')for(j in 1:3)err <- err+X[,j]*normal.genes(length(studies))[si,]
 if(scenario=='dominant_heterogeneity')err <- err*ifelse(g$study_id==dominant,4,1)
 if(scenario=='spatial_50km')err <- Lgroup%*%err
 if(scenario=='heavy_tails')err <- qt(pmin(1-1e-12,pmax(1e-12,pnorm(err))),df=3)/sqrt(3)
 positive <- mu+X%*%bp+err;positive[fraction==0]<-NA_real_
 a <- matrix(NA_real_,J,4)
 for(k in 1:K) {
  ix <- ((k-1)*12+1):(k*12)
  a[ix[1:6],] <- fit.component(fraction[,k],'detection')
  a[ix[7:12],] <- fit.component(positive[,k],'positive')
 }
 truth <- if(signal=='mixed')ifelse(idx$gene %in% markers,ifelse(idx$component=='detection',log(2),.5),0) else rep(0,J)
 b <- a[,1];se <- a[,2];df <- a[,3];p <- a[,4]
 width <- 2*qt(.975,df)*se
 valid <- apply(is.finite(a),1,all)&se>1e-10&df>0&is.finite(width)
 valid[is.na(valid)] <- FALSE
 p[!valid]<-1
 coverage <- abs(b-truth)<=width/2
 reject <- p<.05
 # Sufficient statistics; invalid fits explicitly contribute zero and remain in denominators.
 values <- cbind(valid=as.numeric(valid),b=b,b2=b*b,se=se,se2=se*se,df=df,coverage=as.numeric(coverage),reject=as.numeric(reject),width=width,truth=truth)
 values[!valid,2:9] <- 0
 fdr <- lapply(c('conventional','CR2'),function(method){
  ii <- idx$method==method; r <- p.adjust(p[ii],method='BH')<=.05;tr <- truth[ii]!=0
  data.frame(scenario,signal,replicate,method,discoveries=sum(r),false=sum(r&!tr),true=sum(r&tr),
    n_nonnull=sum(tr),FDP=sum(r&!tr)/max(1,sum(r)),power=if(any(tr))sum(r&tr)/sum(tr) else NA_real_,valid=sum(valid[ii]))
 })
 list(sums=values,fdr=do.call(rbind,fdr))
}
scenarios <- c('independent','study_dependence','dominant_heterogeneity','spatial_50km','heavy_tails')
if(mode=='benchmark') {
 tm <- system.time(for(i in 1:5)invisible(one('independent','null',i,902000+i)))
 print(tm);quit(save='no')
}
if(mode=='validate') {
 set.seed(390318)
 stopifnot(max(abs(diag(Lgroup%*%t(Lgroup))-1))<1e-10,
           max(abs(diag(Lrun%*%t(Lrun))-1))<1e-10)
 y <- mu[,1]+rnorm(N)
 full <- lm(y~0+study+x1+x2+x3,data=g)
 drop <- fit.component(y,'positive')
 z <- coef_test(full,vcov=vcovCR(full,cluster=g$study,type='CR2'),test='Satterthwaite',coefs=c('x1','x2','x3'))
 stopifnot(max(abs(drop[4:6,1]-z$beta))<1e-8,max(abs(drop[4:6,2]-z$SE))<1e-8,max(abs(drop[4:6,3]-z$df_Satt))<1e-8)
 p <- plogis(alpha[,1]+.75*g$recovery)
 draws <- W%*%(matrix(runif(U*10000),U,10000)<p[run.group])
 stopifnot(max(abs(rowMeans(draws)-p))<.025)
 a <- one('independent','mixed',1,390319);b <- one('independent','mixed',1,390319)
 stopifnot(identical(a,b),all(a$fdr$n_nonnull==length(markers)*6),all(table(idx$method)==K*6),
           all(a$sums[,'valid'] %in% c(0,1)),all(is.finite(a$sums)))
 writeLines('PASS: positive-model singleton invariance; spatial unit marginals; fractional-response marginal means; fixed current-panel test family; prespecified available signal markers; seed reproducibility.',file.path(out,'final_implementation_validation.txt'))
 print('Final implementation checks passed');quit(save='no')
}
tasks <- expand.grid(scenario=seq_along(scenarios),signal=c('null','mixed'),chunk=1:10,stringsAsFactors=FALSE)
dir.create(file.path(out,'chunks'),showWarnings=FALSE)
run.chunk <- function(i) {
 task <- tasks[i,];n <- if(task$signal=='null')100 else 50
 reps <- ((task$chunk-1)*n+1):(task$chunk*n)
 dest <- file.path(out,'chunks',sprintf('chunk_%03d.rds',i))
 if(file.exists(dest))return(i)
 sumv <- NULL;fd <- vector('list',n)
 for(j in seq_along(reps)) {
  seed <- 2026110307L+task$scenario*10000L+ifelse(task$signal=='mixed',5000L,0L)+reps[j]
  v <- one(scenarios[task$scenario],task$signal,reps[j],seed)
  if(is.null(sumv))sumv<-v$sums else sumv<-sumv+v$sums
  fd[[j]]<-v$fdr
 }
 metrics <- cbind(scenario=scenarios[task$scenario],signal=task$signal,idx,replicates=n,as.data.frame(sumv))
 saveRDS(list(metrics=metrics,fdr=do.call(rbind,fd)),dest,compress='gzip')
 cat(sprintf('Completed chunk %d/100: %s %s\n',i,scenarios[task$scenario],task$signal))
 i
}
start <- Sys.time()
res <- parallel::mclapply(seq_len(nrow(tasks)),run.chunk,mc.cores=as.integer(Sys.getenv('GENE_ATLAS_WORKERS','2')),mc.preschedule=FALSE)
stopifnot(all(vapply(res,is.numeric,logical(1))))
parts <- lapply(seq_len(nrow(tasks)),function(i)readRDS(file.path(out,'chunks',sprintf('chunk_%03d.rds',i))))
write.table(do.call(rbind,lapply(parts,`[[`,'metrics')),file.path(out,'chunk_sufficient_statistics.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(do.call(rbind,lapply(parts,`[[`,'fdr')),file.path(out,'replicate_FDR.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
writeLines(c(capture.output(sessionInfo()),paste('Elapsed minutes:',as.numeric(difftime(Sys.time(),start,units='mins')))),file.path(out,'final_R_session.txt'))
cat('FINAL CALIBRATION COMPLETE\n')
