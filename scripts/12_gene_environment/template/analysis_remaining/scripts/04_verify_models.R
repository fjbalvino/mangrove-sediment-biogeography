args<-commandArgs(TRUE);o<-file.path(args[1],'calibration_final')
library(clubSandwich)
g<-read.delim(file.path(o,'design.tsv'));g$study<-factor(g$study_id)
cases<-read.delim(file.path(o,'verification_cases.tsv'));ans<-list()
for(i in unique(cases$case)) {
 y<-cases$y[cases$case==i];component<-cases$component[cases$case==i][1]
 ok<-is.finite(y);counts<-table(g$study[ok]);active<-names(counts)[counts>=2]
 if(component=='detection'){bounds<-tapply(y,g$study,function(v)all(v==0)|all(v==1));active<-setdiff(active,names(bounds)[bounds])}
 dd<-droplevels(g[ok&g$study%in%active,]);dd$y<-y[ok&g$study%in%active]
 p<-nlevels(dd$study)+if(component=='detection')4 else 3
 z<-matrix(NA_real_,6,4)
 if(nrow(dd)-p>1&&nlevels(dd$study)>=2) {
  fit<-tryCatch(suppressWarnings(if(component=='detection')glm(y~0+study+x1+x2+x3+recovery,dd,family=quasibinomial(),control=glm.control(epsilon=1e-11,maxit=60)) else lm(y~0+study+x1+x2+x3,dd)),error=function(e)NULL)
  if(!is.null(fit)&&fit$rank==p&&(component=='positive'||(fit$converged&&all(fitted(fit)>1e-8&fitted(fit)<1-1e-8)))) {
   nm<-c('x1','x2','x3');bb<-coef(fit)[nm];se<-sqrt(diag(vcov(fit)))[nm];df<-rep(df.residual(fit),3)
   z[1:3,]<-cbind(bb,se,df,2*pt(-abs(bb/se),df))
   zz<-tryCatch(suppressWarnings(coef_test(fit,vcov=vcovCR(fit,cluster=dd$study,type='CR2'),test='Satterthwaite',coefs=nm)),error=function(e)NULL)
   if(!is.null(zz))z[4:6,]<-cbind(zz$beta,zz$SE,zz$df_Satt,zz$p_Satt)
  }
 }
 ans[[length(ans)+1]]<-data.frame(case=i,component,method=rep(c('conventional','CR2'),each=3),predictor=rep(1:3,2),beta=z[,1],se=z[,2],df=z[,3],p=z[,4])
}
write.table(do.call(rbind,ans),file.path(o,'verification_R.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
