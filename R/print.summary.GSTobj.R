print.summary.GSTobj<-function(x,...){


if(!is.null(x$type))cat(paste("type: ",x$type,"\n\n"))

if(!is.null(x$cb.r))
cat(paste("repeated lower confidence bound: ",round(x$cb.r,digits=3),"\n\n"))
if(!is.null(x$cb.so))
cat(paste("stage-wise adjusted lower confidence bound: ",round(x$cb.so,digits=3),"\n\n"))
if(!is.null(x$pvalue.r))
cat(paste("repeated p-value: ",round(x$pvalue.r,digits=3),"\n\n"))
if(!is.null(x$pvalue.so))
cat(paste("stage-wise adjusted p-value: ",round(x$pvalue.so,digits=3),"\n\n"))

if(!is.null(x$est.ml))
cat(paste("maximum likelihood estimate: ",round(x$est.ml,digits=3),"\n\n"))
if(!is.null(x$est.mu))
cat(paste("median unbiased estimate: ",round(x$est.mu,digits=3),"\n\n"))

if(!is.null(x$GSD)){

cat(paste(length(x$GSD$a)," stage group sequential design"))
cat(paste("\n",expression(alpha),": ",round(x$GSD$al,digits=3),"\t SF: ",x$GSD$SF,"\t phi: ",x$GSD$phi,if(!is.null(x$Imax))paste("\t Imax: ",round(x$Imax,digits=2)),if(!is.null(x$delta))if(x$delta!=0)paste("\t delta: ",round(x$delta,digits=2))))
cat(paste("\n\nBoundaries:","\t "))
for(i in 1:length(x$GSD$a))cat(paste(round(x$GSD$b[i],digits=3),"\t "))
cat(paste("\n           ","\t "))
for(i in 1:length(x$GSD$a))cat(paste(round(x$GSD$a[i],digits=3),"\t "))
cat(paste("\nInformation:","\t "))
for(i in 1:length(x$GSD$a))cat(paste(round(x$GSD$t[i],digits=3),"\t "))
if(!is.null(x$GSD$alab)){
  cat(paste("\n\nals:","\t\t "))
  for(i in 1:length(x$GSD$a))cat(paste(round(x$GSD$als[i],digits=3),"\t "))
  cat(paste("\nalab:","\t\t "))
  for(i in 1:length(x$GSD$a))cat(paste(round(x$GSD$alab[i],digits=3),"\t "))
  cat("\n\n")
  }
cat("\n\ngroup sequential design outcome:\n")
cat(paste("\n\tT: ",x$GSDo$T,"\t z: ",round(x$GSDo$z,digits=3),"\n\n"))
}
}
