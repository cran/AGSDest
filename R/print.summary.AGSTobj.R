print.summary.AGSTobj<-function(x,...){

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

cat(paste("Primary trial: \n\n"))
cat(paste(length(x[[1]]$t)," stage group sequential design"))
cat(paste("\n",expression(alpha),": ",round(x[[1]]$al,digits=3),"\t SF: ",x[[1]]$SF,"\t phi: ",x[[1]]$phi,if(!is.null(x[[1]]$Imax))paste("\t Imax: ",round(x[[1]]$Imax,digits=2)),if(!is.null(x[[1]]$delta))if(x[[1]]$delta!=0)paste("\t delta: ",round(x[[1]]$delta,digits=2))))
cat(paste("\n\n\tBoundaries:","\t "))
for(i in 1:length(x[[1]]$a))cat(paste(round(x[[1]]$b[i],digits=3),"\t "))
cat(paste("\n           ","\t\t "))
for(i in 1:length(x[[1]]$a))cat(paste(round(x[[1]]$a[i],digits=3),"\t "))
cat(paste("\n\tInformation:","\t "))
for(i in 1:length(x[[1]]$a))cat(paste(round(x[[1]]$t[i],digits=3),"\t "))
if(!is.null(x[[1]]$alab)){
  cat(paste("\n\n\tals:","\t\t "))
  for(i in 1:length(x[[1]]$a))cat(paste(round(x[[1]]$als[i],digits=3),"\t "))
  cat(paste("\n\talab:","\t\t "))
  for(i in 1:length(x[[1]]$a))cat(paste(round(x[[1]]$alab[i],digits=3),"\t "))
  }
cat("\n\n")
cat(paste("interim data: \n"))
cat(paste("\n\tT: ",x[[2]]$T,"\t z: ",round(x[[2]]$z,digits=3)))

cat(paste("\n\nSecondary trial: \n\n"))
cat(paste(length(x[[3]]$t)," stage group sequential design"))
cat(paste("\n crp: ",round(x[[3]]$al,digits=3),"\t SF: ",x[[3]]$SF,"\t phi: ",x[[3]]$phi,if(!is.null(x[[3]]$Imax))paste("\t Imax: ",round(x[[3]]$Imax,digits=2)),if(!is.null(x[[3]]$delta))if(x[[3]]$delta!=0)paste("\t delta: ",round(x[[3]]$delta,digits=2))))
cat(paste("\n\n\tBoundaries:","\t "))
for(i in 1:length(x[[3]]$a))cat(paste(round(x[[3]]$b[i],digits=3),"\t "))
cat(paste("\n           ","\t\t "))
for(i in 1:length(x[[3]]$a))cat(paste(round(x[[3]]$a[i],digits=3),"\t "))
cat(paste("\n\tInformation:","\t "))
for(i in 1:length(x[[3]]$a))cat(paste(round(x[[3]]$t[i],digits=3),"\t "))

cat(paste("\n\nSecondary trial outcome: \n"))
cat(paste("\n\tT: ",x[[4]]$T,"\t z: ",round(x[[4]]$z,digits=3),"\n\n"))

}

