print.GSTobj<-function(x,...){


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


if(is.null(x)){
  cat(paste("No design adaptation performed.\n"))
  }
else{

if(!is.null(x$t)){

cat(paste(length(x$t)," stage group sequential design"))
cat(paste("\n",expression(alpha),": ",round(x$al,digits=3),"\t SF: ",x$SF,"\t phi: ",x$phi,if(!is.null(x$Imax))paste("\t Imax: ",round(x$Imax,digits=2)),if(!is.null(x$delta))if(x$delta!=0)paste("\t delta: ",round(x$delta,digits=2))))
cat(paste("\n\nBoundaries:","\t "))
for(i in 1:length(x$t))cat(paste(round(x$b[i],digits=3),"\t "))
cat(paste("\n           ","\t "))
for(i in 1:length(x$t))cat(paste(round(x$a[i],digits=3),"\t "))
cat(paste("\nInformation:","\t "))
for(i in 1:length(x$t))cat(paste(round(x$t[i],digits=3),"\t "))
cat("\n");
if(!is.null(x$alab)){
  cat(paste("\n\nals:","\t\t "))
  for(i in 1:length(x$t))cat(paste(round(x$als[i],digits=3),"\t "))
  cat(paste("\nalab:","\t\t "))
  for(i in 1:length(x$t))cat(paste(round(x$alab[i],digits=3),"\t "))
  cat("\n\n")
  }
}
}


if(is.null(x[[1]])){
  cat(paste("No design adaptation performed.\n"))
  }
else{

if(!is.numeric(x[[1]])){
if(!is.null(x[[1]]$t)){

cat(paste(length(x[[1]]$t)," stage group sequential design"))
cat(paste("\n",expression(alpha),": ",round(x[[1]]$al,digits=3),"\t SF: ",x[[1]]$SF,"\t phi: ",x[[1]]$phi,if(!is.null(x[[1]]$Imax))paste("\t Imax: ",round(x[[1]]$Imax,digits=2)),if(!is.null(x[[1]]$delta))if(x[[1]]$delta!=0)paste("\t delta: ",round(x[[1]]$delta,digits=2))))
cat(paste("\n\nBoundaries:","\t "))
for(i in 1:length(x[[1]]$t))cat(paste(round(x[[1]]$b[i],digits=3),"\t "))
cat(paste("\n           ","\t "))
for(i in 1:length(x[[1]]$t))cat(paste(round(x[[1]]$a[i],digits=3),"\t "))
cat(paste("\nInformation:","\t "))
for(i in 1:length(x[[1]]$t))cat(paste(round(x[[1]]$t[i],digits=3),"\t "))
if(!is.null(x[[1]]$alab)){
  cat(paste("\n\nals:","\t\t "))
  for(i in 1:length(x[[1]]$t))cat(paste(round(x[[1]]$als[i],digits=3),"\t "))
  cat(paste("\nalab:","\t\t "))
  for(i in 1:length(x[[1]]$t))cat(paste(round(x[[1]]$alab[i],digits=3),"\t "))
  cat("\n\n")
  }
if(!is.null(x[[2]]$T)){
  cat("\n\ngroup sequential design outcome:\n")
  cat(paste("\n\tT: ",x[[2]]$T,"\t z: ",round(x[[2]]$z,digits=3),"\n\n"))
  }
}
}
}
}
