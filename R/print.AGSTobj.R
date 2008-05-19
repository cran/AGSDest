print.AGSTobj<-function(x,...){


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
if(!is.null(x$est.cons))
cat(paste("conservative estimate: ",round(x$est.cons,digits=3),"\n\n"))

cat(paste("Primary trial: \n\n"))
cat(paste(length(x[[1]]$t)," stage group sequential design"))
cat(paste("\n",expression(alpha),": ",round(x[[1]]$al,digits=3),"  SF: ",x[[1]]$SF,"  phi: ",x[[1]]$phi,if(!is.null(x[[1]]$Imax))paste("  Imax: ",round(x[[1]]$Imax,digits=2)),if(!is.null(x[[1]]$delta))if(x[[1]]$delta!=0)paste("  delta: ",round(x[[1]]$delta,digits=2)),if(!is.null(x[[1]]$cp))if(x[[1]]$cp!=0)paste("  cp: ",round(x[[1]]$cp,digits=2))))

cat("\n")
mat=matrix(c(round(x[[1]]$b,digits=3),round(x[[1]]$a,digits=3),round(x[[1]]$t,digits=3)),ncol=length(x[[1]]$t),byrow=TRUE)
colnames(mat)=rep("",length(x[[1]]$t))
rownames(mat)=c("Upper bounds","Lower bounds","Information fraction")
print(mat)

cat("\n");
if(!is.null(x[[1]]$alab)){
  mat_als=matrix(c(round(x[[1]]$als,digits=3),round(x[[1]]$alab[1:length(x[[1]]$t)],digits=3)),ncol=length(x[[1]]$t),byrow=TRUE)
  colnames(mat_als)=rep("",length(x[[1]]$t))
  rownames(mat_als)=c("als","alab")
  print(mat_als)
  cat("\n")
  }

cat("\n\n")
cat(paste("interim data: \n"))
cat(paste("\n\tT: ",x[[2]]$T,"  z: ",round(x[[2]]$z,digits=3)))

cat(paste("\n\nSecondary trial: \n\n"))
cat(paste(length(x[[3]]$t)," stage group sequential design"))
cat(paste("\n cer: ",round(x[[3]]$al,digits=3),"  SF: ",x[[3]]$SF,"  phi: ",x[[3]]$phi,if(!is.null(x[[3]]$Imax))paste("  Imax: ",round(x[[3]]$Imax,digits=2)),if(!is.null(x[[3]]$delta))if(x[[3]]$delta!=0)paste("  delta: ",round(x[[3]]$delta,digits=2)),if(!is.null(x[[3]]$cp))if(x[[3]]$cp!=0)paste("  cp: ",round(x[[3]]$cp,digits=2))))


cat("\n")
mat=matrix(c(round(x[[3]]$b,digits=3),round(x[[3]]$a,digits=3),round(x[[3]]$t,digits=3)),ncol=length(x[[3]]$t),byrow=TRUE)
colnames(mat)=rep("",length(x[[3]]$t))
rownames(mat)=c("Upper bounds","Lower bounds","Information fraction")
print(mat)

if(!is.null(x[[4]]$T) && x[[4]]$T!=0) {
  cat(paste("\n\nSecondary trial outcome: \n"))
  cat(paste("\n\tT: ",x[[4]]$T,"  z: ",round(x[[4]]$z,digits=3),"\n\n"))
  }
else cat(paste("\n"))
}

