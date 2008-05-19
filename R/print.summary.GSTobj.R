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
if(!is.null(x$est.cons))
cat(paste("conservative estimate: ",round(x$est.cons,digits=3),"\n\n"))

if(!is.null(x$GSD)){

cat(paste(length(x$GSD$a)," stage group sequential design"))
cat(paste("\n",expression(alpha),": ",round(x$GSD$al,digits=3),"  SF: ",x$GSD$SF,"  phi: ",x$GSD$phi,if(!is.null(x$Imax))paste("  Imax: ",round(x$Imax,digits=2)),if(!is.null(x$delta))if(x$delta!=0)paste("  delta: ",round(x$delta,digits=2)),if(!is.null(x$cp))if(x$cp!=0)paste("  cp: ",round(x$cp,digits=2))))

cat("\n")
mat=matrix(c(round(x$GSD$b,digits=3),round(x$GSD$a,digits=3),round(x$GSD$t,digits=3)),ncol=length(x$GSD$t),byrow=TRUE)
colnames(mat)=rep("",length(x$GSD$t))
rownames(mat)=c("Upper bounds","Lower bounds","Information fraction")
print(mat)

cat("\n");
if(!is.null(x$GSD$alab)){
  mat_als=matrix(c(round(x$GSD$als,digits=3),round(x$GSD$alab,digits=3)),ncol=length(x$GSD$t),byrow=TRUE)
  colnames(mat_als)=rep("",length(x$GSD$t))
  rownames(mat_als)=c("als","alab")
  print(mat_als)
  cat("\n\n")
  }

cat("\n\ngroup sequential design outcome:\n")
cat(paste("\n\tT: ",x$GSDo$T,"  z: ",round(x$GSDo$z,digits=3),"\n\n"))
}
}
