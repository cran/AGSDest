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
if(!is.null(x$est.cons))
cat(paste("conservative estimate: ",round(x$est.cons,digits=3),"\n\n"))


if(is.null(x)){
  cat(paste("No design adaptation performed.\n"))
  }
else{

if(!is.null(x$t)){

cat(paste(length(x$t)," stage group sequential design"))
cat(paste("\n",expression(alpha),": ",round(x$al,digits=3),"  SF: ",x$SF,"  phi: ",x$phi,if(!is.null(x$Imax))paste("  Imax: ",round(x$Imax,digits=2)),if(!is.null(x$delta))if(x$delta!=0)paste("  delta: ",round(x$delta,digits=2)),if(!is.null(x$cp))if(x$cp!=0)paste("  cp: ",round(x$cp,digits=2))))

cat("\n")
mat=matrix(c(round(x$b,digits=3),round(x$a,digits=3),round(x$t,digits=3)),ncol=length(x$t),byrow=TRUE)
colnames(mat)=rep("",length(x$t))
rownames(mat)=c("Upper bounds","Lower bounds","Information fraction")
print(mat)

cat("\n");
if(!is.null(x$alab)){
  mat_als=matrix(c(round(x$als,digits=3),round(x$alab[1:length(x$t)],digits=3)),ncol=length(x$t),byrow=TRUE)
  colnames(mat_als)=rep("",length(x$t))
  rownames(mat_als)=c("als","alab")
  print(mat_als)
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
cat(paste("\n",expression(alpha),": ",round(x[[1]]$al,digits=3),"  SF: ",x[[1]]$SF,"  phi: ",x[[1]]$phi,if(!is.null(x[[1]]$Imax))paste("  Imax: ",round(x[[1]]$Imax,digits=2)),if(!is.null(x[[1]]$delta))if(x[[1]]$delta!=0)paste("  delta: ",round(x[[1]]$delta,digits=2)),if(!is.null(x[[1]]$cp))if(x[[1]]$cp!=0)paste("  cp: ",round(x[[1]]$cp,digits=2))))

cat("\n")
mat=matrix(c(round(x[[1]]$b,digits=3),round(x[[1]]$a,digits=3),round(x[[1]]$t,digits=3)),ncol=length(x[[1]]$t),byrow=TRUE)
colnames(mat)=rep("",length(x[[1]]$t))
rownames(mat)=c("Upper bounds","Lower bounds","Information fraction")
print(mat)

if(!is.null(x[[1]]$alab)){
  mat_als=matrix(c(round(x[[1]]$als,digits=3),round(x[[1]]$alab[1:length(x[[1]]$t)],digits=3)),ncol=length(x[[1]]$t),byrow=TRUE)
  colnames(mat_als)=rep("",length(x[[1]]$t))
  rownames(mat_als)=c("als","alab")
  print(mat_als)
  cat("\n")
  }

if(!is.null(x[[2]]$T)){
  cat("\ngroup sequential design outcome:\n")
  cat(paste("\n\tT: ",x[[2]]$T,"  z: ",round(x[[2]]$z,digits=3),"\n\n"))
  }
}
}
}
}
