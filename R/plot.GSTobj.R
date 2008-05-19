
plot.GSTobj<-function(x,main="GSD",print.pdf=FALSE,...){

#par.old=par(no.readonly = TRUE)
#par(mfrow=c(1,1),oma = c(0,12,0,0))

if(is.null(x)){
     cat(paste("No design adaptation performed.\n"))
     return(0)
     }
else{
 if(!is.null(x$t)){
    
  if(mode(print.pdf)=="logical" & print.pdf==FALSE){nf <- layout(matrix(c(1,2),1,2,byrow=TRUE), c(2.2,3), c(3,3));layout.show(nf)}
  if(mode(print.pdf)=="logical" & print.pdf==TRUE)pdf("GST.pdf")
  if(mode(print.pdf)=="character")pdf(print.pdf)
 
  plot(1,1, type="n",ann=FALSE, axes=FALSE, col="black")
  box(col="black")

  text(1,1,paste(length(x$t)," stage GSD,\n",switch(x$SF,"OBF","Pocock","Power family","Hwang-Shih-DeCani","self specified","NA")," boundaries \nat level ",round(x$al,digits=3),
  if(!is.null(x$Imax))paste(",\nImax=",round(x$Imax,digits=2)),
  if(!is.null(x$delta))if(x$delta!=0)paste(",\ndelta=",round(x$delta,digits=2)),
  paste("\n"),
  if(!is.null(x$pvalue.r))paste("\npvalue.r = ",round(x$pvalue.r,digits=3)),
  if(!is.null(x$pvalue.so))paste("\npvalue.so = ",round(x$pvalue.so,digits=3)),
  if(!is.null(x$cb.r))paste("\ncb.r = ",round(x$cb.r,digits=3)),
  if(!is.null(x$cb.so))paste("\ncb.so = ",round(x$cb.so,digits=3)),
  if(!is.null(x$est.ml))paste("\nest.ml = ",round(x$est.ml,digits=3)),
  if(!is.null(x$est.mu))paste("\nest.mu = ",round(x$est.mu,digits=3))),cex=0.8)

  
  plot(x$t,x$b,axes=FALSE,ylim=c(0,max(x$b)+1.5),xlab="Cumulative Information Fraction",ylab="Wald Teststatistic",cex.lab=1.25,type="l",
  main=main,
  lwd=1.5,cex.main=0.9,cex.lab=0.9,cex=0.9,cex.axis=0.9,cex.sub=0.9,...)
  
  axis(2,cex.axis=0.9)
  axis(1, at=x$t, labels=(round(x$t,digits=3)),cex.axis=0.9)
  points(x$t,x$b,lt=1,pch=19)

  #legend(x$t[1],min(x$b),c("Boundaries"),lty=c(1),pch=c(19),cex=0.8,...)


  if(mode(print.pdf)=="logical"& print.pdf==TRUE){dev.off();cat("File GST.pdf is created\n")}
  if(mode(print.pdf)=="character"){dev.off();cat("File",print.pdf,"is created\n")}
  }
}

if(is.null(x[[1]])){
    cat(paste("No design adaptation performed.\n"))
    }
else{
  if(!is.numeric(x[[1]])){
  if(!is.null(x[[1]]$t)){
   
  if(mode(print.pdf)=="logical" & print.pdf==FALSE){nf <- layout(matrix(c(1,2),1,2,byrow=TRUE), c(2.2,3), c(3,3));layout.show(nf)}
  if(mode(print.pdf)=="logical"& print.pdf==TRUE)pdf("GST.pdf")
  if(mode(print.pdf)=="character")pdf(print.pdf)
  
  plot(1,1, type="n",ann=FALSE, axes=FALSE, col="black")
  box(col="black")

  text(1,1,paste(length(x[[1]]$t)," stage GSD,\n",switch(x[[1]]$SF,"OBF","Pocock","Power family","Hwang-Shih-DeCani","self specified","NA")," boundaries \nat level ",round(x[[1]]$al,digits=3),
  if(!is.null(x[[1]]$Imax))paste(",\nImax=",round(x[[1]]$Imax,digits=2)),
  if(!is.null(x[[1]]$delta))if(x[[1]]$delta!=0)paste(",\ndelta=",round(x[[1]]$delta,digits=2)),
  if(!is.null(x[[2]]$T))paste("\n\nTrial stops at look T = ",x[[2]]$T,"\nwith z = ",round(x[[2]]$z,digits=3)),
  paste("\n"),
  if(!is.null(x$pvalue.r))paste("\npvalue.r = ",round(x$pvalue.r,digits=3)),
  if(!is.null(x$pvalue.so))paste("\npvalue.so = ",round(x$pvalue.so,digits=3)),
  if(!is.null(x$cb.r))paste("\ncb.r = ",round(x$cb.r,digits=3)),
  if(!is.null(x$cb.so))paste("\ncb.so = ",round(x$cb.so,digits=3)),
  if(!is.null(x$est.ml))paste("\nest.ml = ",round(x$est.ml,digits=3)),
  if(!is.null(x$est.mu))paste("\nest.mu = ",round(x$est.mu,digits=3))),cex=0.8)

  
  plot(x[[1]]$t,x[[1]]$b,axes=FALSE,ylim=c(0,max(x[[1]]$b)+1.5),xlab="Cumulative Information Fraction",ylab="Wald Teststatistic",cex.lab=1.25,type="l",
  main=main,
  lwd=1.5,cex.main=0.9,cex.lab=0.9,cex=0.9,cex.axis=0.9,cex.sub=0.9,...)
  
  axis(2,cex.axis=0.9)
  axis(1, at=x[[1]]$t, labels=(round(x[[1]]$t,digits=3)),cex.axis=0.9)
  points(x[[1]]$t,x[[1]]$b,lt=1,pch=19)

  if(!is.null(x[[2]]$T)||x[[2]]$T==0){
    points(x[[1]]$t[x[[2]]$T],x[[2]]$z,pch=21)
    #legend(x[[1]]$t[1],min(x[[1]]$b,x[[2]]$z)-0.3,c("Standardized test statistic","Boundaries"),lty=c(0,1),pch=c(21,19),cex=0.8,...)
    }
  #else legend(x[[1]]$t[1],min(x[[1]]$b),c("Boundaries"),lty=c(1),pch=c(19),cex=0.8,...)


  if(mode(print.pdf)=="logical" & print.pdf==TRUE){dev.off();cat("File GST.pdf is created\n")}
  if(mode(print.pdf)=="character"){dev.off();cat("File",print.pdf,"is created\n")}
  }
  }
}
}