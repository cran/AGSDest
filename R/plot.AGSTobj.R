
plot.AGSTobj<-function(x,main=c("primary trial","secondary trial"),print.pdf=FALSE,...){

if(!is.null(x[[1]]$t) && !is.null(x[[2]]$T) && !is.null(x[[3]]$t)){
#par.old=par(no.readonly = TRUE)
#par(mfrow=c(2,1),oma = c(0,12,0,0))




if(mode(print.pdf)=="logical" & print.pdf==FALSE){nf <- layout(matrix(c(1,2,3,4),2,2,byrow=TRUE), c(2,3), c(3,3));layout.show(nf)}
if(mode(print.pdf)=="logical" & print.pdf==TRUE)pdf("AGST.pdf")
if(mode(print.pdf)=="character")pdf(print.pdf)

  plot(1,1, type="n",ann=FALSE, axes=FALSE, col="black")
  box(col="black")

  text(1,1,paste(length(x[[1]]$t)," stage primary trial,\n",switch(x[[1]]$SF,"OBF","Pocock","Power family","Hwang-Shih-DeCani","self specified","NA")," boundaries\n at level ",x[[1]]$al,
  if(!is.null(x[[1]]$Imax))paste(",\nImax=",round(x[[1]]$Imax,digits=2)),
  if(!is.null(x[[1]]$delta))if(x[[1]]$delta!=0)paste(",\ndelta=",round(x[[1]]$delta,digits=2)),
  if(!is.null(x[[2]]$T))paste("\n\nAdaptation at stage T=",x[[2]]$T), " \nwith z = ",round(x[[2]]$z,digits=3)),cex=0.8)


  plot(x[[1]]$t,x[[1]]$b,axes=FALSE,ylim=c(0,max(x[[1]]$b)+1.5),xlab="Cumulative Information Fraction",ylab="Wald Teststatistic",cex.lab=1.25,type="l",
  main=main[1],
  lwd=1.5,cex.main=0.9,cex.lab=0.9,cex=0.9,cex.axis=0.9,cex.sub=0.9,...)

  axis(2,cex.axis=0.9)
  axis(1, at=x[[1]]$t, labels=round(x[[1]]$t,digits=3),cex.axis=0.9)
  points(x[[1]]$t,x[[1]]$b,lt=1,pch=19)

  points(x[[1]]$t[x[[2]]$T],x[[2]]$z,pch=21)




#if(print.pdf){dev.off();cat("\nFile pT.pdf is created");pdf("sT.pdf")}
#else{print("Click to switch to next plot");par(ask=TRUE)}

#if(locator(1)){

  plot(1,1, type="n",ann=FALSE, axes=FALSE, col="black")
  box(col="black")

  text(1,1,paste(length(x[[3]]$t)," stage secondary trial,\n",switch(x[[3]]$SF,"OBF","Pocock","Power family","Hwang-Shih-DeCani","self specified","NA")," boundaries\n at level ",round(x[[3]]$al,digits=3),
  if(!is.null(x[[3]]$Imax))paste(",\nImax=",round(x[[3]]$Imax,digits=2)),
  if(!is.null(x[[3]]$delta))if(x[[3]]$delta!=0)paste(",\ndelta=",round(x[[3]]$delta,digits=2)),
  if(!is.null(x[[4]]$T) && x[[4]]$T!=0){paste("\n\nTrial stops at stage T=",x[[4]]$T,"\nwith z = ",round(x[[4]]$z,digits=3))},
  paste("\n"),
  #if(!is.null(x$pvalue.r) || !is.null(x$pvalue.so))paste("\n"),
  if(!is.null(x$pvalue.r))paste("\npvalue.r = ",round(x$pvalue.r,digits=3)),
  if(!is.null(x$pvalue.so))paste("\npvalue.so = ",round(x$pvalue.so,digits=3)),
  #if(!is.null(x$cb.r) || !is.null(x$cb.so))paste("\n"),
  if(!is.null(x$cb.r))paste("\ncb.r = ",round(x$cb.r,digits=3)),
  if(!is.null(x$cb.so))paste("\ncb.so = ",round(x$cb.so,digits=3)),  
  #if(!is.null(x$est.ml) || !is.null(x$est.mu))paste("\n"),
  if(!is.null(x$est.ml))paste("\nest.ml = ",round(x$est.ml,digits=3)),
  if(!is.null(x$est.mu))paste("\nest.mu = ",round(x$est.mu,digits=3))),cex=0.8)


  plot(x[[3]]$t,x[[3]]$b,axes=FALSE,ylim=c(0,max(x[[3]]$b)+1.5),xlab="Cumulative Information Fraction",ylab="Wald Teststatistic",cex.lab=1.25,type="l",
  main=main[2],
  lwd=1.5,cex.main=0.9,cex.lab=0.9,cex=0.9,cex.axis=0.9,cex.sub=0.9,...)

  axis(2,cex.axis=0.9)
  axis(1, at=x[[3]]$t, labels=round(x[[3]]$t,digits=3),cex.axis=0.9)
  points(x[[3]]$t,x[[3]]$b,lt=1,pch=19)

  if(!is.null(x[[4]]$T) && x[[4]]$T!=0)points(x[[3]]$t[x[[4]]$T],x[[4]]$z,pch=21)


#  }

if(mode(print.pdf)=="logical" & print.pdf==TRUE){dev.off();cat("File AGST.pdf is created\n")}
if(mode(print.pdf)=="character"){dev.off();cat("File",print.pdf,"is created\n")}


#par(par.old)
}
}