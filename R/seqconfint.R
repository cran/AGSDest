`seqconfint` <-
function(object, type="b",level=NULL){


if(class(object)=="GSTobj"){
  if(is.null(level))level=object$GSD$al

  if(type=="r"){
    if(level==0.5)return(list(est.cons=cb.r.gsd(object$GSD,object$GSDo,level)))
    else return(list(cb.r=cb.r.gsd(object$GSD,object$GSDo,level)))
  }
  if(type=="so"){
    if(object$GSDo$z<object$GSD$b[object$GSDo$T] & object$GSDo$T<object$GSD$K){
      cat("cb.so : z < b[T]; Stopping rule NOT met.\n")
      }
    else{ 
	if(level==0.5)return(list(est.mu=cb.so.gsd(object$GSD,object$GSDo,level)))
	else return(list(cb.so=cb.so.gsd(object$GSD,object$GSDo,level)))
    }
  }
  if(type=="b"){
    if(object$GSDo$z<object$GSD$b[object$GSDo$T] & object$GSDo$T<object$GSD$K){
      cat("cb.so : z < b[T]; Stopping rule NOT met. Only the repeated confidence bound is calculated.\n")
      if(level==0.5)return(list(est.cons=cb.r.gsd(object$GSD,object$GSDo,level)))
      else return(list(cb.r=cb.r.gsd(object$GSD,object$GSDo,level)))
      }
    else{
	if(level==0.5)return(list(est.mu=cb.so.gsd(object$GSD,object$GSDo,level),
                  est.cons=cb.r.gsd(object$GSD,object$GSDo,level)))
	else return(list(cb.so=cb.so.gsd(object$GSD,object$GSDo,level),
                  cb.r=cb.r.gsd(object$GSD,object$GSDo,level)))
    }
  }
}
if(class(object)=="AGSTobj"){
  if(is.null(level))level=object$pT$al

  if(type=="r"){
    if(level==0.5)return(list(est.cons=cb.r.ad(object$pT,object$iD,object$sT,object$sTo,level)))
    else return(list(cb.r=cb.r.ad(object$pT,object$iD,object$sT,object$sTo,level)))
  }
  if(type=="so"){
    if(object$sTo$z<object$sT$b[object$sTo$T] & object$sTo$T<object$sT$K){
      cat("cb.so : z < b[T]; Stopping rule NOT met.\n")
      }
    else{ 
	if(level==0.5)return(list(est.mu=cb.so.ad(object$pT,object$iD,object$sT,object$sTo,level)))
	else return(list(cb.so=cb.so.ad(object$pT,object$iD,object$sT,object$sTo,level)))
    }
  }
  if(type=="b"){
    if(object$sTo$z<object$sT$b[object$sTo$T] & object$sTo$T<object$sT$K){
      cat("cb.so : z < b[T]; Stopping rule NOT met. Only the repeated confidence bound is calculated.\n")
      if(level==0.5)return(list(est.cons=cb.r.ad(object$pT,object$iD,object$sT,object$sTo,level)))
      else return(list(cb.r=cb.r.ad(object$pT,object$iD,object$sT,object$sTo,level)))
      }
    else{
	if(level==0.5)return(list(est.mu=cb.so.ad(object$pT,object$iD,object$sT,object$sTo,level),est.cons=cb.r.ad(object$pT,object$iD,object$sT,object$sTo,level)))
	else return(list(cb.so=cb.so.ad(object$pT,object$iD,object$sT,object$sTo,level),cb.r=cb.r.ad(object$pT,object$iD,object$sT,object$sTo,level)))
    }
  }
}
}

