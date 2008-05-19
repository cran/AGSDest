summary.AGSTobj<-function(object,ctype="b",ptype="b",etype="b",overwrite=FALSE,...){


if(class(object)=="AGSTobj"){

if(ctype!="n"){
  if(ctype=="r"){
    if(is.null(object$cb.r) || object$cb.r==0 || (!is.null(object$cb.r) && overwrite)){
      object$cb.r=cb.r.ad(object$pT,object$iD,object$sT,object$sTo)
    }    
  }
  if(ctype=="so"){
    if(is.null(object$cb.so) || object$cb.so==0 || (!is.null(object$cb.so) && overwrite)){
      if(object$sTo$z<object$sT$b[object$sTo$T] & object$sTo$T<object$sT$K){
        cat("cb.so : z < b[T]; Stopping rule NOT met.\n")
      }
      else object$cb.so=cb.so.ad(object$pT,object$iD,object$sT,object$sTo)
      }
  }
  if(ctype=="b"){
    if(is.null(object$cb.so) || object$cb.so==0 || (!is.null(object$cb.so) && overwrite)){
      if(object$sTo$z<object$sT$b[object$sTo$T] & object$sTo$T<object$sT$K){
        cat("cb.so : z < b[T]; Stopping rule NOT met. Only the repeated confidence bound is calculated.\n")
        }
      else object$cb.so=cb.so.ad(object$pT,object$iD,object$sT,object$sTo)
      }
    if(is.null(object$cb.r) || object$cb.r==0 || (!is.null(object$cb.r) && overwrite)){
      object$cb.r=cb.r.ad(object$pT,object$iD,object$sT,object$sTo)
    }
  }
}
  
if(ptype!="n"){  
  if(ptype=="r"){
    if(is.null(object$pvalue.r) || object$pvalue.r==0  || (!is.null(object$pvalue.r) && overwrite)){
      object$pvalue.r=P.r.ad(h=0,object$pT,object$iD,object$sT,object$sTo)
    }
  }
  if(ptype=="so"){
    if(is.null(object$pvalue.so) || object$pvalue.so==0   || (!is.null(object$pvalue.so) && overwrite)){
      if(object$sTo$z<object$sT$b[object$sTo$T]){
        cat("pvalue.so : z < b[T]; Stopping rule NOT met.\n")
      }
      else object$pvalue.so=P.so.ad(h=0,object$pT,object$iD,object$sT,object$sTo)
    }
  }
  if(ptype=="b"){
    if(is.null(object$pvalue.so) || object$pvalue.so==0   || (!is.null(object$pvalue.so) && overwrite)){
      if(object$sTo$z<object$sT$b[object$sTo$T]){
        cat("pvalue.so : z < b[T]; Stopping rule NOT met. Only the repeated p-value is calculated.\n")
        }
      else object$pvalue.so=P.so.ad(h=0,object$pT,object$iD,object$sT,object$sTo)
    }
    if(is.null(object$pvalue.r) || object$pvalue.r==0  || (!is.null(object$pvalue.r) && overwrite)){
      object$pvalue.r=P.r.ad(h=0,object$pT,object$iD,object$sT,object$sTo)
    }
  }                                                                                                 
}

if(etype!="n"){  
  if(etype=="ml"){
    if(is.null(object$est.ml) || object$est.ml==0  || (!is.null(object$est.ml) && overwrite)){
      object$est.ml=(object$pT$t[object$iD$T]*object$pT$Imax*object$iD$z*sqrt(object$pT$t[object$iD$T]*object$pT$Imax)+object$sT$t[object$sTo$T]*object$sT$Imax*object$sTo$z*sqrt(object$sT$t[object$sTo$T]*object$sT$Imax))/(object$pT$t[object$iD$T]*object$pT$Imax+object$sT$t[object$sTo$T]*object$sT$Imax)
    }
  }
  if(etype=="mu"){
    if(is.null(object$est.mu) || object$est.mu==0   || (!is.null(object$est.mu) && overwrite)){
      if(object$sTo$z<object$sT$b[object$sTo$T] & object$sTo$T<object$sT$K){
        cat("est.mu : z < b[T]; Stopping rule NOT met.\n")
      }
      else object$est.mu=cb.so.ad(object$pT,object$iD,object$sT,object$sTo,level=0.5)
    }
  }
  if(etype=="b"){
    if(is.null(object$est.mu) || object$est.mu==0   || (!is.null(object$est.mu) && overwrite)){
      if(object$sTo$z<object$sT$b[object$sTo$T] & object$sTo$T<object$sT$K){
        cat("est.mu : z < b[T]; Stopping rule NOT met. Only the maximum likelihood estimate is calculated.\n")
        }
      else object$est.mu=cb.so.ad(object$pT,object$iD,object$sT,object$sTo,level=0.5)
    }
    if(is.null(object$est.ml) || object$est.ml==0  || (!is.null(object$est.ml) && overwrite)){
      object$est.ml=(object$pT$t[object$iD$T]*object$pT$Imax*object$iD$z*sqrt(object$pT$t[object$iD$T]*object$pT$Imax)+object$sT$t[object$sTo$T]*object$sT$Imax*object$sTo$z*sqrt(object$sT$t[object$sTo$T]*object$sT$Imax))/(object$pT$t[object$iD$T]*object$pT$Imax+object$sT$t[object$sTo$T]*object$sT$Imax)
    }
  }                                                                                                 
}

return(object)
}
}

