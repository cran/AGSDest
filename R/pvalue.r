`pvalue` <-
function(object, type="b"){

if(class(object)=="GSTobj"){
  if(!is.null(object$GSDo$z)|!is.null(object$GSDo$T)){
  if(type=="r"){
    return(list(pvalue.r=P.r.gsd(h=0,object$GSD,object$GSDo)))
  }
  if(type=="so"){
    if(object$GSDo$z<object$GSD$b[object$GSDo$T] & object$GSDo$T<object$GSD$K){
      cat("pvalue.so : z < b[T]; Stopping rule NOT met.\n")
      }
    else{ return(list(pvalue.so=P.so.gsd(h=0,object$GSD,object$GSDo)))
          }
  }
  if(type=="b"){
    if(object$GSDo$z<object$GSD$b[object$GSDo$T] & object$GSDo$T<object$GSD$K){
      cat("pvalue.so : z < b[T]; Stopping rule NOT met. Only the repeated p-value is calculated.\n")
      return(list(pvalue.r=P.r.gsd(h=0,object$GSD,object$GSDo)))
      }
    else {return(list(pvalue.so=P.so.gsd(h=0,object$GSD,object$GSDo),pvalue.r=P.r.gsd(h=0,object$GSD,object$GSDo)))
          }
  }
  }
  else{
    print("interim data missing")
    }
}
if(class(object)=="AGSTobj"){
  if(type=="r"){
    return(list(pvalue.r=P.r.ad(h=0,object$pT,object$iD,object$sT,object$sTo)))
  }
  if(type=="so"){
    if(object$sTo$z<object$sT$b[object$sTo$T] & object$sTo$T<object$sT$K){
      cat("pvalue.so : z < b[T]; Stopping rule NOT met.\n")
      }
    else {return(list(pvalue.so=P.so.ad(h=0,object$pT,object$iD,object$sT,object$sTo)))
          }
  }
  if(type=="b"){
    if(object$sTo$z<object$sT$b[object$sTo$T] & object$sTo$T<object$sT$K){
      cat("pvalue.so : z < b[T]; Stopping rule NOT met. Only the repeated p-value is calculated.\n")
      return(list(pvalue.r=P.r.ad(h=0,object$pT,object$iD,object$sT,object$sTo)))
      }
    else {return(list(pvalue.so=P.so.ad(h=0,object$pT,object$iD,object$sT,object$sTo),
                      pvalue.r=P.r.ad(h=0,object$pT,object$iD,object$sT,object$sTo)))
          }
  }
}
}
