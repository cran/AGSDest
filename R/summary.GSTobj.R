summary.GSTobj<-function(object,ctype="b",ptype="b",etype="a",overwrite=FALSE,...){


if(class(object)=="GSTobj"){

if(!is.null(object$GSDo$z) & (!is.null(object$GSDo$T))){
if(ctype!="n"){  
  if(ctype=="r"){
    if(is.null(object$cb.r) || object$cb.r==0 || (!is.null(object$cb.r) && overwrite)){
      object$cb.r=cb.r.gsd(object$GSD,object$GSDo)
    }
  }
  if(ctype=="so"){
    if(is.null(object$cb.so) || object$cb.so==0 ||  (!is.null(object$cb.so) && overwrite)){
      if(object$GSDo$z<object$GSD$b[object$GSDo$T] & object$GSDo$T<object$GSD$K){
        cat("cb.so : z < b[T]; Stopping rule NOT met.\n")
        }
      else object$cb.so=cb.so.gsd(object$GSD,object$GSDo)
    }
  }
  if(ctype=="b"){
    if(is.null(object$cb.so) || object$cb.so==0  || (!is.null(object$cb.so) && overwrite)){
      if(object$GSDo$z<object$GSD$b[object$GSDo$T] & object$GSDo$T<object$GSD$K){
        cat("cb.so : z < b[T]; Stopping rule NOT met. Only the repeated confidence bound is calculated.\n")
        }
      else object$cb.so=cb.so.gsd(object$GSD,object$GSDo)
      }
    if(is.null(object$cb.r) || object$cb.r==0  || (!is.null(object$cb.r) && overwrite))object$cb.r=cb.r.gsd(object$GSD,object$GSDo)
  }
}
  
if(ptype!="n"){  
  if(ptype=="r"){
    if(is.null(object$pvalue.r) || object$pvalue.r==0  || (!is.null(object$pvalue.r) && overwrite))
    object$pvalue.r=P.r.gsd(h=0,object$GSD,object$GSDo)
  }
  if(ptype=="so"){
    if(is.null(object$pvalue.so) || object$pvalue.so==0   || (!is.null(object$pvalue.so) && overwrite)){
      if(object$GSDo$z<object$GSD$b[object$GSDo$T] & object$GSDo$T<object$GSD$K){
        cat("pvalue.so : z < b[T]; Stopping rule NOT met.\n")
        }
      else object$pvalue.so=P.so.gsd(h=0,object$GSD,object$GSDo)
    }
  }
  if(ptype=="b"){
    if(is.null(object$pvalue.so) || object$pvalue.so==0   || (!is.null(object$pvalue.so) && overwrite)){
      if(object$GSDo$z<object$GSD$b[object$GSDo$T] & object$GSDo$T<object$GSD$K){
        cat("pvalue.so : z < b[T]; Stopping rule NOT met. Only the repeated p-value is calculated.\n")
        }
      else object$pvalue.so=P.so.gsd(h=0,object$GSD,object$GSDo)      
      }
    if(is.null(object$pvalue.r) || object$pvalue.r==0   || (!is.null(object$pvalue.r) && overwrite))object$pvalue.r=P.r.gsd(h=0,object$GSD,object$GSDo)
  }
}


if(etype!="n"){  
  if(etype=="ml"){
    if(is.null(object$est.ml) || object$est.ml==0  || (!is.null(object$est.ml) && overwrite))
    object$est.ml=object$GSDo$z/sqrt(object$GSD$t[object$GSDo$T]*object$GSD$Imax)
  }
  if(etype=="mu"){
    if(is.null(object$est.mu) || object$est.mu==0   || (!is.null(object$est.mu) && overwrite)){
      if(object$GSDo$z<object$GSD$b[object$GSDo$T]){
        cat("est.mu : z < b[T]; Stopping rule NOT met.\n")
        }
      else object$est.mu=cb.so.gsd(object$GSD,object$GSDo,level=0.5)
    }
  }
  if(etype=="cons"){
    if(is.null(object$est.cons) || object$est.cons==0  || (!is.null(object$est.cons) && overwrite))
    object$est.cons=cb.r.gsd(object$GSD,object$GSDo,level=0.5)
  }

  if(etype=="a"){
    if(is.null(object$est.mu) || object$est.mu==0   || (!is.null(object$est.mu) && overwrite)){
      if(object$GSDo$z<object$GSD$b[object$GSDo$T]){
        cat("est.mu : z < b[T]; Stopping rule NOT met. Only the maximum likelihood estimate is calculated.\n")
        }
      else object$est.mu=cb.so.gsd(object$GSD,object$GSDo,level=0.5)      
    }
    if(is.null(object$est.ml) || object$est.ml==0   || (!is.null(object$est.ml) && overwrite))object$est.ml=object$GSDo$z/sqrt(object$GSD$t[object$GSDo$T]*object$GSD$Imax)
    if(is.null(object$est.cons) || object$est.cons==0  || (!is.null(object$est.cons) && overwrite)) object$est.cons=cb.r.gsd(object$GSD,object$GSDo,level=0.5)
  }
}


return(object)
}
else{print("Missing interim data")}
}
}
