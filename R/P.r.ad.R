`P.r.ad` <-
function(h=0,pT,iD,sT,sTo,prec=0.001){
      K<-length(pT$a)
      P2=P.r.gsd(h=0,GSD=sT,GSDo=sTo)
      res=bisearch(f=function(u) {cer.r(u,K,pT,iD)-P2},c(0,1),signfl=-1, signfu=1,tol=prec)$root
      
      return(res)
      }

