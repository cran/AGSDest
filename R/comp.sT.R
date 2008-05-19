`comp.sT` <-
function(pT,iD){
             sT=list(K=K<-length(pT$t)-iD$T,
                  ce0=cerr(h=0,pT=pT,iD=iD),
                  a=pT$a[1:K],
                  b=pbounds(h=0,pT,iD),
                  t=pT$t[(iD$T+1):length(pT$t)]*pT$Imax-pT$t[iD$T]*pT$Imax)
             class(sT)<-"GSD"
             return(sT)
           }

