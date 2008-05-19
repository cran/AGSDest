`comp.z2` <-
function(pT,pTo,iD){
               K<-length(pT$t);
               if(iD$T>=pTo$T) { print("Error: Cannot compute z2 for pTo$T<=iD$T") ; rep(NaN,K)
               } else {
                 (pTo$z*sqrt(pT$t[pTo$T]*pT$Imax)-ifelse(iD$T==0,0,iD$z*sqrt(pT$t[iD$T]*pT$Imax)))/
                 sqrt(pT$t[pTo$T]*pT$Imax-ifelse(iD$T==0,0,pT$t[iD$T]*pT$Imax))
                }
          }

