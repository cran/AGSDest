`comp.z1` <-
function(pT,iD,sT,sTo){
                 (iD$z*sqrt(pT$t[iD$T]*pT$Imax)+sTo$z*sqrt(sT$t[sTo$T]*sT$Imax))/sqrt(pT$t[iD$T]*pT$Imax+sT$t[sTo$T]*sT$Imax)
          }

