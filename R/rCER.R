`rCER` <-
function(h=0,pT,iD,level=NULL){
          pT$K<-length(pT$t);
          if(ifelse(!is.null(level),pT$al!=level,FALSE)){
             pT$al <- level;
             pT$b  <- compBounds(t=1:pT$K/pT$K, t2 = pT$t*pT$Imax, iuse = pT$SF, asf = NULL, 
                      alpha = level, phi = ifelse(is.null(pT$phi),0,pT$phi), 
                      ztrun = 8)#$upper.bounds
            }
            CO(h=0,pT=pT,iD=list(T=iD$T,z=iD$z-h*sqrt(pT$t[iD$T]*pT$Imax)))
         }

