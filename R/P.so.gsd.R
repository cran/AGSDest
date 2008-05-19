`P.so.gsd` <-
function(h=0,GSD,GSDo){
        if((GSDo$T<length(GSD$t))&(GSDo$z<GSD$b[GSDo$T])){ stop("stopping rule NOT met")
        } else {
          if(GSDo$T==1){ res=1-pnorm(GSDo$z-h*sqrt(GSD$t[1]*GSD$Imax))
          } else {
            res=seqmon(a=GSD$a[1:GSDo$T],
                   b=c(GSD$b[1:(GSDo$T-1)],GSDo$z)-h*sqrt(GSD$t[1:GSDo$T]*GSD$Imax),
                   t=GSD$t[1:GSDo$T]*GSD$Imax,int=500*array(c(1),GSDo$T)
                   )[2*GSDo$T]
          }
        }
      return(res)      
      }

