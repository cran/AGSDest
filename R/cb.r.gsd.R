`cb.r.gsd` <-
function(GSD,GSDo,level=NULL){
      if(ifelse(!is.null(level),level!=GSD$al,FALSE)){
                 GSD$al <- level;
                 GSD$b  <- compBounds(t=GSD$t/GSD$t[GSD$K], t2 = GSD$t*GSD$Imax, iuse = GSD$SF, asf = NULL, 
                                   alpha = GSD$al, phi =ifelse(is.null(GSD$phi),0,GSD$phi), 
                                   ztrun = 8)
      }
      res=(GSDo$z-GSD$b[GSDo$T])/sqrt(GSD$t[GSDo$T]*GSD$Imax)
      return(res)
      }

