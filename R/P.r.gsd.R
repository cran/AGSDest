`P.r.gsd` <-
function(h=0,GSD,GSDo,prec=0.001){
        K <- length(GSD$t);
        GSD$K<-length(GSD$t);
        if(ifelse(is.null(GSD$SF),TRUE,GSD$SF==5)) {
         b <- GSDo$z-h*sqrt(GSD$t[GSDo$T]*GSD$Imax);
         if(K==1){ 1-pnorm(b)
          } else {
            res=seqmon(a=GSD$a,
                   b=GSD$b*(b/GSD$b[GSDo$T]),
                   t=(GSD$t*GSD$Imax)/(GSD$t[K]*GSD$Imax),int=500*array(c(1),K)
                   )[2*K]
          }
       } else {
         if(GSD$SF==7){
            h1 <- (b-GSD$b[GSDo$T])/sqrt(GSD$t[GSDo$T]*GSD$Imax);
            res=seqmon(a=GSD$a,
                   b=GSD$b+h1*sqrt(GSD$t),
                   t=(GSD$t*GSD$Imax)/(GSD$t[K]*GSD$Imax),int=500*array(c(1),K)
                   )[2*K]
         } else {
           res=bisearch(function(u) GSDo$z-h*sqrt(GSD$t[GSDo$T]*GSD$Imax)-
                            compBounds(t=1:GSD$K/GSD$K, t2 = GSD$t*GSD$Imax, iuse = GSD$SF, asf = NULL, 
                                   alpha = u, phi = ifelse(is.null(GSD$phi),0,GSD$phi), 
                                   ztrun = 8)[GSDo$T],
                        c(0,1),signfl=-1, signfu=1,tol=prec)$root
         }
       }
      return(res)
      
      }

