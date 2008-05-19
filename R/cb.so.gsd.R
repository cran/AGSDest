`cb.so.gsd` <-
function(GSD,GSDo,level=NULL){
             GSD$K <- length(GSD$t)
             if(ifelse(!is.null(level),GSD$al!=level,FALSE)) {
                GSD$al<-level; GSD$alab<-comp.alab(GSD=GSD)  } else {
                if(is.null(GSD$alab)) GSD$alab <- comp.alab(GSD=GSD) 
             };
             if(GSDo$T==1) res=((GSDo$z-qnorm(1-GSD$al))/sqrt(GSD$t[GSDo$T]*GSD$Imax)) else {
               if(GSDo$z==GSD$b[GSDo$T]) res=(GSD$alab[GSDo$T]) else {             
                if(GSDo$z>GSD$b[GSDo$T]) hl<-GSD$alab[GSDo$T] else {
                   if(is.null(GSD$als)) GSD$als<-comp.als(GSD=GSD);
                   hl <- (GSDo$z-qnorm(1-GSD$al+GSD$als[GSD$K-1]))/sqrt(GSD$t[GSDo$T]*GSD$Imax)
                };                       
                res=uniroot( f=function(h){ seqmon(a=GSD$a[1:GSDo$T],
                                                b=pbounds(h=h,pT=list(t=GSD$t[1:GSDo$T],b=c(GSD$b[1:(GSDo$T-1)],GSDo$z),Imax=GSD$Imax)),
                                                t=GSD$t[1:GSDo$T]*GSD$Imax,int=500*array(c(1),GSDo$T))[2*GSDo$T]-GSD$al},
                                         interval=c(hl,GSD$alab[GSDo$T-1])
                          )$root
                         
                }
             }   
          return(res)
          }

