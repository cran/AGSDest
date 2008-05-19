`comp.als` <-
function(GSD){ 
               K <- length(GSD$t)
               if(K==1){ 1-pnorm(GSD$b)
               } else {
                 seqmon(a=GSD$a,b=GSD$b,t=GSD$t*GSD$Imax,int=500*array(c(1),K))[(K+1):(2*K)]
               }
             }

