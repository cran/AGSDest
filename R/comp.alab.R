`comp.alab` <-
function(GSD,pprec=1e-05){
              K <- length(GSD$a);
              alab <- numeric(K+1);
              for(k1 in 1:(K-1)){
                 if(k1==1){ alab[1] <-(GSD$b[1]-qnorm(1-GSD$al))/sqrt(GSD$t[1]*GSD$Imax)
                 } else {
                   alab[k1] <- uniroot(function(x) seqmon(a=GSD$a[1:k1],
                                                         b=GSD$b[1:k1]-x*sqrt(GSD$t[1:k1]*GSD$Imax),
                                                         t=GSD$t[1:k1]*GSD$Imax,
                                                         int=500*array(c(1),k1))[2*k1]-GSD$al,
                                       c(0,alab[k1-1]))$root
                 }
              }
              alab[K]   <-0;
              #lower critical boundary for delta's where the exceeding probability 
              #(estimated by Bonferroni) is smaller than pprec
              alab[K+1] <-min((GSD$b[1:(K-1)]-qnorm(1-pprec/(1:(K-1))))/sqrt(GSD$t[1:(K-1)]*GSD$Imax));
              alab[K+2] <-pprec;
              alab
             }

