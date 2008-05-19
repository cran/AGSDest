`A` <-
function(h,k=length(pT$t),pT,iD){
                  if(k-iD$T<=1) 0 else {
                      if(k-iD$T==2) 1-pnorm(pbounds(h,pT,iD)[1]) else {
                        seqmon(a=pT$a[(iD$T+1):(k-1)],
                               b=pbounds(h,pT,iD)[1:(k-1-iD$T)],
                               t=pT$t[(iD$T+1):(k-1)]*pT$Imax-pT$t[iD$T]*pT$Imax,
                               int=500*array(c(1),k-1-iD$T))[2*(k-1-iD$T)]
                     }                   
                   } 
              }

