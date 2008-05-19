`cer.r` <-
function(u,K,pT,iD){
                    if(K-iD$T<=0){ print("Cannot compute CO for T>=K"); 0}
                    b=compBounds(t=1:K/K, t2 = pT$t*pT$Imax, iuse = pT$SF, asf = NULL,
                                   alpha = u, phi = ifelse(is.null(pT$phi),0,pT$phi),
                                   ztrun = 8)#$upper.bounds

                    ub <- pbounds(h=0,pT=list(t=pT$t[1:K],b=b,Imax=pT$Imax),iD=iD);
                    if(K-iD$T==1){ 1-pnorm(ub)
                    } else {
                      seqmon(a=pT$a[(iD$T+1):K],
                              b=ub[1:(K-iD$T)],t=pT$t[(iD$T+1):K]*pT$Imax-pT$t[iD$T]*pT$Imax,
                              int=500*array(c(1),K-iD$T))[2*(K-iD$T)]
                              }
                    }

