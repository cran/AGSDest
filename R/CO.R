`CO` <-
function(h,pT,x=pT$b[length(pT$t)],iD){
         if(h==0){ k<-length(pT$t) } else {
            if(is.null(pT$alab)) pT$alab <-comp.alab(GSD=pT);
            k<-j.alab(h,GSD=pT)
         };
         if(k-iD$T<=0){ print("Cannot compute CO for T>=k"); 0 
         } else {
           if(x==Inf){ A(h,k,pT,iD)
           } else {
             ub <- pbounds(h=h,pT=list(t=pT$t[1:k],b=c(pT$b[1:(k-1)],x),Imax=pT$Imax),iD=iD);
             if(k-iD$T==1){ 1-pnorm(ub)
             } else {
               seqmon(a=pT$a[(iD$T+1):k],
                      b=ub[1:(k-iD$T)],t=pT$t[(iD$T+1):k]*pT$Imax-pT$t[iD$T]*pT$Imax,
                      int=500*array(c(1),k-iD$T))[2*(k-iD$T)]
             }
           }
         }
      }

