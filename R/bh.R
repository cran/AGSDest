`bh` <-
function(h,pT,xl=NULL,xu=NULL){
        if(h==0){ bh <-pT$b[length(pT$t)]
        } else {
         if(is.null(pT$alab)) pT$alab <-comp.alab(GSD=pT);
         if(h%in%pT$alab[1:(length(pT$t)-1)]){ Inf
          } else {
            k  <- j.alab(h,GSD=pT);
            if(is.null(xl)) xl <- qnorm(1-pT$al)+h*sqrt(pT$t[k]*pT$Imax);
            if((k==1)|(h<=pT$alab[length(pT$t)+1])){ xl
            } else {
                if(k<length(pT$t)) xl <-max(xl,pT$b[k]);
                if(is.null(xu)) {
                   al0 <-ifelse(k==2,1-pnorm(pT$b[(k-1)]-h*sqrt(pT$t[(k-1)]*pT$Imax)),
                                seqmon(a=pT$a[1:(k-1)],b=pT$b[1:(k-1)]-h*sqrt(pT$t[1:(k-1)]*pT$Imax),t=pT$t[1:(k-1)]*pT$Imax,
                                int=500*array(c(1),k-1))[2*(k-1)]);
                   al1 <-(pT$al-al0)/(1-al0)
                   xu  <-qnorm(1-al1)+h*sqrt(pT$t[k]*pT$Imax)
                }
                uniroot(f=function(x)  seqmon(a=pT$a[1:k],b=c(pT$b[1:(k-1)],x)-h*sqrt(pT$t[1:k]*pT$Imax),t=pT$t[1:k]*pT$Imax,
                                            int=500*array(c(1),k))[2*k]-pT$al, interval=c(xl,xu))$root
              }
            }
         }
        }

