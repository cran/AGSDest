`compBounds` <-
function(t=1:length(t2)/length(t2), t2, iuse, asf = NULL,alpha, phi,ztrun){
#   if(whatSpendingFunctionIsUsed==1 || whatSpendingFunctionIsUsed==3){
   K<-length(t)
   n<-length(t)
   if(iuse==3){
      bmin <- qnorm(1-alpha);
      bmax <- qnorm(1-(alpha/n))/min(1,n^(phi-0.5));
      if(n>1){b <- uniroot(f=function(x)seqmon(a=rep(-ztrun,n),
                           b=x*(1:n)^(phi-0.5),
                           t=t2,
                           int=500*array(1,n))[2*n]-alpha
                           ,interval=c(bmin,bmax))$root;
        b<-b*(1:n)^(phi-0.5)
        }
      else b <- qnorm(1-alpha)
   }
   else {
   if((iuse==4)){
    if(phi == 0){print("phi must be unequal 0");}
    else{
    if(n>1){
     levsp <- alpha*(1-exp(-phi*t2/t2[n]))/(1-exp(-phi));
      b    <- numeric(n);
      b[1] <- qnorm(1-levsp[1]);
      for(i in 2:n){
        bmin <- qnorm(1-levsp[i]);
        bmax <- qnorm(1-levsp[i]+levsp[i-1]);
        if(abs(bmin-bmax)<=10^(-4)) b[i] <- bmin else {
             b[i] <- uniroot( function(x) seqmon(a=rep(-ztrun,i),
                                                 b=c(b[1:(i-1)],x),
                                                 t=t2[1:i],
                                                 int=500*array(1,i))[2*i]-levsp[i],
                                                 interval=c(bmin,bmax))$root
        }
      }
    } else  b <- qnorm(1-alpha)
   }
   } else {
     b=bounds(t=1:K/K, t2 = t2, iuse = iuse, asf = NULL,alpha = alpha, phi = phi,ztrun = 8)$upper.bounds
     }
   }
 return(b=b);
}

