`adapt` <-
function(pT,iD,SF,phi,cp,theta=iD$z/(pT$t[iD$T]*pT$Imax),I2min,I2max,swImax,delta=pT$delta,weight=0.5,exact=FALSE){
  cerror<-cer(pT=pT,iD=iD)
  K2min=round(I2min/swImax)
  K2max=round(I2max/swImax)
  if(!is.null(delta) & delta!=0){
    if(!is.null(weight))theta=(1-weight)*delta+weight*theta
    }
  
  if(theta<=0) {cat(paste("No design adaptation performed.\n"))
                sT=NULL
                class(sT)<-"GSTobj"
                return(sT)
                } 
  else {
       I2<-min(I2max,max(I2min,(qnorm(cp)-qnorm(cerror))^2/theta^2))
       if(K2max==K2min){K2<-K2max;} else {K2<-min(K2max,max(K2min,ceiling(I2/swImax)));}
                                   
         if(exact){
            sT=plan.GST(K=K2,SF=SF,phi=phi,alpha=cerror,delta=theta,pow=cp,compute.alab=FALSE,compute.als=FALSE)
            while(sT$Imax/K2>swImax & K2<8){
              K2=K2+1
              sT=plan.GST(K=K2,SF=SF,phi=phi,alpha=cerror,delta=theta,pow=cp,compute.alab=FALSE,compute.als=FALSE)
              }
            }
         else {
         sT=plan.GST(K=K2,Imax=I2,SF=SF,phi=phi,alpha=cerror,compute.alab=FALSE,compute.als=FALSE)
         sT$delta=theta
         } 
         class(sT)<-"GSTobj"
         return(sT)
     } 
}

