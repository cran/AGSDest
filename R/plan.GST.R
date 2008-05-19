`plan.GST` <-
function(K,t=(1:K)/K,Imax=NULL,SF,phi,alpha,delta=NULL,pow=NULL,compute.alab=TRUE,compute.als=TRUE){

if(!is.null(pow))if(pow>1)pow=pow/100

if(SF==4 & phi == 0){print("phi must be unequal 0");}
else{
if((!is.null(delta) & !is.null(pow)) | !is.null(Imax) ){

  if (!is.null(delta) & !is.null(pow) & !is.null(Imax))print("\nGST is planned only on the basis of Imax.\n");
  t=t[order(t)]/max(t)
  if(is.null(Imax))Imax=1
  b<-compBounds(t=t, t2 = t*Imax, iuse = SF, asf = NULL,
          alpha = alpha, phi = phi,ztrun = 8)#$upper.bounds
  a=rep(-8,K)
  
  if(!is.null(delta)){  
    hmin=qnorm(1-alpha)-qnorm(1-pow)
    hmax=b[K]-qnorm(1-pow)

    if(K==1) h=hmax
    else 
      h=uniroot(function(x) seqmon(a=a[1:K],b=pbounds(h=x,pT=list(t=t,b=b,Imax=Imax),iD=list(T=0)),t=t[1:K],
                                int=500*array(c(1),K))[2*K]-pow,c(hmin,hmax))$root          

    Imax=((h/delta)^2)
    }
  else delta=0
  b<-compBounds(t=(1:K)/K, t2 = Imax*t, iuse = SF, asf = NULL,
          alpha = alpha, phi = phi,ztrun = 8)

  GSD=list(K=K,al=alpha,t=t,SF=SF,phi=phi,a=a,b=b,Imax=Imax,delta=delta)
          
  if(compute.alab)GSD$alab=comp.alab(GSD=GSD)          
  if(compute.als)GSD$als=comp.als(GSD=GSD)           

  class(GSD)<-"GSTobj"
  return(GSD)
  }
else print("Please specify either delta and pow or Imax");
}
}

