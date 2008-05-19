`adapt` <-
function(pT,iD,SF,phi,cp,theta=iD$z/(pT$t[iD$T]*pT$Imax),I2min,I2max,swImax,delta=pT$delta,weight=1,warn=TRUE){
  cerror<-cer(pT=pT,iD=iD)
  K2min=ceiling(I2min/swImax)
  K2max=ceiling(I2max/swImax)
  if(!is.null(delta) & delta!=0){
    if(!is.null(weight))theta=(1-weight)*delta+weight*theta
    }
  
  if(theta<=0) {cat(paste("No design adaptation performed.\n"))
                sT=NULL
                class(sT)<-"GSTobj"
                return(sT)
                } 
  else {
	I2<-(qnorm(cp)-qnorm(cerror))^2/theta^2
	if(I2<=swImax){
		cat("\nI2<=swImax, 1 stage trial.\n")		
		sT=list(K=1,Imax=I2,SF=SF,phi=phi,t=1,a=-8,b=1-qnorm(cerror),alpha=cerror,compute.alab=FALSE,compute.als=FALSE)
		cpower=cp
		sT$delta=theta
	}
	else{
		K2<-K2min
		I2=K2*swImax
		while(K2<=K2max){
			sT=plan.GST(K=K2,SF=SF,phi=phi,alpha=cerror,Imax=I2,compute.alab=FALSE,compute.als=FALSE)			
			sT$delta=theta
			cpower=cp(sT)
			if(cpower>cp){
				sT=plan.GST(K=ifelse(K2>K2min,K2-1,K2),SF=SF,phi=phi,alpha=cerror,delta=theta,pow=cp,compute.alab=FALSE,compute.als=FALSE)	
				cpower=cp(sT)
				break
			}
			else if(cpower<=cp & K2>=K2max){
				if(warn)print("Power may be lower than planned.")
				break
			}
			else{
				K2=K2+1
				I2=K2*swImax
			}
		}
	}
	sT$cp=cpower
	class(sT)<-"GSTobj"
        return(sT)
  } 
}

