`cb.r.ad` <-
function(pT,iD,sT,sTo,level=NULL){
################################################################################
# sT$SF = 6, b_(j,u) according to stage wise ordering (exact second stage p-value; no flexible monitoring in secondary trial)
# sT$SF = 7, b_{j,u} according the ordering of the repeated confidence bound (not implemented yet)
################################################################################
          pT$K<-length(pT$t);
          sT$K<-length(sT$t);
          if(ifelse(!is.null(level),level!=pT$al,FALSE)){
                 pT$al <- level;
                 pT$b  <- compBounds(t=pT$t/pT$t[pT$K], t2 = pT$t*pT$Imax, iuse = pT$SF, asf = NULL, 
                                   alpha = pT$al, phi =ifelse(is.null(pT$phi),0,pT$phi), 
                                   ztrun = 8)#$upper.bounds;
                  sT$ce0 <- rCER(0,pT,iD);
                  if(sT$K>1){
                    sT$b  <-compBounds(t=sT$t/sT$t[sT$K], t2 = sT$t*sT$Imax, iuse = sT$SF, asf = NULL, 
                                   alpha = sT$ce0, phi =ifelse(is.null(sT$phi),0,sT$phi), 
                                   ztrun = 8)#$upper.bounds;
                   } else { sT$b <- qnorm(1-rCER(h=0,pT,iD)) }
          };
          
          if(sTo$z>=sT$b[sTo$T]){ hl<-0; hu<-(sTo$z-min(qnorm(1-pT$al),sT$b[sTo$T]))/sqrt(sT$t[sTo$T]*sT$Imax)
             }else{ hu<-0; hl<-(sTo$z-sT$b[sTo$T])/sqrt(sT$t[sTo$T]*sT$Imax) };
          if(sT$K==1) {
             res=uniroot(function(h){ sTo$z-h*sqrt(sT$t[sTo$T]*sT$Imax)-qnorm(1-rCER(h,pT,iD))
                                },c(hl,hu)
             )$root
          } else {
            if(ifelse(!is.null(sT$SF),sT$SF==6,FALSE)){
               res=uniroot(function(h){P.so.gsd(h,GSD=sT)-rCER(h,pT,iD)},c(hl,hu))$root
            } else  {
                 res=uniroot(function(h){
                       sTo$z-h*sqrt(sT$t[sTo$T]*sT$Imax)-
                         compBounds(t=sT$t/sT$t[sT$K], t2 = sT$t*sT$Imax, iuse = sT$SF, asf = NULL, 
                                   alpha = rCER(h,pT,iD), phi =ifelse(is.null(sT$phi),0,sT$phi), 
                                   ztrun = 8)[sTo$T]#$upper.bounds[sTo$T]                         
                       },c(hl,hu)
                  )$root
               }
            }
return(res)
}

