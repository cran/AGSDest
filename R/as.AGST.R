`as.AGST` <-
function(pT,iD,sT,sTo=NULL){
              if(is.null(sT)){
                            cat(paste("No design adaptation and hence no secondary trial; a GST object is created.\n\n"))
                            GST<-as.GST(GSD=pT,GSDo=iD)
                            return(GST)
                            }
              if(is.null(sTo))AGST<-list(pT=pT,iD=iD,sT=sT,sTo=list(T=0,z=0))
              else AGST<-list(pT=pT,iD=iD,sT=sT,sTo=sTo)
              class(AGST)<-"AGSTobj"
              return(AGST)             
              }

