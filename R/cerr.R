`cerr` <-
function(h=0,pT,iD,bhh=NULL){
          if(is.null(pT$alab)) pT$alab<-comp.alab(GSD=pT);  
          if(h==0){
              CO(h=0,pT=pT,iD=iD)
          } else {
            if(h%in%pT$alab[1:(length(pT$t)-1)]){ A(h,k=j.alab(h,GSD=pT)+1,pT,iD)
             } else {
             if(is.null(bhh)) bhh <- bh(h,pT);
             CO(h=h,pT=pT,x=bhh,iD=iD) }
          }
        }

