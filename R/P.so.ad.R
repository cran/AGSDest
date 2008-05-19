`P.so.ad` <-
function(h=0,pT,iD,sT,sTo,prec=0.0001){
cer <- sapply((iD$T+1):length(pT$a), function(j){pT_new<-list(a=pT$a[1:j],b=pT$b[1:j],t=pT$t[1:j],Imax=pT$Imax);
                                          ifelse(j>iD$T,CO(h=0,pT=pT_new,iD=iD),0)})
P2=P.so.gsd(h=0,sT,sTo)
k=sum(ifelse(P2>cer,1,0))+iD$T+1

bkl=qnorm(1-pT$als[k])
bku=1-qnorm((pT$als[k]+pT$als[k-1])/2+pT$als[k-1])


pT_new=list(a=pT$a[1:k],b=pT$b[1:k],t=pT$t[1:k],al=pT$als[k],Imax=pT$Imax)
pT_new$b[k]=bkl
ceu=cerr(h=0,pT=pT_new,iD)
u=0

while(P2>ceu & u==0){
      pT_new$b[k]= pT_new$b[k]+(bkl+bku)/2
      bku=bkl
      bkl=bkl+(bkl+bku)/2
      ceu=cerr(h=0,pT=pT_new,iD)
      if(abs(sword(h=0,k,pT,bku)-pT$als[k+1])<=prec){u=sword(h=0,k,pT,bku);return(u)}
      }
bku=uniroot(f=function(bk) cerr(h=0,pT=list(a=pT$a[1:k],b=c(pT$b[1:k-1],bk),t=pT$t[1:k],al=pT$als[k],Imax=pT$Imax),iD)-P2,interval=c(bkl,bku))$root 
res=sword(h=0,k,pT,bku)

return(res)
}

