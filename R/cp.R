`cp`<-function(GSD){
if(GSD$K>1)return(seqmon(a=GSD$a[1:GSD$K],b=pbounds(h=GSD$delta,pT=list(t=GSD$t,b=GSD$b,Imax=GSD$Imax),iD=list(T=0)),t=GSD$t[1:GSD$K],int=500*array(c(1),GSD$K))[2*GSD$K])
}
