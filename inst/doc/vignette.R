### R code from vignette source 'vignette.Rnw'

###################################################
### code chunk number 1: vignette.Rnw:125-128
###################################################
library(AGSDest)
GSD<-plan.GST(K=3,SF=4,phi=-4,alpha=0.025,delta=5,pow=0.8,compute.alab=TRUE,compute.als=TRUE)
GSD


###################################################
### code chunk number 2: vignette.Rnw:164-166
###################################################
GST<-as.GST(GSD=GSD,GSDo=list(T=2, z=1.088))
GST


###################################################
### code chunk number 3: vignette.Rnw:170-171
###################################################
pvalue(GST,type="r")


###################################################
### code chunk number 4: vignette.Rnw:184-185
###################################################
pvalue(as.GST(GSD,list(T=2,z=2.63)),type='so')


###################################################
### code chunk number 5: vignette.Rnw:201-202
###################################################
seqconfint(GST,type='r')


###################################################
### code chunk number 6: vignette.Rnw:214-215
###################################################
seqconfint(as.GST(GSD,list(T=2,z=2.63)),type='so')


###################################################
### code chunk number 7: vignette.Rnw:223-224
###################################################
seqconfint(as.GST(GSD,list(T=2,z=2.63)),type="so",level=0.5)


###################################################
### code chunk number 8: vignette.Rnw:231-232
###################################################
seqconfint(as.GST(GSD,list(T=2,z=1.088)),type='r',level=0.5)


###################################################
### code chunk number 9: vignette.Rnw:261-264
###################################################
GSD1<-as.GST(GSD,list(T=2,z=2.63))
GSD1<-summary(GSD1,ctype='so',ptype='so',etype='n')
GSD1


###################################################
### code chunk number 10: vignette.Rnw:378-380
###################################################
pT=GSD
iD=list(T=1, z=0.731)


###################################################
### code chunk number 11: vignette.Rnw:384-385
###################################################
cer(pT,iD)


###################################################
### code chunk number 12: vignette.Rnw:389-390
###################################################
swImax=200/(4*20^2)


###################################################
### code chunk number 13: vignette.Rnw:394-400
###################################################
I2min=2*swImax
I2max=5*swImax

sT=adapt(pT=pT,iD=iD,SF=1,phi=0,cp=0.9,theta=4,I2min=I2min,I2max=I2max,swImax=swImax)
sT
AGSD<-as.AGST(pT,iD,sT)


###################################################
### code chunk number 14: vignette.Rnw:442-445
###################################################
sTo=list(T=2,z=1.532)
AGSD<-as.AGST(pT,iD,sT,sTo)
pvalue(AGSD,type='r')


###################################################
### code chunk number 15: vignette.Rnw:468-470
###################################################
AGSD1<-as.AGST(pT,iD,sT,list(T=3,z=2.73))
pvalue(AGSD1,type='so')


###################################################
### code chunk number 16: vignette.Rnw:494-495
###################################################
seqconfint(AGSD,type='r')


###################################################
### code chunk number 17: vignette.Rnw:525-526
###################################################
seqconfint(AGSD1,type='so')


###################################################
### code chunk number 18: vignette.Rnw:534-535
###################################################
seqconfint(AGSD1,type="so",level=0.5)


###################################################
### code chunk number 19: vignette.Rnw:541-542
###################################################
seqconfint(AGSD,type="r",level=0.5)


###################################################
### code chunk number 20: vignette.Rnw:571-573
###################################################
AGSD1<-summary(AGSD1,ctype='so',ptype='so',etype='a')
AGSD1


