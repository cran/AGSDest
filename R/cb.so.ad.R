`cb.so.ad` <-
function(pT,iD,sT,sTo,level=NULL){

error=0;

#if(!is.loaded(symbol.C("AGSDest"))) {
#lib.file <- file.path(paste("AGSDest", .Platform$dynlib.ext, sep=""))
#dyn.load(paste("../library/AGSDest/libs/",lib.file, sep=""))
#cat(" -Loaded ", lib.file, "\n")
#}

if(is.null(level)){
  level=pT$al
  }
else{theta_1<-comp.alab(GSD=list(a=pT$a,b=pT$b,t=pT$t,al=level,Imax=pT$Imax))
    }

if(iD$T>length(pT$a)){print("iD$T > number stages pT");error=1}
if(sTo$T>length(sT$a)){print("sTo$T > number stages sT");error=1}

k_1<-length(pT$a);

if(is.null(pT$a)){print("pT$a is missing");error=1}
else a_1 <- pT$a
if(is.null(pT$b)){print("pT$b is missing");error=1}
else b_1 <- pT$b
if(is.null(pT$t)){print("pT$t is missing");error=1}
else t_1 <- pT$t*pT$Imax

if(is.null(pT$alab)){print("pT$alab is missing");error=1}
else theta_1<-pT$alab

k_2<-length(sT$a);
if(is.null(sT$a)){print("sT$a is missing");error=1}
else a_2 <- sT$a
if(is.null(sT$b)){print("sT$b is missing");error=1}
else b_2 <- sT$b
if(is.null(sT$t)){print("sT$t is missing");error=1}
else t_2 <- sT$t*sT$Imax

if(is.null(iD$T)){print("iD$T is missing");error=1}
else T_1 <- iD$T
if(is.null(iD$z)){print("iD$z is missing");error=1}
else z_1 <- iD$z

if(is.null(sTo$T)){print("sTo$T is missing");error=1}
else T_2 <- sTo$T
if(is.null(sTo$z)){print("sTo$z is missing");error=1}
else zT_2 <- sTo$z



e<-c(0,0);



out=0

if(error==0){

out <- .C(mainf,
          k = as.integer(k_1),
          a = as.numeric(a_1),
          b = as.numeric(b_1),
          t = as.numeric(t_1),
          alab = as.numeric(theta_1),
          al = as.numeric(level),
          T = as.integer(T_1),
          z = as.numeric(z_1),
          k2 = as.integer(k_2),
          a2 = as.numeric(a_2),
          b2 = as.numeric(b_2),
          t2 = as.numeric(t_2),
          T2 = as.integer(T_2),
          zT = as.numeric(zT_2),
          erg = as.numeric(e))
  return(out$erg[2]);
}
else {
  print("Values are missing.")
  return(0);
  }
}

