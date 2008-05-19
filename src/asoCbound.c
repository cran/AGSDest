#include <R.h>             //specific R library
#include <Rdefines.h>      //specific R library
#include <Rmath.h>

#include <stdlib.h>
#include <stdio.h>

#define INF -100      //pre-defined value for minus infinity
#define pprec 0.0001  //pre-defined value for the precision


typedef double *var;

//Struct of the primary trial object
typedef struct{
    double level;       //level
    double *a;          //lower bound
    double *b;          //upper bound
    double *t;          //information
    int k;              //number of stages
    double *theta;      //alpha_globbings
    double bonf_algl;   //bonferroni estimate
    }pT_obj;

//Struct of the secondary trial object
typedef struct{
    double *a;          //lower bound
    double *b;          //upper bound
    double *t;          //information
    int k;              //number of stages
    double zT;          //z-statistic at stage T
    int T;              //stage where the trial stops
    //int OD;             //old design
    }sT_obj;

//Struct of the interim data object
typedef struct{
    int L;              //stage of the adaption
    double z;           //z-statistic at stage L
    }iD_obj;

//Struct of the H object
//describtion:  contains all variables for the calculation of the algorithm
typedef struct{
       double hl;       //lower interval bound
       double hu;       //upper interval bound
       double b1;       //lower estimation of the threshold boundary (8)
       double b2;       //upper estimation of the threshold boundary (8)
       double Al;       //value of A at the lower interval bound
       double Au;       //value of A at the upper interval bound
       double A;        //conditional rejection probability
       double B1;       //upper bound for Bk
       double B2;       //lower bound for Bk
       double p2;       //p-value of secondary trial
       int test;        //result of the test
       }H_list;

H_list initial_H(H_list H);

iD_obj initial(iD_obj iD);

var alloc_var(int sizei);

pT_obj alloc_seqmon_obj(pT_obj pT,int size_i);

double seqmon(double *a,double *b,double *t,int k,double theta);

double seqmon_b(double *a,double *b,double *t,int k,double theta);

double alpha_glob (double x, void *params);

int  alpha_globbing(pT_obj pT);

int  theta_interval(pT_obj pT,double h);

double cerror(double theta_c,pT_obj pT,iD_obj iD,int calc_b);

double epsilon(double theta_a,pT_obj pT,iD_obj iD,double x);

double B(double theta_a,pT_obj pT,iD_obj iD,double x);

double A(double theta_a,pT_obj pT,iD_obj iD);

sT_obj alloc_sT_obj(sT_obj s,int size_i);

double p2(sT_obj sT,double theta);

double sword(double theta,int k,pT_obj pT,double x);

void mainf(int *k,var a,var b,var t,var theta,double *level,int *T,double *z,
int *k2,var a2,var b2,var t2,int *T2,double *zT ,double *erg);


//all variables of the H_list object are initialised to 0. For detail of the
//H_list object see file functions.h
H_list initial_H(H_list H){
       H.hl=0;
       H.hu=0;
       H.b1=0;
       H.b2=0;
       H.Al=0;
       H.Au=0;
       H.A=0;
       H.B1=0;
       H.B2=0;
       H.p2=0;
       H.test=0;
       return H;
       }

//the iD_obj object is initialised to 0. For detail of the iD_obj object see
//file functions.h
iD_obj initial(iD_obj iD){
       iD.L=0;
       iD.z=0;
       return(iD);
       }

//memory is allocated for the variables of the secondary trial
sT_obj alloc_sT_obj(sT_obj s,int size_i){
    int i=0;
    s.a=alloc_var(size_i);
    s.b=alloc_var(size_i);
    s.t=alloc_var(size_i);
    //init to 0
    for(i=0;i<size_i;i++){
            s.a[i]=0;
            s.b[i]=0;
            s.t[i]=0;
            }
    return s;
    }


//********************************   p2   **************************************
//description:    Calculates the p-value of the secondary trial
//
//input:          sT           Object for the secondary trial
//                theta        Value of theta
//
//output:         erg          p-value of the secondary trial

double p2(sT_obj sT,double theta){
    int i=0;
    double erg=0;
    var b_new=alloc_var(sT.k);

    for(i=0;i<sT.k;i++)b_new[i]=sT.b[i];

    int new_T=sT.T;
    new_T--;
    if(sT.zT<b_new[new_T] && sT.T<sT.k){Rprintf("Error: zT < b[T]");return 0;}
    else {
         //the secondary trial stops at the first stage
         if(sT.T==1){
                     free(b_new);
                     return 1-pnorm(sT.zT-theta*sqrt(sT.t[0]),0,1,1,0);
                     }
         else{
              b_new[new_T]=sT.zT;
              erg=seqmon(sT.a,b_new,sT.t,sT.T,theta);
              }
         }
    free(b_new);
    return erg;
    }

//*****************************   theta_interval   *****************************
//description:    Calculates which interval contains theta
//
//input:          pT       Object of the primary trial
//                h        Value of theta
//
//output:         i        number of the interval that contains theta

int  theta_interval(pT_obj pT,double h){
  int i;
  for(i=0;i<pT.k;i++){
          if(h>=pT.theta[i]/*+0.00001*/){i++;return i;} //
          }
  return i;
}


//***********************************   B   ************************************
//description:    Calculates the value of B(h)
//
//input:          theta_a      Value of theta
//                pt           Object for the primary trial
//                iD           Object for the interim data
//                x            Value of the last critical boundary
//
//output:         B            value of B(h)

double B(double theta_a,pT_obj pT,iD_obj iD,double x){
int i=0;
double B=0;
int interval=theta_interval(pT,theta_a);

int last=interval;
last--;
pT_obj pT_new;

pT_new.k=interval;

pT_new=alloc_seqmon_obj(pT_new,interval);
//  if (!pT_new.a || !pT_new.b|| !pT_new.t|| !pT_new.theta )
//    {
//      Rprintf("No memory free for double-components (B)\n");
//      exit(1);
//    }

for(i=0;i<interval;i++){
        pT_new.b[i]=pT.b[i];
        pT_new.a[i]=pT.a[i];
        pT_new.t[i]=pT.t[i];
        if(i==last)pT_new.b[i]=x;
        }

B=cerror(theta_a,pT_new,iD,1);

//free memory
free(pT_new.a);free(pT_new.b);free(pT_new.t);free(pT_new.theta);

return B;
}

//***********************************   A   ************************************
//description:      Calculates the value of A(h)
//
//input:            theta_a      Value of theta
//                  pt           Object for the primary trial
//                  iD           Object for the interim data
//
//output:           A            value of A(h)

double A(double theta_a,pT_obj pT,iD_obj iD){
int i=0;
double A=0;
int interval=theta_interval(pT,theta_a);

if(interval-iD.L<=1)return 0;

else{

interval--;

pT_obj pT_new;
pT_new.k=interval;

pT_new=alloc_seqmon_obj(pT_new,interval);
  if (!pT_new.a || !pT_new.b|| !pT_new.t|| !pT_new.theta )
    {
      Rprintf("No memory free for double-components (A)\n");
      exit(1);
    }

for(i=0;i<interval;i++){
        pT_new.b[i]=pT.b[i];
        pT_new.a[i]=pT.a[i];
        pT_new.t[i]=pT.t[i];

        }


A=cerror(theta_a,pT_new,iD,0);

//free memory
free(pT_new.a);free(pT_new.b);free(pT_new.t);free(pT_new.theta);

return A;
}
}

//*******************************   cerror   ***********************************
//description:    Calculates the conditional error probability
//
//input: theta_c    Value of theta
//       pT         Object for the primary Trial
//       iD         Object for the interim Data
//       calc_b     if 0   claculate seqmon
//                  if 1   claculate seqmon_b
//
//output: e         conditional error probability

double cerror(double theta_c,pT_obj pT,iD_obj iD,int calc_b){
int i=0;
double e=0;
int k_new=pT.k;
k_new-=iD.L;
int L=iD.L;
L--;

if(k_new<1){
            Rprintf("cannot compute cerror for L>K ");
            return 0;
            }

var a_new=alloc_var(k_new);
var b_new=alloc_var(k_new);
var t_new=alloc_var(k_new);

int c=0;
for(i=iD.L;i<pT.k;i++){
        b_new[c]=(pT.b[i]*sqrt(pT.t[i])-iD.z*sqrt(pT.t[L]))/sqrt(pT.t[i]-pT.t[L]);
        a_new[c]=pT.a[i];
        t_new[c]=pT.t[i]-pT.t[L];
        c++;
        }

if(!calc_b)
e=seqmon(a_new,b_new,t_new,k_new,theta_c);
else
e=seqmon_b(a_new,b_new,t_new,k_new,theta_c);

//free memory
free(a_new);free(b_new);free(t_new);

return e;
}


//******************************   alloc_var   *********************************
//description:      Allocates unused space for an array
//
//input:            sizei     size of the allocated array
//
//output:           V         the new allocated array

var alloc_var(int sizei){
int i=0;
var V = (double *) calloc(sizei , sizeof(double));
//init to 0
for(i=0;i<sizei;i++){
            V[i]=0;
            }
return V;
}


//**************************   alloc_seqmon_obj   ******************************
//description:      Allocates unused space for all variables of a pT_obj
//
//input:            pT        Object for the primary Trial
//                  sizei     size of the allocated array
//
//output:           pT        A pT_obj with the allocated variables

pT_obj alloc_seqmon_obj(pT_obj pT,int size_i){
    int i=0;
    int newsize=size_i+1;
    pT.a=alloc_var(size_i);
    pT.b=alloc_var(size_i);
    pT.t=alloc_var(size_i);
    pT.theta=alloc_var(newsize);
    //init to 0
    for(i=0;i<size_i;i++){
            pT.a[i]=0;
            pT.b[i]=0;
            pT.t[i]=0;
            pT.theta[i]=0;
            }
    pT.theta[newsize-1]=0;
    return pT;
    }


//*******************************   seqmon   ***********************************
//description:        computes the probabilities of crossing boundaries in
//                    a group sequential clinical trial. It implements the
//                    Armitage-McPherson and Rowe (1969) algorithm using
//                    the method described in Schoenfeld D. (2001).
//
//input:              a          lower boundaries
//                    b          upper boundaries
//                    t          information
//                    k          number of stages
//                    theta      value of theta
//
//output:             pU_last    probability of crossing the boundaries

double seqmon(double *a,double *b,double *t,int k,double theta){

  int c=1;
  int interval=500;
  var M=alloc_var(interval);
  var Ma=alloc_var(interval);
  var Maf=alloc_var(interval);
  var x0=alloc_var(interval);
  var d=alloc_var(k);
  var pU=alloc_var(k);
  var pL=alloc_var(k);
  var VU=alloc_var(interval);
  var VL=alloc_var(interval);
  var x=alloc_var(interval);
  int i=0;
  int j=0;
  int w=0;

  var b_new=alloc_var(k);

  if (!x || !VU|| !VL|| !x0 || !Maf ||!Ma||!M ||!d||!pU||!pL||!b_new)
    {
      Rprintf("No memory free for double-components (seqmon)\n");
      exit(1);
    }

  double pUd,pLd,VUf,VLf,VUM=0,VLM=0;

	for(i=0;i<k;i++){
     if(theta) b_new[i]=b[i]-theta*sqrt(t[i]);
     else b_new[i]=b[i];
     d[i]=(b_new[i]-a[i])/(interval);
     }

  c=1;

  for(i=0;i<interval;i++){
      x0[i]=a[0]+(c-0.5)*(d[0]);
      M[i]=(d[0]/sqrt(2*M_PI))*exp(-(pow((sqrt(t[0])*x0[i]),2))/(2*t[0]));
      c++;
      }

  pUd=-(sqrt(t[0])*b_new[0])/sqrt(t[0]);
  //pnorm: cumulative distribution functions
  //For details see Rmath.h
  pU[0]=pnorm(pUd,0,1,1,0);


  pLd=(sqrt(t[0])*a[0])/sqrt(t[0]);
  pL[0]=pnorm(pLd,0,1,1,0);

  int count=0;
  int last=k;
  last--;
  for(i=1;i<k;i++){
    c=1;
     for(j=0;j<interval;j++){
         VUf=-(sqrt(t[i])*b_new[i] + (-sqrt(t[count])*x0[j]))/sqrt(t[i]-t[count]);
         VU[j]=pnorm(VUf,0,1,1,0);
         VLf=(sqrt(t[i])*a[i] + (-sqrt(t[count])*x0[j]))/sqrt(t[i]-t[count]);
         VL[j]=pnorm(VLf,0,1,1,0);
         VLM+=VL[j]*M[j];
         VUM+=VU[j]*M[j];
         x[j]=a[i] + (c - 0.5 ) * d[i];
         c++;
         }


     pL[i]=pL[count]+VLM;
     pU[i]=pU[count]+VUM;

     VLM=0.0;
     VUM=0.0;
   if(i!=last){
     for(w=0;w<interval;w++){
         for(j=0;j<interval;j++){
            Ma[j]=(d[i]*sqrt(t[i])/(sqrt(2.0*M_PI)*sqrt(t[i] - t[count])))*
                    expl(-pow((sqrt(t[i]) * x[w] - sqrt(t[count]) * x0[j]),2)/(2 * (t[i] - t[count])));
         	Maf[w]+=Ma[j]*M[j];
            }
     		}
     for(j=0;j<interval;j++){
         M[j]=Maf[j];
         x0[j]=x[j];
         Maf[j]=0.0;
         Ma[j]=0.0;
         }
     }
     count++;
   }



for(i=0;i<k;i++){
                 if(!pU[i] || !pL[i]){
                           free(M);free(Ma);free(Maf);free(x0);free(VU);free(VL);free(x);free(d);
                           free(pU);free(pL);free(b_new);
                           return -1;
                           }
                 }

double pU_last=pU[last];

free(M);free(Ma);free(Maf);free(x0);free(VU);free(VL);free(x);free(d);
free(pU);free(pL);free(b_new);

return pU_last;

}


//*******************************   seqmon_b   *********************************
//description:        It's a modified version of seqmon and calculates the
//                    probability of crossing only the last upper boundary
//
//input:              a          lower boundaries
//                    b          upper boundaries
//                    t          information
//                    k          number of stages
//                    theta      theta
//
//output:             pU_last    probability of crossing the last boundary


double seqmon_b(double *a,double *b,double *t,int k,double theta){

  int c=1;
  int interval=500;
  var M=alloc_var(interval);
  var Ma=alloc_var(interval);
  var Maf=alloc_var(interval);
  var x0=alloc_var(interval);
  var d=alloc_var(k);
  var VU=alloc_var(interval);
  var x=alloc_var(interval);
  var b_new=alloc_var(k);

  double pU,pUd;
  int i=0;
  int j=0;
  int w=0;

  if (!x || !VU|| !x0 || !Maf ||!Ma||!M ||!d||!b_new)
    {
      Rprintf("No memory free for double-components (seqmon_b)\n");
      exit(1);
    }

  double VUf=0,VUM=0;

	for(i=0;i<k;i++){
     if(theta) b_new[i]=b[i]-theta*sqrt(t[i]);
     else b_new[i]=b[i];
     d[i]=(b_new[i]-a[i])/(interval);
     }

  c=1;



  for(i=0;i<interval;i++){
      x0[i]=a[0]+(c-0.5)*(d[0]);
      M[i]=(d[0]/sqrt(2*M_PI))*exp(-(pow((sqrt(t[0])*x0[i]),2))/(2*t[0]));
      c++;
      }

  pUd=-(sqrt(t[0])*b_new[0])/sqrt(t[0]);
  pU=pnorm(pUd,0,1,1,0);

  int count=0;
  int last=k;

  if(k==1){
           free(M);free(Ma);free(Maf);free(x0);free(VU);free(x);free(d);
           free(b_new);
           return pU;
           }

  last--;
  for(i=1;i<k;i++){
     c=1;
     for(j=0;j<interval;j++){
          x[j]=a[i] + (c - 0.5 ) * d[i];
          c++;
          }
   if(i!=last){
     for(w=0;w<interval;w++){
         for(j=0;j<interval;j++){
            Ma[j]=(d[i]*sqrt(t[i])/(sqrt(2.0*M_PI)*sqrt(t[i] - t[count])))*
                    exp(-pow((sqrt(t[i]) * x[w] - sqrt(t[count]) * x0[j]),2)/(2 * (t[i] - t[count])));
         	Maf[w]+=Ma[j]*M[j];
            }
     		}
     for(j=0;j<interval;j++){
         M[j]=Maf[j];
         x0[j]=x[j];
         Maf[j]=0.0;
         Ma[j]=0.0;
         }
     }
   if(i==last){
     for(j=0;j<interval;j++){
         VUf=-(sqrt(t[i])*b_new[i] + (-sqrt(t[count])*x0[j]))/sqrt(t[i]-t[count]);
         VU[j]=pnorm(VUf,0,1,1,0);
         VUM+=VU[j]*M[j];
         }
     pU=VUM;
     }
     count++;
   }


VUf=VUM=0;
free(M);free(Ma);free(Maf);free(x0);free(VU);free(x);free(d);
free(b_new);

if(!pU)return -1;
else return pU;

}








//*******************************   seqmonc   ***********************************
//description:        computes the probabilities of crossing boundaries in
//                    a group sequential clinical trial. It implements the
//                    Armitage-McPherson and Rowe (1969) algorithm using
//                    the method described in Schoenfeld D. (2001).
//
//input:              a          lower boundaries
//                    b          upper boundaries
//                    t          information
//                    k          number of stages
//                    theta      value of theta
//
//output:             pU_last    probability of crossing the boundaries

void seqmonc(double *a,double *b,double *t,int k,double theta,double *erg){

  int c=1;
  int interval=500;
  var M=alloc_var(interval);
  var Ma=alloc_var(interval);
  var Maf=alloc_var(interval);
  var x0=alloc_var(interval);
  var d=alloc_var(k);
  var pU=alloc_var(k);
  var pL=alloc_var(k);
  var VU=alloc_var(interval);
  var VL=alloc_var(interval);
  var x=alloc_var(interval);
  int i=0;
  int j=0;
  int w=0;

  var b_new=alloc_var(k);

  if (!x || !VU|| !VL|| !x0 || !Maf ||!Ma||!M ||!d||!pU||!pL||!b_new)
    {
      Rprintf("No memory free for double-components (seqmon)\n");
      exit(1);
    }

  double pUd,pLd,VUf,VLf,VUM=0,VLM=0;

	for(i=0;i<k;i++){
     if(theta) b_new[i]=b[i]-theta*sqrt(t[i]);
     else b_new[i]=b[i];
     d[i]=(b_new[i]-a[i])/(interval);
     }

  c=1;

  for(i=0;i<interval;i++){
      x0[i]=a[0]+(c-0.5)*(d[0]);
      M[i]=(d[0]/sqrt(2*M_PI))*exp(-(pow((sqrt(t[0])*x0[i]),2))/(2*t[0]));
      c++;
      }

  pUd=-(sqrt(t[0])*b_new[0])/sqrt(t[0]);
  //pnorm: cumulative distribution functions
  //For details see Rmath.h
  pU[0]=pnorm(pUd,0,1,1,0);


  pLd=(sqrt(t[0])*a[0])/sqrt(t[0]);
  pL[0]=pnorm(pLd,0,1,1,0);

  int count=0;
  int last=k;
  last--;
  for(i=1;i<k;i++){
    c=1;
     for(j=0;j<interval;j++){
         VUf=-(sqrt(t[i])*b_new[i] + (-sqrt(t[count])*x0[j]))/sqrt(t[i]-t[count]);
         VU[j]=pnorm(VUf,0,1,1,0);
         VLf=(sqrt(t[i])*a[i] + (-sqrt(t[count])*x0[j]))/sqrt(t[i]-t[count]);
         VL[j]=pnorm(VLf,0,1,1,0);
         VLM+=VL[j]*M[j];
         VUM+=VU[j]*M[j];
         x[j]=a[i] + (c - 0.5 ) * d[i];
         c++;
         }


     pL[i]=pL[count]+VLM;
     pU[i]=pU[count]+VUM;

     VLM=0.0;
     VUM=0.0;
   if(i!=last){
     for(w=0;w<interval;w++){
         for(j=0;j<interval;j++){
            Ma[j]=(d[i]*sqrt(t[i])/(sqrt(2.0*M_PI)*sqrt(t[i] - t[count])))*
                    expl(-pow((sqrt(t[i]) * x[w] - sqrt(t[count]) * x0[j]),2)/(2 * (t[i] - t[count])));
         	Maf[w]+=Ma[j]*M[j];
            }
     		}
     for(j=0;j<interval;j++){
         M[j]=Maf[j];
         x0[j]=x[j];
         Maf[j]=0.0;
         Ma[j]=0.0;
         }
     }
     count++;
   }



for(i=0;i<k;i++){
                 if(!pU[i] || !pL[i]){
                           free(M);free(Ma);free(Maf);free(x0);free(VU);free(VL);free(x);free(d);
                           free(pU);free(pL);free(b_new);
                           //return -1;
                           break;
                           }
                 }

//double pU_last=pU[last];

for(i=0;i<k;i++)
erg[i]=pL[i];
for(i=0;i<k;i++)
erg[i+k]=pU[i];


free(M);free(Ma);free(Maf);free(x0);free(VU);free(VL);free(x);free(d);
free(pU);free(pL);free(b_new);

}












//****************************sword*****************************************
//description:    Calculates condtitonal rejection probability with last
//                critical boundary equal x
//
//input:          theta        value of theta
//                k            number of stages
//                pt           object for the primary trial
//                x            value of the last critical boundary
//
//output:         erg


double sword(double theta,int k,pT_obj pT,double x){
       int i=0;
       int new_k=k;
       double erg=0;
       new_k--;
         if(k==1){ return 1-pnorm(x-theta*sqrt(pT.t[new_k]),0,1,1,0);}
         else {
              var b_new=alloc_var(k);
              var t_new=alloc_var(k);

              for(i=0;i<k;i++){
                        b_new[i]=pT.b[i];
                        t_new[i]=pT.t[i];
                        }
              b_new[new_k]=x;

              erg=seqmon(pT.a,b_new,t_new,k,theta);
              free(b_new);
              free(t_new);
              return erg;
               }
}








/******************************************************************************/
//The function testint calculates if the tested intervals is rejected, accepted
//or that there is no decision possible.
//input: H      obejct of H_list
//       pT     object of the primary trial
//       sT     object of the secondary trial
//       iD     object of the interim data
//output: H
//        the output of the test is saved in
//        H.test(2=reject; 1=accept; 0=no decision)
/******************************************************************************/

H_list testint(H_list H,pT_obj pT,iD_obj iD,sT_obj sT){
    int intervalhl=0;
    int intervalhu=0;
    int interval=0;
    int last=0;
    int check_a=0;
    int check_b=0;
    int i=0;

    //theta_interval calculates the stage containing hl and hu
    intervalhl=theta_interval(pT,H.hl);
    intervalhu=theta_interval(pT,H.hu);

    //if hl and hu are not in the same interval we have to abort
    if(intervalhl!=intervalhu && H.hl>0 && intervalhu>iD.L){
    Rprintf("hl and hu not in the same alpha-absorbing interval or interval < iD L");
    H.test=-1;
    return H;
    }

    last=intervalhu;
    last--;

    //Value of A at the upper bound hu
    if(H.Au==0){H.Au=A(H.hu,pT,iD);}

    //P-value of the secondary trial
    if(H.p2==0)H.p2=p2(sT,H.hu);
//****************************************************************************
//parameter values so small that we can neglect rejection boundaries from the
//first K-1 stages
    if(H.hu==pT.bonf_algl){
       //qnorm inverse cumulative distribution functions.
       //For details see GNU Scientific Library
       H.b1=H.b2=qnorm(1-pT.level,0,1,1,0)+H.hu*sqrt(pT.t[last]);
       if(!H.B1)H.B1=B(H.hu,pT,iD,H.b1);
      if(H.p2<=H.B1){H.test=2;}
      else{
           if(H.p2>H.B1+H.Au)H.test=1;
           else H.test=0;
           }
      return H;
      }

//****************************************************************************
//parameter values not small enough that we can neglect rejection boundaries
    else {
         //Value of A at the lower bound hu
         if(!H.Al){
                   if(H.hl==INF)H.Al=0; //if hl=-inf
                   else H.Al=A(H.hl,pT,iD);
                   }

         interval=last;
         interval--;

//****************************************************************************
//hu is equal to the upper boundary of the tested interval

    if(fabs(pT.theta[interval]-H.hu)<0.00001 && H.hu>0.00001)
      {
      if(H.p2<=H.Al){H.test=2; // reject interval
            }
      else {
           if(H.p2<=H.Au){H.test=1;//accept interval
                }
           else {H.test=0;//no decision possible
                }
           }
      return H;
      }
//***************************************************************************

         var a_new=alloc_var(last);
         var b_new=alloc_var(last);
         var t_new=alloc_var(last);

         if(!H.A){
                  for(i=0;i<last;i++){
                        a_new[i]=pT.a[i];
                        b_new[i]=pT.b[i];
                        t_new[i]=pT.t[i];
                        }
                 H.A=seqmon(a_new,b_new,t_new,last,H.hu);
             }

         free(a_new);
         free(b_new);
         free(t_new);



//lower and upper estimations for the boundary bk(h) as in (8)
if(!H.b1)H.b1=qnorm(1-pT.level,0,1,1,0)+H.hu*sqrt(pT.t[last]);
if(!H.b2)H.b2=qnorm(1-((pT.level-H.A)/(1-H.A)),0,1,1,0)+H.hu*sqrt(pT.t[last]);

//Values of B at the upper and lower bound b1 and b2
if(!H.B1)H.B1=B(H.hu,pT,iD,H.b1);
if(!H.B2)H.B2=B(H.hu,pT,iD,H.b2);


//Algorithm
for(i=0;i<100;i++){

//step (S1)
if(!check_a){
    //case (a)
    if(H.p2<=H.Al+H.B2){
                        H.test=2;//reject interval
                        check_a=1;
                        break;
                        }
    //case (b)
    if(H.p2>H.Al+H.B1){
                            check_a=1;//exclude case (a) and pass
                                      //directly to (S2)
                            }
    //case (c)
    if(H.Al+H.B2<H.p2 && H.p2<=H.Al+H.B1){
                            check_a=0;//Without excluding case (a), pass
                                      //directly to (S2)
                            }
    }
//step (S2)
if(!check_b){
         //case (a)
         if(H.p2>H.Au+H.B1){
                            check_b=1;
                            H.test=1;//accept the interval
                            break;
                            }
         //case (b)
         if(H.Al+H.B1<H.p2 && H.p2<=H.Au+H.B2){//refine interval
                           H.test=0;
                           check_b=1;
                           break;
                           }
         else{
              //calculate rejection probability
              if(sword(H.hu,intervalhl,pT,(H.b1+H.b2)/2)>pT.level){
                  H.b1=(H.b1+H.b2)/2;
                  H.B1=B(H.hu,pT,iD,H.b1);
                  }
              else{
                  H.b2=(H.b1+H.b2)/2;
                  H.B2=B(H.hu,pT,iD,H.b2);
                  }
         }
}
}
if(!H.test && i>98){
          //maximum number of 100 iterations reached without convergence
          Rprintf("maximum number of 100 iterations reached without convergence\n");
          H.test=-1;
          return H;
          }
else return H;
}

}



//description: The function main calculates from a given group sequential plan
//             with adaptive design the associated intervals, passes them to
//             the function testint and returns the lower confidence bound
//
//input:   parameters for the primary trial:
//             k          number of stages
//             a          lower boundaries
//             b          upper boundaries
//             t          information
//             theta      theta
//             level      alpha
//         parameters for the interim data
//             L          stage of the adaption
//             z          z-statistic at stage L
//         parameters for the secondary trial
//             k2         number of stages
//             a2         lower boundaries
//             b2         upper boundaries
//             t2         information
//             T2         stage where the trial stops
//             zT         z-statistic at stage T2
//
//output:      erg        lower bound of the calculated confidence interval

void mainf(int *k,var a,var b,var t,var theta,double *level,int *T,double *z,
int *k2,var a2,var b2,var t2,int *T2,double *zT ,double *erg){


int i=0;

int *L=T;

pT_obj pT;//declare primary trial object pT
iD_obj iD;//declare interim data object iD
iD=initial(iD);
sT_obj sT;//declare secondary trial object sT

pT.k=k[0];//stages pT

pT=alloc_seqmon_obj(pT,pT.k);
//  if (!pT.a || !pT.b|| !pT.t|| !pT.theta )
//    {
//      Rprintf("No memory available for double-components!\n");
//      exit(1);
//    }

pT.theta=alloc_var(pT.k+1);
for(i=0;i<pT.k;i++){
        pT.a[i]=a[i];//lower boundaries pT
        pT.b[i]=b[i];//upper boundaries pT
        pT.t[i]=t[i];//information pT
        pT.theta[i]=theta[i];//alpha_globbings
        }


sT.k=k2[0];//stages sT
sT=alloc_sT_obj(sT,sT.k);
  if (!sT.a || !sT.b|| !sT.t )
    {
      Rprintf("No memory available for double-components!\n");
      exit(1);
    }

for(i=0;i<sT.k;i++){
        sT.a[i]=a2[i];//lower boundaries pT
        sT.b[i]=b2[i];//upper boundaries pT
        sT.t[i]=t2[i];//information sT
        }


iD.L=L[0];//stage of the design adaptations

iD.z=z[0];//z-statistic at stage L

sT.zT=zT[0];//z-statistic at stage T2

sT.T=T2[0];//stage where the secondary trials stops

pT.level=level[0];//alpha


int j=0;
double min_bonf=0;
int count =0;
int start=pT.k;
start-=2;

int reject=0;
int thetac=pT.k;
thetac--;
int interval=0;
int testit=0;
int old_interval=0;


//initial object Hl
H_list Hl[pT.k][100];
for(i=0;i<pT.k;i++){
        for(j=0;j<100;j++){
        Hl[i][j]=initial_H(Hl[i][j]);
        }
}


double malgl=0;//mean of the alpha-globbing constants (6)
for(i=0;i<pT.k;i++){
        malgl+=pT.theta[i];
        }
malgl=malgl/pT.k;



//*****************************   Bonferroni estimate  *************************

pT.bonf_algl=100;

for(i=0;i<pT.k-1;i++){
        min_bonf=(pT.b[i]-qnorm(1-pprec/i,0,1,1,0))/sqrt(pT.t[i]);
        if(min_bonf<=pT.bonf_algl)pT.bonf_algl=min_bonf;
        }

//******************************************************************************


for(i=start;i>=0;i--){

count=0;
interval=0;

//bounds for the first interval
if(i==start){
     Hl[i][0].hl=INF;
     Hl[i][0].hu=pT.theta[i]-0.00001;
     }
//bounds for the rest of the intervals
else{
     Hl[i][interval].hl=pT.theta[thetac]+0.00001;
     Hl[i][interval].hu=pT.theta[i]-0.00001;
     }


Hl[i][interval]=testint(Hl[i][interval],pT,iD,sT);
count++;

//Rprintf("1 hl: %f, hu: %f, test: %d \n",Hl[i][interval].hl,Hl[i][interval].hu,Hl[i][interval].test);

if(Hl[i][interval].test==2){//the whole interval can be rejected
   //Rprintf("1 hl: %f, hu: %f, test: %d \n",Hl[i][interval].hl,Hl[i][interval].hu,Hl[i][interval].test);
  }
else {//the interval is accepted or there is no decision possible
  //if(Hl[i][interval].test==-1){Rprintf("-1");break;}
  old_interval=interval;
  interval++; //increment list-counter
  Hl[i][interval].Al=Hl[i][old_interval].Al;
  Hl[i][interval].hl=Hl[i][old_interval].hl;
  if(i==start)Hl[i][interval].hu=Hl[i][old_interval].hu-malgl;//only in the case
                                                     //of the first interval (7)
  else Hl[i][interval].hu=(Hl[i][old_interval].hl+Hl[i][old_interval].hu)/2;
  Hl[i][interval]=testint(Hl[i][interval],pT,iD,sT);//new subinterval is tested
  count++;

  //Rprintf("2 hl: %f, hu: %f, test: %d, int: %d \n",Hl[i][interval].hl,Hl[i][interval].hu,Hl[i][interval].test,interval);

  reject=0;

  //continues as long as we don't have the pre-defined precision
  while(fabs(Hl[i][interval].hu-Hl[i][interval].hl)>=pprec &&
        testit!=1){

        if(Hl[i][interval].test==2){ //reject the interval
             reject++;
             if(reject>1){
                          old_interval=interval;
                          old_interval-=reject;
                          }
             else{
             old_interval=interval;
             old_interval--;
             }

             //in the case of the first interval
             if(i==start){
                          Hl[i][interval].hl=Hl[i][interval].hu;
                          Hl[i][interval].hu=Hl[i][old_interval].hu;
                          Hl[i][interval].p2=0;
                          Hl[i][interval].Al=0;
                          Hl[i][interval].Au=0;
                          Hl[i][interval].B1=0;
                          Hl[i][interval].B2=0;
                          Hl[i][interval].b1=0;
                          Hl[i][interval].b2=0;
                          }
             //in all other cases we can save a couple of values
             //to save computation time
             else {
                  Hl[i][interval].hl=Hl[i][interval].hu;
                  Hl[i][interval].hu=Hl[i][old_interval].hu;
                  Hl[i][interval].B1=Hl[i][old_interval].B1;
                  Hl[i][interval].B2=Hl[i][old_interval].B2;
                  Hl[i][interval].b1=Hl[i][old_interval].b1;
                  Hl[i][interval].b2=Hl[i][old_interval].b2;
                  Hl[i][interval].Au=Hl[i][old_interval].Au;
                  Hl[i][interval].Al=0;
                  Hl[i][interval].p2=0;
                  }
             }
        else{// akzept or no decision
            if(reject>0)reject--;
             old_interval=interval;
            interval++;//increment list-counter


            //in the case of the first interval
            if(i==start && Hl[i][interval].hl==INF){
                 Hl[i][interval].hu=Hl[i][old_interval].hu-malgl;
                 Hl[i][interval].p2=0;
                 }
            else{
            Hl[i][interval].hl=Hl[i][old_interval].hl;
            Hl[i][interval].hu=(Hl[i][old_interval].hl+
                                Hl[i][old_interval].hu)/2;
            Hl[i][interval].Al=Hl[i][old_interval].Al;
            Hl[i][interval].p2=0;
            }
        }


        //in the case of a new interval
        if(Hl[i][interval].hl==Hl[i][interval].hu &&
            Hl[i][interval].test==2 && Hl[i][interval].hu==pT.theta[i]){
                                     break;
                                     }

        //test new defined interval
        Hl[i][interval]=testint(Hl[i][interval],pT,iD,sT);
        count++;

        //Rprintf("3 hl: %f, hu: %f, test: %d, int: %d \n",Hl[i][interval].hl,Hl[i][interval].hu,Hl[i][interval].test,interval);


        if(fabs(pT.theta[i]-Hl[i][interval].hu)<0.0001 &&
            Hl[i][interval].test==2)break;

        //if interval is accepted and the pre-defined precision is
        //reached
        if(Hl[i][interval].test!=2 && fabs(Hl[i][interval].hu-Hl[i][interval].hl)<pprec)
             testit=1;
        else testit=0;


  }//end while
  if(fabs(Hl[i][interval].hu-Hl[i][interval].hl)<pprec){
                                                  break;
                                                 }

}//end else

thetac--;
}//end for

//free memory
free(pT.a);
free(pT.b);
free(pT.t);
free(pT.theta);

free(sT.a);
free(sT.b);
free(sT.t);

//save results in erg
erg[0]=Hl[i][interval].hl;
erg[1]=Hl[i][interval].hu;


}//end main

