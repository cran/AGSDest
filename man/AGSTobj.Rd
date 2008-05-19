\name{AGSTobj}
\alias{AGSTobj}
\alias{print.AGSTobj}
\alias{plot.AGSTobj}
\alias{summary.AGSTobj}
\alias{print.summary.AGSTobj}
\title{Adaptive group sequential trial object (AGSTobj)}
\description{
  The \code{AGSTobj} includes design and outcome of primary and secondary trial.
}
\usage{
AGSTobj(x, ...)
\method{print}{AGSTobj}(x, ...)
\method{plot}{AGSTobj}(x,main=c("primary trial","secondary trial"),print.pdf=FALSE,...)
\method{summary}{AGSTobj}(object,ctype="b",ptype="b",etype="a",overwrite=FALSE,...)
\method{print}{summary.AGSTobj}(x, ...)
}
\arguments{
  \item{x}{ object of the \code{class} \code{AGSTobj}}
  \item{object}{ object of the \code{class} \code{AGSTobj}}
  \item{main}{ Title of the plots (default: first plot: "primary trial"; second plot: "secondary trial")}
  \item{print.pdf}{ option; if TRUE a pdf file is created. Instead of setting print.pdf to TRUE, the user can specify a character string giving the name or the path of the file.}
  \item{ctype}{ confidence type: repeated "r", stage-wise ordering "so", both "b" or none "n" (default: "b")}
  \item{ptype}{ p-value type: repeated "r", stage-wise ordering "so", both "b" or none "n" (default: "b")}
  \item{etype}{ point estimate: maximum likelihood "ml", median unbiased "mu", all "a" or none "n" (default: "a")}
  \item{overwrite}{ option; if TRUE all old values are deleted and new values are calculated (default: FALSE) }
  \item{\dots}{ additional arguments. }
}
\details{
  A \code{AGSTobj} object is designed.
  
  The function \code{summary} returns an object of \code{class} \code{AGSTobj}.
  
  \code{ctype} defines the type of confidence interval that is calculated.
  \tabular{ll}{  
    \code{"r"} \tab Repeated confidence bound for a GSD with design adaptations\cr
    \code{"so"} \tab Confidence bound for a GSD with design adaptation based on the stage-wise ordering\cr
    \code{"b"} \tab both: repeated confidence bound and confidence bound based on the stage-wise ordering for a GSD with design adaptations\cr
    \code{"n"} \tab no confidence bound is calculated\cr
  }
  The calculated confidence bounds are saved as:  
  \tabular{ll}{
  \code{cb.r} \tab repeated confidence bound \cr
  \code{cb.so} \tab confidence bound based on the stage-wise ordering \cr
  }
    
  \code{ptype} defines the type of p-value that is calculated.
  \tabular{ll}{
    \code{"r"}\tab Repeated p-value for a GSD with design adaptations \cr
    \code{"so"}\tab Stage-wise adjusted p-value for a GSD with design adaptations \cr
    \code{"b"}\tab both: repeated and stage-wise adjusted p-value for a GSD with design adaptations\cr
    \code{"n"} \tab no p-value is calculated\cr
  }
  The calculated p-values are saved as:  
  \tabular{ll}{
  \code{pvalue.r}\tab repeated p-value \cr
  \code{pvalue.so}\tab stage-wise adjusted p-value\cr
  }
  
  \code{etype} defines the type of point estimate
    \tabular{ll}{
    \code{"ml"}\tab maximum likelihood estimate (ignoring the sequential and adaptive nature of the design)\cr
    \code{"mu"}\tab median unbiased estimate (stage-wise lower confidence bound at level 0.5) for a GSD with design adaptations\cr
    \code{"cons"}\tab conservative estimate (repeated lower confidence bound at level 0.5) for a GSD with design adaptations\cr
    \code{"a"}\tab all: maximum likelihood, median unbiased and conservative point estimate for a GSD with design adaptations\cr
    \code{"n"}\tab No point estimate is calculated\cr                                                                                                             
    }
  
  The calculated point estimates are saved as:
  \tabular{ll}{
  \code{est.ml}\tab Maximum likelihood estimate\cr
  \code{est.mu}\tab Median unbiased estimate\cr
  \code{est.cons}\tab Conservative estimate\cr
  }
  
  The stage-wise adjusted confidence bound, p-value and the median unbiased point estimate can only be calculated at the stage where the trial stops and are only valid if the stopping rule is met.

  The repeated confidence bound, repeated p-value, conservative estimate and maximum likelihood estimate can be calculated at every stage of the trial and
  not just at the stage where the trial stops and are also valid if the stopping rule is not met. 
  For calculating the repeated confidence bounds or p-values the user has to specify \code{sTo} (secondary trial outcome) in the object \code{AGSTobj} (see example below). If the stopping rule is not met in object \code{sTo} then stage-wise adjusted confidence bounds and p-values will not be computed while a warning message is given when their computation have erroneously been specified. 
}
\value{
  An object of class \code{AGSTobj}, which is basically a list with the elements

\item{cb.so}{ confidence bound based on the stage-wise ordering (stage-wise adjusted confidence bound)}

\item{cb.r}{ repeated confidence bound}

\item{pvalue.so}{ p-value based on the stage-wise ordering (stage-wise adjusted p-value)}

\item{pvalue.r}{ repeated p-value}

\item{est.ml}{ maximum likelihood estimate}

\item{est.mu}{ median unbiased point estimate}

\item{est.cons}{ conservative point estimate}

\item{pT}{}
  
  \item{K}{ number of stages}

  \item{al}{ alpha (type I error rate)}

  \item{a}{ lower critical bounds of primary group sequential design (are currently always set to -8)}

  \item{b}{ upper critical bounds of primary group sequential design}

  \item{t}{ vector with cumulative information fraction}

  \item{SF}{ spending function (for details see below)}

  \item{phi}{ parameter of spending function when SF=3 or 4 (for details see below)}

  \item{alab}{ alpha-absorbing parameter values of primary group sequential design}

  \item{als}{ alpha-values ''spent'' at each stage of primary group sequential design}
  
  \item{Imax}{ maximum information number}
  
  \item{delta}{ effect size used for planning the primary trial}

  \item{cp}{ conditional power for planning the primary trial}

\item{iD}{}

  \item{L}{ stage of the adaptation}
  
  \item{z}{ z-statistic at adaptive interim analysis}

\item{sT}{}
  
  \item{K}{ number of stages }

  \item{al}{ conditional rejection probability}

  \item{a}{ lower critical bounds of secondary group sequential design (are currently always set to -8)}

  \item{b}{ upper critical bounds of secondary group sequential design}

  \item{t}{ vector with cumulative information fraction}

  \item{SF}{ spending function (for details see below)}

  \item{phi}{ parameter of spending function when SF=3 or 4 (for details see below)}
  
  \item{Imax}{ maximum information number}
  
  \item{delta}{ effect size used for planning the secondary trial}

  \item{cp}{ conditional power for planning the secondary trial}

\item{sTo}{}

  \item{T}{ stage where trial stops }
  
  \item{z}{ z-statistic at stage where trial stops}
 
}
\author{ Niklas Hack \email{niklas.hack@meduniwien.ac.at} and Werner Brannath \email{werner.brannath@meduniwien.ac.at} }
\note{ 
 The \code{AGSTobj} should always have the same ordering and names as given in the list above or as given in the example. 
 
 1. \code{pt}, 2. \code{iD}, 3. \code{sT}, 4. \code{sTo}
  
 \code{SF} defines the spending function.
   
    \code{SF} = 1 O'Brien and Fleming type spending function of Lan and DeMets (1983)
    
    \code{SF} = 2 Pocock type spending function of Lan and DeMets (1983)
    
    \code{SF} = 3 Power family (\eqn{c_\alpha* t^\phi}); phi must be greater than 0.
    
    \code{SF} = 4 Hwang-Shih-DeCani family;\eqn{(1-e^{-\phi t})/(1-e^{-\phi})}, where \code{phi} cannot be 0.

 A value of \code{SF}=3 corresponds to the power family. Here, the spending function is \eqn{t^{\phi}}, 
 where phi must be greater than 0. A value of \code{SF}=4 corresponds to the Hwang-Shih-DeCani family, 
 with the spending function \eqn{(1-e^{-\phi t})/(1-e^{-\phi})}, where phi cannot be 0.
 
 If a path is specified for \code{print.pdf}, all \\ must be changed to /. If a filename is specified the ending of the file must be (.pdf).

 In the current version the vector of lower bounds \code{a} should be set to rep(-8,K)
}
\seealso{ \code{\link{AGSTobj}}, \code{\link{print.AGSTobj}}, \code{\link{plot.AGSTobj}}, \code{\link{summary.AGSTobj}}}
\examples{

pT=plan.GST(K=3,SF=4,phi=-4,alpha=0.05,delta=6,pow=0.9,compute.alab=TRUE,compute.als=TRUE)

iD=list(T=1, z=1.090728)

swImax=0.0625

I2min=3*swImax
I2max=3*swImax

sT=adapt(pT=pT,iD=iD,SF=1,phi=0,cp=0.8,theta=5,I2min,I2max,swImax)

sTo=list(T=2, z=2.393)

AGST<-as.AGST(pT=pT,iD=iD,sT=sT,sTo=sTo)

AGST
plot(AGST)

AGST<-summary(AGST)
plot(AGST)

##The repeated confidence interval and p-value at an earlier stage than the one where the trial stops (T=3).

summary(as.AGST(pT,iD,sT,sTo=list(T=1,z=1.7)),ctype="r",ptype="r")

\dontrun{
##If the stage-wise adjusted confidence interval is calculated at this stage, 
##the function returns an error message

summary(as.AGST(pT,iD,sT,sTo=list(T=1,z=1.7)),ctype="so",ptype="so")
}
}
\keyword{datasets}
