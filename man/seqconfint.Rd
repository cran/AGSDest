\name{seqconfint}
\alias{seqconfint}
\title{ Calculates confidence interval }
\description{
  Calculates the repeated confidence bound or the confidence bound based on the stage-wise ordering of a GSD or a AGSD
}
\usage{
seqconfint(object, type = "b", level = NULL)
}
\arguments{
  \item{object}{ object of the \code{class} \code{GSTobj} or of the \code{class} \code{AGSTobj}}
  \item{type}{ confidence type: repeated "r", stage-wise ordering "so" or both "b" (default: "b")}
  \item{level}{ type I error rate (default: NULL) }
}
\details{
  \code{object} can be an object of the \code{class} \code{GSTobj} or an object of the \code{class} \code{AGSTobj}. 
  The function identifies the \code{class} of the object and calculates the corresponding confidence interval (classical or adaptive).
  
  If \code{object} has \code{class} \code{GSTobj}, then a confidence bound for a classical GSD is calculated.
  \code{type} defines the type of confidence interval that is calculated
  \tabular{ll}{
    \code{"r"}\tab Repeated confidence bound for a classical GSD\cr
    \code{"so"}\tab Confidence bound for a classical GSD based on the stage-wise ordering\cr
    \code{"b"}\tab both: repeated confidence bound and confidence bound based on the stage-wise for a GSD\cr
  }
  If \code{object} has \code{class} \code{AGSTobj}, then a confidence bound for a GSD with design adaptation is calculated.
  \code{type} defines the type of confidence interval that is calculated
  \tabular{ll}{
    \code{"r"}\tab Repeated confidence bound for a GSD with design adaptations\cr
    \code{"so"}\tab Confidence bound for a GSD with design adaptation based on the stage-wise ordering\cr
    \code{"b"}\tab both: repeated confidence bound and confidence bound based on the stage-wise for a GSD with design adaptations\cr
  }
  By setting \code{level} to the value 0.5 the conservative point estimate is calculated. Default is the \code{level} of the primary trial.

}
\value{
  The function \code{seqconfint} returns according to the \code{class} of \code{object} the classical or adaptive confidence bound.
  If \code{object} has \code{class} \code{GSTobj} the classical confidence bound is calculated. If the 
  parameter value has the \code{class} \code{AGSTobj} the adaptive confidence bound is calculated.
  
  The calculated confidence bounds are saved as:
  
  \item{cb.r}{ repeated confidence bound}
  \item{cb.so}{ confidence bound based on the stage-wise ordering}

  If the \code{level} is set to 0.5, the calculated point estimates are:

  \item{est.mu}{ Median unbiased point estimate, based on the stage-wise ordering}
  \item{est.cons}{Flexible, but conservative repeated point estimate}

}
\note{
The stage-wise adjusted confidence interval can only be calculated at the stage where the trial stops and is only valid if the stopping rule is met.

The repeated confidence interval can be calculated at every stage of the trial and
not just at the stage where the trial stops and is also valid if the stopping rule is not met.

For calculating the sequential confidence intervals at stage \code{T} the user has to specify the outcome \code{GSDo} in the object \code{GSTobj}
or \code{sTo} (secondary trial outcome) in the object \code{AGSTobj}. A trial outcome is a list of the form 
\code{list=(T=stage of interim analysis, z = interim z-statistic)}; see the example below.
}
\references{ 
Brannath, W, Mehta, CR, Posch, M (2008) ''Exact confidence bounds following
adaptive group sequential tests'', \emph{Biometrics} accepted.

Jennison, C, Turnbull, BW (1989) ''Repeated confidence intervals for group
sequential clinical trials'', \emph{Contr. Clin. Trials}, 5, 33-45.

Mehta, CR, Bauer, P, Posch, M, Brannath, W (2007) ''Repeated confidence
intervals for adaptive group sequential trials'', \emph{Statistics in Medicine}, 26, 5422-5433.

Mueller, HH, Schaefer, H (2001) ''Adaptive group sequential design for clinical
trials: Combining the advantages of adaptive and of classical group sequential
approaches'', \emph{Biometrics}, 57, 886-891. 

Tsiatis,AA, Rosner,GL, Mehta,CR (1984) ''Exact confidence intervals
following a group sequential test'', \emph{Biometrics}, 40, 797-804.
}
\author{ Niklas Hack \email{niklas.hack@meduniwien.ac.at} and Werner Brannath \email{werner.brannath@meduniwien.ac.at} }
\seealso{\code{\link{AGSTobj}}, \code{\link{GSTobj}} }
\examples{
##The following calculates the repeated confidence bound of a group sequential trial

GSD=plan.GST(K=4,SF=1,phi=0,alpha=0.025,delta=6,pow=0.8,compute.alab=TRUE,compute.als=TRUE)

GST<-as.GST(GSD=GSD,GSDo=list(T=2, z=3.1))
seqconfint(GST,type="r")

##The confidence bound based on the stage-wise ordering of a group sequential trial is calculated by

seqconfint(GST,type="so")

##The repeated confidence interval at the earlier stage T=1 where the
##trial stopping rule is not met.

seqconfint(as.GST(GSD,GSDo=list(T=1,z=0.7)),type="r")

##The repeated confidence bound and the confidence bound 
##based on the stage-wise ordering of a group sequential trial 
##after a design adaptation is calculated by

pT=plan.GST(K=3,SF=4,phi=-4,alpha=0.05,delta=6,pow=0.9,compute.alab=TRUE,compute.als=TRUE)

iD=list(T=1, z=1.090728)

swImax=0.0625

I2min=3*swImax
I2max=3*swImax

sT=adapt(pT=pT,iD=iD,SF=1,phi=0,cp=0.8,theta=5,I2min,I2max,swImax)

sTo=list(T=2, z=2.393)

AGST<-as.AGST(pT=pT,iD=iD,sT=sT,sTo=sTo)
seqconfint(AGST)

##The repeated confidence interval at the earlier stage T=2 where the
##trial stopping rule is not met.

seqconfint(as.AGST(pT,iD,sT,sTo=list(T=2,z=1.7)),type="r")

\dontrun{
##If the stage-wise adjusted confidence interval is calculated at this stage, 
##the function returns an error message

seqconfint(as.AGST(pT,iD,sT,sTo=list(T=2,z=1.7)),type="so")
}
}
\keyword{methods}
