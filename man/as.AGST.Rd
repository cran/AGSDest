\name{as.AGST}
\alias{as.AGST}
\title{ as Adaptive Group Sequential Trial }
\description{
  Function \code{as.AGST} builds an adaptive group sequential trial object 
}
\usage{
as.AGST(pT,iD,sT,sTo)
}
\arguments{
  \item{pT}{ object of the \code{class} \code{GSTobj}; primary trial design }
  \item{iD}{ interim data; a list with the variables \code{T} and \code{z}; list(T = stage of interim analysis, z = interim z-statistic)}
  \item{sT}{ object of the \code{class} \code{GSTobj}; secondary trial design}
  \item{sTo}{ secondary trial outcome; a list with the variables \code{T} and \code{z}; list(T = stage where trial stops, z = z-statistic at stage where trial stops)}
}
\value{
  \code{as.AGST} returns a list containing the \code{pT}, \code{iD}, \code{sT} and \code{sTo} with \code{class}=\code{AGSTobj}
}
\author{ Niklas Hack \email{niklas.hack@meduniwien.ac.at} and Werner Brannath \email{werner.brannath@meduniwien.ac.at} }
\seealso{ \code{\link{AGSTobj}} }
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
}
\keyword{methods}