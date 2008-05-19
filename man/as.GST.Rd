\name{as.GST}
\alias{as.GST}
\title{ as Group Sequential Trial }
\description{
  Function \code{as.GST} builds a group sequential trial object 
}
\usage{
as.GST(GSD,GSDo)
}
\arguments{
  \item{GSD}{ object of the \code{class} \code{GSTobj}; group sequential design}
  \item{GSDo}{ group sequential design outcome;  a list with the variables \code{T} and \code{z}; list(T = stage where trial stops, z = z-statistic at stage where trial stops)}
}
\value{
  \code{as.GST} returns a list containing the \code{GSD} and \code{GSDo} with \code{class}=\code{GSTobj}
}
\author{ Niklas Hack \email{niklas.hack@meduniwien.ac.at} and Werner Brannath \email{werner.brannath@meduniwien.ac.at} }
\seealso{ \code{\link{GSTobj}} }
\examples{

GSD=plan.GST(K=4,SF=1,phi=0,alpha=0.025,delta=6,pow=0.8,compute.alab=TRUE,compute.als=TRUE)

GSDo=list(T=2, z=3.1)

GST=as.GST(GSD=GSD,GSDo=GSDo)

GST
}
\keyword{methods}