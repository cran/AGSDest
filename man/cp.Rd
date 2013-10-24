\name{cp}
\alias{cp}
\title{ coditional power of a GSD }
\description{
  \code{cp} is a function that computes the conditional power of a GSD.
}
\usage{
cp(GSD)
}
\arguments{
  \item{GSD}{ object of the \code{class} \code{GSTobj} or list with the following elements: K = number of stages, 
  a = vector with futility boundaries (not supported yet), b = rejection boundaries, t = vector with information fractions, Imax = maximum information number,
  delta = effect size used for planning the trial; see example blow.}
}
\value{
  \code{cp} returns the conditional power of a GSD.  
}
\references{ 
 O'Brien, PC, Fleming, TR (1979) ''A multiple testing procedure for clinical trials'',
\emph{Biometrics}, 35 , 549-556 

Schoenfeld, D (2001) ''A simple Algorithm for Designing Group Sequential Clinical Trials'',
\emph{Biometrics}, 27, 972-974 
}
\author{ Niklas Hack \email{niklas.hack@meduniwien.ac.at} and Werner Brannath \email{werner.brannath@meduniwien.ac.at} }
\seealso{ \code{\link{GSTobj}}}
\examples{
##The following calculates the conditional power of a GSD. 

GSD<-list(K=4,a=rep(-8,4),b=c(4.333,2.963,2.359,2.014),
t=c(0.25,0.5,0.75,1),Imax=0.22,delta=4)

cp(GSD)
}
\keyword{methods}
