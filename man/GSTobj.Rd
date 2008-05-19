\name{GSTobj}
\alias{GSTobj}
\alias{print.GSTobj}
\alias{plot.GSTobj}
\alias{summary.GSTobj}
\alias{print.summary.GSTobj}
\title{Group sequential trial object (GSTobj)}
\description{
  The \code{GSTobj} includes design and outcome of primary trial.
}
\usage{
GSTobj(x, ...)
\method{print}{GSTobj}(x, ...)
\method{plot}{GSTobj}(x,main="GSD",print.pdf=FALSE, ...)
\method{summary}{GSTobj}(object,ctype="b",ptype="b",etype="a",overwrite=FALSE,...)
\method{print}{summary.GSTobj}(x, ...)
}
\arguments{
  \item{x}{ object of the \code{class} \code{GSTobj} }
  \item{object}{ object of the \code{class} \code{GSTobj}}
  \item{main}{ Title of the plots (default: "GSD")}
  \item{print.pdf}{ option; if TRUE a pdf file is created. Instead of setting print.pdf to TRUE, the user can specify a character string giving the name or the path of the file.}
  \item{ctype}{ confidence type: repeated "r", stage-wise ordering "so", both "b" or none "n" (default: "b")}
  \item{ptype}{ p-value type: repeated "r", stage-wise ordering "so", both "b" or none "n" (default: "b")}
  \item{etype}{ point estimate: maximum likelihood "ml", median unbiased "mu", all "a" or none "n" (default: "a")}
  \item{overwrite}{ option; if TRUE all old values are deleted and new values are calculated (default: FALSE) }
  \item{\dots}{ additional arguments. }
}
\details{
  A \code{GSTobj} object is designed.
  
  The function \code{summary} returns an object of \code{class} \code{GSTobj}.
  
  \code{ctype} defines the type of confidence interval that is calculated.
  \tabular{ll}{
    \code{"r"}\tab Repeated confidence bound for a classical GSD\cr
    \code{"so"}\tab Confidence bound for a classical GSD based on the stage-wise ordering\cr
    \code{"b"}\tab both: repeated confidence bound and confidence bound based on the stage-wise for a classical GSD\cr
    \code{"n"}\tab no confidence bound is calculated\cr
  }
  The calculated confidence bounds are saved as:
  \tabular{ll}{
  \code{cb.r}\tab repeated confidence bound\cr
  \code{cb.so}\tab confidence bound based on the stage-wise ordering\cr
  }
  
  \code{ptype} defines the type of p-value that is calculated.
    \tabular{ll}{
    \code{"r"}\tab Repeated p-value for a classical GSD\cr
    \code{"so"}\tab Stage-wise adjusted p-value for a classical GSD \cr
    \code{"b"}\tab both: repeated and stage-wise adjusted p-value for a classical GSD\cr
    \code{"n"}\tab no p-value is calculated\cr
  }
  
  The calculated p-values are saved as:
  \tabular{ll}{
  \code{pvalue.r}\tab repeated p-value\cr
  \code{pvalue.so}\tab stage-wise adjusted p-value\cr
  }
  
  \code{etype} defines the type of point estimate
    \tabular{ll}{
    \code{"ml"}\tab maximum likelihood estimate (ignoring the sequential nature of the design)\cr
    \code{"mu"}\tab median unbiased estimate (stage-wise lower confidence bound at level 0.5) for a classical GSD\cr
    \code{"cons"}\tab Conservative estimate (repeated lower confidence bound at level 0.5) for a classical GSD\cr
    \code{"a"}\tab all: maximum likelihood, median unbiased and conservative point estimate for a classical GSD\cr
    \code{"n"}\tab No point estimate is calculated\cr                                                                                                             
    }
  
  The calculated point estimates are saved as:
  \tabular{ll}{
  \code{est.ml}\tab Maximum likelihood estimate\cr
  \code{est.mu}\tab Median unbiased estimate\cr
  \code{est.cons}\tab Conservative estimate\cr
  }
  
  The stage-wise adjusted confidence interval and p-value and the median unbiased point estimate can only be calculated at the stage where the trial stops and is only valid if the stopping rule is met.

  The repeated confidence interval and repeated p-value, conservative estimate and maximum likelihood estimate can be calculated at every stage of the trial and
  not just at the stage where the trial stops and is also valid if the stopping rule is not met. 
  For calculating the repeated confidence interval or p-value at any stage of the trial the user has to specify the outcome \code{GSDo} in the object \code{GSTobj} (see example below).
}
\value{
  An object of class \code{GSTobj}, is basically a list with the elements

\item{cb.so}{ confidence bound based on the stage-wise ordering}

\item{cb.r}{ repeated confidence bound}

\item{pvalue.so}{ stage-wise adjusted p-value}

\item{pvalue.r}{ repeated p-value}

\item{est.ml}{ maximum likelihood estimate}

\item{est.mu}{ median unbiased point estimate}

\item{est.cons}{ conservative point estimate}

\item{GSD}{}

  \item{K}{ number of stages}

  \item{al}{ alpha (type I error rate)}

  \item{a}{ lower critical bounds of group sequential design (are currently always set to -8)}

  \item{b}{ upper critical bounds of group sequential design }

  \item{t}{ vector with cumulative information fraction}

  \item{SF}{ spending function (for details see below)}

  \item{phi}{ parameter of spending function when SF=3 or 4 (for details see below)}

  \item{alab}{ alpha-absorbing parameter values of group sequential design}

  \item{als}{ alpha-values ''spent'' at each stage of group sequential design}
  
  \item{Imax}{ maximum information number}
  
  \item{delta}{ effect size used for planning the primary trial}

  \item{cp}{ conditionla power of the trial}
  
\item{GSDo}{}

  \item{T}{ stage where trial stops}
  
  \item{z}{ z-statistic at stage where trial stops}

}
\author{ Niklas Hack \email{niklas.hack@meduniwien.ac.at} and Werner Brannath \email{werner.brannath@meduniwien.ac.at} }
\note{ 
\code{SF} defines the spending function.

    \code{SF} = 1 O'Brien and Fleming type spending function of Lan and DeMets (1983)

    \code{SF} = 2 Pocock type spending function of Lan and DeMets (1983)

    \code{SF} = 3 Power family (\eqn{c_\alpha* t^\phi}). phi must be greater than 0.

    \code{SF} = 4 Hwang-Shih-DeCani family.\eqn{(1-e^{-\phi t})/(1-e^{-\phi})}, where \code{phi} cannot be 0.

A value of \code{SF}=3 corresponds to the power family. Here, the spending function is \eqn{t^{\phi}},
where phi must be greater than 0. A value of \code{SF}=4 corresponds to the Hwang-Shih-DeCani family,
with the spending function \eqn{(1-e^{-\phi t})/(1-e^{-\phi})}, where phi cannot be 0.

If a path is specified for \code{print.pdf}, all \ must be changed to /. If a filename is specified the ending of the file must be (.pdf).

In the current version \code{a} should be set to rep(-8,K)
}
\seealso{ \code{\link{GSTobj}}, \code{\link{print.GSTobj}}, \code{\link{plot.GSTobj}}, \code{\link{summary.GSTobj}}}
\examples{

GSD=plan.GST(K=4,SF=1,phi=0,alpha=0.025,delta=6,pow=0.8,compute.alab=TRUE,compute.als=TRUE)

GST<-as.GST(GSD=GSD,GSDo=list(T=2, z=3.1))

GST
plot(GST)

GST<-summary(GST)
plot(GST)

##The repeated confidence interval, p-value and maximum likelihood estimate 
##at the earlier stage T=1 where the trial stopping rule is not met.

summary(as.GST(GSD,GSDo=list(T=1,z=0.7)),ctype="r",ptype="r",etype="ml")

\dontrun{
##If e.g. the stage-wise adjusted confidence interval is calculated at this stage, 
##the function returns an error message

summary(as.GST(GSD,GSDo=list(T=1,z=0.7)),ctype="so",etype="mu")
}
}
\keyword{datasets}