\name{plotState}
\alias{plotState}
\title{Changepoint State Plot}
\description{Plot the posterior probability that each time point is in each state.}

\usage{
   plotState(mcmcout, main="Posterior Regime Probability", ylab=expression(paste("Pr(", S[t], "= k |", Y[t], ")")),
   legend.control = NULL, cex = 0.8, lwd = 1.2, start=1)
}

\arguments{

\item{mcmcout}{The \code{mcmc} object containing the posterior density sample from a changepoint model.  Note that this must have a \code{prob.state} attribute.}

\item{main}{Title of the plot.} 

\item{ylab}{Label for the y-axis.}

\item{legend.control}{Control the location of the legend.  It is necessary to pass both the x and y locations; i.e., \code{c(x,y)}.}

\item{cex}{Control point size.}

\item{lwd}{Line width parameter.}

\item{start}{The time of the first observation to be shown in the time series plot.}

}

\author{Jong Hee Park, \email{jhp@uchicago.edu}, \url{http://home.uchicago.edu/~jhp/}.}

\keyword{hplot}

\seealso{\code{\link{MCMCpoissonChange}}, \code{\link{MCMCbinaryChange}}}
