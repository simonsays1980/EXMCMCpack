\name{plotChangepoint}
\alias{plotChangepoint}
\title{Posterior Density of Regime Change Plot}
\description{Plot the posterior density of regime change.}

\usage{
   plotChangepoint(mcmcout, main="Posterior Density of Regime Change Probabilities", 
   xlab = "Time", ylab = "", verbose = FALSE, start=1, overlay=FALSE)
}

\arguments{

\item{mcmcout}{The \code{mcmc} object containing the posterior density sample from a changepoint model.  Note that this must have a \code{prob.state} attribute.}

\item{main}{Title of the plot} 

\item{xlab}{Label for the x-axis.}

\item{ylab}{Label for the y-axis.}

\item{verbose}{If verbose=TRUE, expected changepoints are printed.}

\item{start}{The time of the first observation to be shown in the time series plot.}

\item{overlay}{If overlay=TRUE, the probability of each regime change is drawn separately, which will be useful to draw multiple plots in one screen. See the example in \code{MCMCpoissonChange}. Otherwise, multiple plots of regime change probabilities will be drawn.}
}

\author{Jong Hee Park, \email{jhp@uchicago.edu}, \url{http://home.uchicago.edu/~jhp/}.}

\keyword{hplot}

\seealso{\code{\link{MCMCpoissonChange}}, \code{\link{MCMCbinaryChange}}}
