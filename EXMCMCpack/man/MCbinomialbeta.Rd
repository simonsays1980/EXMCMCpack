\name{MCbinomialbeta}
\alias{MCbinomialbeta}
\title{Monte Carlo Simulation from a Binomial Likelihood with a Beta Prior}

\description{
  This function generates a sample from the posterior distribution
  of a binomial likelihood with a Beta prior.
  }
  
\usage{
MCbinomialbeta(y, n, alpha=1, beta=1, mc=1000, ...)
}

\arguments{
    \item{y}{The number of successes in the independent Bernoulli trials.}

    \item{n}{The number of independent Bernoulli trials.}

    \item{alpha}{Beta prior distribution alpha parameter.}

    \item{beta}{Beta prior distribution beta parameter.}

    \item{mc}{The number of Monte Carlo draws to make.}
        
    \item{...}{further arguments to be passed}       
}

\value{
   An mcmc object that contains the posterior sample.  This 
   object can be summarized by functions provided by the coda package.
}

\details{
  \code{MCbinomialbeta} directly simulates from the posterior distribution. 
  This model is designed primarily for instructional use.  \eqn{\pi}{pi} is
  the probability of success for each independent Bernoulli trial.  We assume
  a conjugate Beta prior:
  \deqn{\pi \sim \mathcal{B}eta(\alpha, \beta)}{pi ~ Beta(alpha, beta)}
  \eqn{y} is the number of successes in \eqn{n} trials.
  By default, a uniform prior is used.
  }
  
\examples{
\dontrun{
posterior <- MCbinomialbeta(3,12,mc=5000)
summary(posterior)
plot(posterior)
grid <- seq(0,1,0.01)
plot(grid, dbeta(grid, 1, 1), type="l", col="red", lwd=3, ylim=c(0,3.6),
  xlab="pi", ylab="density")
lines(density(posterior), col="blue", lwd=3)
legend(.75, 3.6, c("prior", "posterior"), lwd=3, col=c("red", "blue"))
}
}

\keyword{models}

\seealso{\code{\link[coda]{plot.mcmc}},
  \code{\link[coda]{summary.mcmc}}}

