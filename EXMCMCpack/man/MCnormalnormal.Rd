\name{MCnormalnormal}
\alias{MCnormalnormal}
\title{Monte Carlo Simulation from a Normal Likelihood (with known variance) with a Normal Prior}

\description{
  This function generates a sample from the posterior distribution
  of a Normal likelihood (with known variance) with a Normal prior.
  }
  
\usage{
MCnormalnormal(y, sigma2, mu0, tau20, mc=1000, ...)
}

\arguments{
    \item{y}{The data.}

    \item{sigma2}{The known variance of y.}

    \item{mu0}{The prior mean of mu.}

    \item{tau20}{The prior variance of mu.}

    \item{mc}{The number of Monte Carlo draws to make.}
        
    \item{...}{further arguments to be passed}
    
}

\value{
   An mcmc object that contains the posterior sample.  This 
   object can be summarized by functions provided by the coda package.
}

\details{
  \code{MCnormalnormal} directly simulates from the posterior distribution. 
  This model is designed primarily for instructional use.  \eqn{\mu}{mu}
  is the parameter of interest of the Normal distribution.
  We assume
  a conjugate normal prior:
  \deqn{\mu \sim \mathcal{N}(\mu_0, \tau^2_0)}{mu ~ N(mu0, tau20)}
  \eqn{y} is a vector of observed data.
  }
  
\examples{
\dontrun{
y <- c(2.65, 1.80, 2.29, 2.11, 2.27, 2.61, 2.49, 0.96, 1.72, 2.40)
posterior <- MCMCpack:::MCnormalnormal(y, 1, 0, 1, 5000)
summary(posterior)
plot(posterior)
grid <- seq(-3,3,0.01)
plot(grid, dnorm(grid, 0, 1), type="l", col="red", lwd=3, ylim=c(0,1.4),
   xlab="mu", ylab="density")
lines(density(posterior), col="blue", lwd=3)
legend(-3, 1.4, c("prior", "posterior"), lwd=3, col=c("red", "blue"))
}
}

\keyword{models}

\seealso{\code{\link[coda]{plot.mcmc}},
  \code{\link[coda]{summary.mcmc}}}

