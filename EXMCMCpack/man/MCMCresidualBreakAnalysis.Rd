\name{MCMCresidualBreakAnalysis}
\alias{MCMCresidualBreakAnalysis}
\title{Break Analysis of Univariate Time Series using Markov Chain Monte Carlo}
\description{
	This function performs a break analysis for univariate time series data using 
	a linear Gaussian changepoint model. The code is written mainly for an internal 
	use in \code{testpanelSubjectBreak}.}
  
\usage{
MCMCresidualBreakAnalysis(resid, m = 1, 
           b0 = 0, B0 = 0.001, c0 = 0.1, d0 = 0.1, a = NULL, b = NULL,
           mcmc = 1000, burnin = 1000,  thin = 1, verbose = 0, 
           seed = NA, beta.start = NA, P.start = NA,
           marginal.likelihood = c("none", "Chib95"), ...) }

\arguments{
    \item{resid}{Univariate time series}

    \item{m}{The number of breaks.}

    \item{b0}{The prior mean of \eqn{\beta}{beta}.  This can either be a 
    scalar or a
    column vector with dimension equal to the number of betas. If this
    takes a scalar  value, then that value will serve as the prior
    mean for all of the betas.} 
    
    \item{B0}{The prior precision of \eqn{\beta}{beta}.  This can either be a 
    scalar or a square matrix with dimensions equal to the number of betas.  
    If this
    takes a scalar value, then that value times an identity matrix serves
    as the prior precision of beta. Default value of 0 is equivalent to
    an improper uniform prior for beta.}
    
    \item{c0}{\eqn{c_0/2}{c0/2} is the shape parameter for the inverse
    Gamma prior on \eqn{\sigma^2}{sigma^2} (the variance of the
    disturbances). The amount of information in the inverse Gamma prior
    is something like that from \eqn{c_0}{c0} pseudo-observations.} 
    
    \item{d0}{\eqn{d_0/2}{d0/2} is the scale parameter for the
    inverse Gamma prior on \eqn{\sigma^2}{sigma^2} (the variance of the
    disturbances). In constructing the inverse Gamma prior,
    \eqn{d_0}{d0} acts like the sum of squared errors from the
    \eqn{c_0}{c0} pseudo-observations.}

     \item{a}{\eqn{a}{a} is the shape1 beta prior for transition probabilities. By default,
    the expected duration is computed and corresponding a and b values are assigned. The expected
    duration is the sample period divided by the number of states.}

    \item{b}{\eqn{b}{b} is the shape2 beta prior for transition probabilities. By default,
    the expected duration is computed and corresponding a and b values are assigned. The expected
    duration is the sample period divided by the number of states.}

    \item{burnin}{The number of burn-in iterations for the sampler.}

    \item{mcmc}{The number of MCMC iterations after burnin.}

    \item{thin}{The thinning interval used in the simulation.  The number of
      MCMC iterations must be divisible by this value.}

    \item{verbose}{A switch which determines whether or not the progress of
    the sampler is printed to the screen.  If \code{verbose} is greater
    than 0 the iteration number, the
    \eqn{\beta}{beta} vector, and the error variance are printed to 
    the screen every \code{verbose}th iteration.}

    \item{seed}{The seed for the random number generator.  If NA, the Mersenne
    Twister generator is used with default seed 12345; if an integer is 
    passed it is used to seed the Mersenne twister.  The user can also
    pass a list of length two to use the L'Ecuyer random number generator,
    which is suitable for parallel computation.  The first element of the
    list is the L'Ecuyer seed, which is a vector of length six or NA (if NA 
    a default seed of \code{rep(12345,6)} is used).  The second element of 
    list is a positive substream number. See the MCMCpack 
    specification for more details.}

    \item{beta.start}{The starting values for the \eqn{\beta}{beta} vector.
    This can either be a scalar or a
    column vector with dimension equal to the number of betas.
    The default value of of NA will use the OLS
    estimate of \eqn{\beta}{beta} as the starting value.  If this is a 
    scalar,  that value will serve as the starting value
    mean for all of the betas.}

    \item{P.start}{The starting values for the transition matrix. 
    A user should provide a square matrix with dimension equal to the number of states. 
    By default, draws from the \code{Beta(0.9, 0.1)} 
    are used to construct a proper transition matrix for each raw except the last raw.}	

    \item{marginal.likelihood}{How should the marginal likelihood be
    calculated? Options are: \code{none} in which case the marginal
    likelihood will not be calculated, and
    \code{Chib95} in which case the method of Chib (1995) is used.}
  
    \item{...}{further arguments to be passed}       
}

\value{
   An mcmc object that contains the posterior sample.  This 
   object can be summarized by functions provided by the coda package.
}

\details{
  \code{MCMCresidualBreakAnalysis} simulates from the posterior distribution using 
  standard Gibbs sampling (a multivariate Normal draw for the betas, and an
  inverse Gamma draw for the conditional error variance).  The simulation
  proper is done in compiled C++ code to maximize efficiency.  Please consult
  the coda documentation for a comprehensive list of functions that can be
  used to analyze the posterior sample.
  
  The model takes the following form:
  \deqn{y_{i} \sim \mathcal{N}(\beta_{m}, \sigma^2_{m}) \;\; m = 1, \ldots, M}{y_i ~ N(beta_m, sigma^2_m), m = 1,...,M.}
 
   We assume standard, semi-conjugate priors:
  \deqn{\beta \sim \mathcal{N}(b_0,B_0^{-1})}{beta ~ N(b0,B0^(-1))}
  And:
  \deqn{\sigma^{-2} \sim \mathcal{G}amma(c_0/2, d_0/2)}{sigma^(-2) ~
    Gamma(c0/2, d0/2)}
  Where \eqn{\beta}{beta} and \eqn{\sigma^{-2}}{sigma^(-2)} are assumed 
  \emph{a priori} independent.  
  
  And:
  \deqn{p_{mm} \sim \mathcal{B}eta(a, b),\;\; m = 1, \ldots, M}{p_mm ~ Beta(a, b), m = 1,...,M.}
  Where \eqn{M}{M} is the number of states.
  }
  
  \references{
   Jong Hee Park, 2011. ``A Unified Method for Dynamic and Cross-Sectional Heterogeneity: 
   Introducing Hidden Markov Panel Models." Working Paper.
   
   Siddhartha Chib. 1998. ``Estimation and comparison of multiple change-point models.''
   \emph{Journal of Econometrics}. 86: 221-241.
}


\examples{
\dontrun{
line   <- list(X = c(-2,-1,0,1,2), Y = c(1,3,3,3,5))
ols <- lm(Y~X)
residual <-   rstandard(ols)
posterior  <- MCMCresidualBreakAnalysis(residual, m = 1, data=line, mcmc=1000, verbose=200)
plotState(posterior)
summary(posterior)
}
}

\keyword{models}

\seealso{\code{\link[coda]{plot.mcmc}},
  \code{\link[coda]{summary.mcmc}}, \code{\link[stats]{lm}}}
