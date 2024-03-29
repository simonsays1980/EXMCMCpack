\name{MCMCprobitChange}
\alias{MCMCprobitChange}

\title{Markov Chain Monte Carlo for a linear Gaussian Multiple Changepoint Model}
\description{
  This function generates a sample from the posterior distribution
  of a linear Gaussian model with multiple changepoints. The function uses
  the Markov chain Monte Carlo method of Chib (1998).
  The user supplies data and priors, and
  a sample from the posterior distribution is returned as an mcmc
  object, which can be subsequently analyzed with functions
  provided in the coda package.
}

\usage{MCMCprobitChange(formula, data=parent.frame(),  m = 1,
        burnin = 10000, mcmc = 10000, thin = 1, verbose = 0,
        seed = NA, beta.start = NA, P.start = NA,
        b0 = NULL, B0 = NULL, a = NULL, b = NULL,
        marginal.likelihood = c("none", "Chib95"), ...)}

\arguments{
   \item{formula}{Model formula.}

    \item{data}{Data frame.}

    \item{m}{The number of changepoints.}

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
    The default value of of NA will use the MLE
    estimate of \eqn{\beta}{beta} as the starting value.  If this is a
    scalar,  that value will serve as the starting value
    mean for all of the betas.}

    \item{P.start}{The starting values for the transition matrix. 
    A user should provide a square matrix with dimension equal to the number of states. 
    By default, draws from the \code{Beta(0.9, 0.1)} 
    are used to construct a proper transition matrix for each raw except the last raw.}	

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

     \item{a}{\eqn{a}{a} is the shape1 beta prior for transition probabilities. By default,
    the expected duration is computed and corresponding a and b values are assigned. The expected
    duration is the sample period divided by the number of states.}

    \item{b}{\eqn{b}{b} is the shape2 beta prior for transition probabilities. By default,
    the expected duration is computed and corresponding a and b values are assigned. The expected
    duration is the sample period divided by the number of states.}

    \item{marginal.likelihood}{How should the marginal likelihood be
    calculated? Options are: \code{none} in which case the marginal
    likelihood will not be calculated, and
    \code{Chib95} in which case the method of Chib (1995) is used.}

    \item{...}{further arguments to be passed}
}

\value{
   An mcmc object that contains the posterior sample.  This
   object can be summarized by functions provided by the coda package.
   The object contains an attribute \code{prob.state} storage matrix that contains the probability of \eqn{state_i}{state_i} for each period, the
   log-likelihood of the model (\code{loglike}), and
   the log-marginal likelihood of the model (\code{logmarglike}).
}

\details{
  \code{MCMCprobitChange} simulates from the posterior distribution of
  a probit regression model with multiple parameter breaks. The simulation is based on Chib (1998) and Park (2011).

  The model takes the following form:
  \deqn{\Pr(y_t = 1) = \Phi(x_i'\beta_m) \;\; m = 1, \ldots, M}{
   Pr(y_t = 1) = Phi(x_i'beta_m)}
  Where \eqn{M}{M} is the number of states, and \eqn{\beta_m}{beta_m} 
  is a parameter when a state is \eqn{m}{m} at \eqn{t}{t}. 

  We assume Gaussian distribution for prior of \eqn{\beta}{beta}:
  \deqn{\beta_m \sim \mathcal{N}(b_0,B_0^{-1}),\;\; m = 1, \ldots, M}{beta_m ~ N(b0,B0^(-1)), m = 1,...,M.}

  And:
  \deqn{p_{mm} \sim \mathcal{B}eta(a, b),\;\; m = 1, \ldots, M}{p_mm ~ Beta(a, b), m = 1,...,M.}
  Where \eqn{M}{M} is the number of states.


}

\author{Jong Hee Park, \email{jhp@uchicago.edu}, \url{http://home.uchicago.edu/~jhp/}.}

\references{
  Jong Hee Park. 2011. ``Changepoint Analysis of Binary and Ordinal Probit Models:
  An Application to Bank Rate Policy Under the Interwar Gold Standard." \emph{Political Analysis}. 19: 188-204.

Andrew D. Martin, Kevin M. Quinn, and Jong Hee Park. 2011.
 ``MCMCpack: Markov Chain Monte Carlo in R.'',
 \emph{Journal of Statistical Software}. 42(9): 1-21.
 \url{http://www.jstatsoft.org/v42/i09/}.

  Siddhartha Chib. 1998. ``Estimation and comparison of multiple change-point models.''
   \emph{Journal of Econometrics}. 86: 221-241.

  Albert, J. H. and S. Chib. 1993. ``Bayesian Analysis of Binary and
  Polychotomous Response Data.'' \emph{J. Amer. Statist. Assoc.} 88, 669-679
}

\examples{
\dontrun{
set.seed(1973)
x1 <- rnorm(300, 0, 1)
true.beta <- c(-.5, .2, 1)
true.alpha <- c(.1, -1., .2)
X <- cbind(1, x1)

## set two true breaks at 100 and 200
true.phi1 <- pnorm(true.alpha[1] + x1[1:100]*true.beta[1])
true.phi2 <- pnorm(true.alpha[2] + x1[101:200]*true.beta[2])
true.phi3 <-  pnorm(true.alpha[3] + x1[201:300]*true.beta[3])

## generate y
y1 <- rbinom(100, 1, true.phi1)
y2 <- rbinom(100, 1, true.phi2)
y3 <- rbinom(100, 1, true.phi3)
Y <- as.ts(c(y1, y2, y3))

## fit multiple models with a varying number of breaks
out0 <- MCMCprobitChange(formula=Y~X-1, data=parent.frame(), m=0,
                         mcmc=1000, burnin=1000, thin=1, verbose=1000, 
                         b0 = 0, B0 = 10, a = 1, b = 1,  marginal.likelihood = c("Chib95"))	
out1 <- MCMCprobitChange(formula=Y~X-1, data=parent.frame(), m=1,
                         mcmc=1000, burnin=1000, thin=1, verbose=1000, 
                         b0 = 0, B0 = 10, a = 1, b = 1,  marginal.likelihood = c("Chib95"))
out2 <- MCMCprobitChange(formula=Y~X-1, data=parent.frame(), m=2,
                         mcmc=1000, burnin=1000, thin=1, verbose=1000, 
                         b0 = 0, B0 = 10, a = 1, b = 1,  marginal.likelihood = c("Chib95"))
out3 <- MCMCprobitChange(formula=Y~X-1, data=parent.frame(), m=3,
                         mcmc=1000, burnin=1000, thin=1, verbose=1000, 
                         b0 = 0, B0 = 10, a = 1, b = 1,  marginal.likelihood = c("Chib95"))

## find the most reasonable one
BayesFactor(out0, out1, out2, out3)

## draw plots using the "right" model
plotState(out2)
plotChangepoint(out2)
}
}

\keyword{models}

\seealso{\code{\link{plotState}}, \code{\link{plotChangepoint}}}
