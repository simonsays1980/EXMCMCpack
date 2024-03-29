\name{MCMCoprobitChange}
\alias{MCMCoprobitChange}

\title{Markov Chain Monte Carlo for Ordered Probit Changepoint Regression Model}
\description{
  This function generates a sample from the posterior distribution
  of an ordered probit regression model with multiple parameter breaks. The function uses
  the Markov chain Monte Carlo method of Chib (1998).
  The user supplies data and priors, and
  a sample from the posterior distribution is returned as an mcmc
  object, which can be subsequently analyzed with functions
  provided in the coda package.
}

\usage{MCMCoprobitChange(formula, data=parent.frame(),  m = 1,
        burnin = 1000, mcmc = 1000, thin = 1, tune = NA, verbose = 0,
        seed = NA, beta.start = NA, gamma.start=NA, P.start = NA,
        b0 = NULL, B0 = NULL, a = NULL, b = NULL,
        marginal.likelihood = c("none", "Chib95"), gamma.fixed=0, ...)}

\arguments{
   \item{formula}{Model formula.}

    \item{data}{Data frame.}

    \item{m}{The number of changepoints.}

    \item{burnin}{The number of burn-in iterations for the sampler.}

    \item{mcmc}{The number of MCMC iterations after burnin.}

    \item{thin}{The thinning interval used in the simulation.  The number of
      MCMC iterations must be divisible by this value.}

    \item{tune}{The tuning parameter for the Metropolis-Hastings
      step. Default of NA corresponds to a choice of 0.05 divided by the
      number of categories in the response variable.}

    \item{verbose}{A switch which determines whether or not the progress of
    the sampler is printed to the screen.  If \code{verbose} is greater
    than 0 the iteration number, the \eqn{\beta}{beta} vector, and the error variance are printed to
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

    \item{gamma.start}{The starting values for the \eqn{\gamma}{gamma} vector.
    This can either be a scalar or a
    column vector with dimension equal to the number of gammas.
    The default value of of NA will use the MLE
    estimate of \eqn{\gamma}{gamma} as the starting value.  If this is a
    scalar,  that value will serve as the starting value
    mean for all of the gammas.}

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
    
    \item{gamma.fixed}{1 if users want to constrain \eqn{\gamma}{gamma} values to be constant. By default, 
     \eqn{\gamma}{gamma} values are allowed to vary across regimes.}

    \item{...}{further arguments to be passed}
}

\value{
   An mcmc object that contains the posterior sample.  This
   object can be summarized by functions provided by the coda package.
   The object contains an attribute \code{prob.state} storage matrix that contains the probability of \eqn{state_i}{state_i} 
   for each period, the log-likelihood of the model (\code{loglike}), and
   the log-marginal likelihood of the model (\code{logmarglike}). 
}

\details{
  \code{MCMCoprobitChange} simulates from the posterior distribution of
  an ordinal probit regression model with multiple parameter breaks. The simulation of latent states is based on 
  the linear approximation method discussed in Park (2011). 

  The model takes the following form:
  \deqn{\Pr(y_t = 1) = \Phi(\gamma_{c, m} - x_i'\beta_m) - \Phi(\gamma_{c-1, m} - x_i'\beta_m)\;\; m = 1, \ldots, M}{
   Pr(y_t = 1) = Phi(gamma_(c, m) - x_i'beta_m) - Phi(gamma_(c-1, m) - x_i'beta)
   }
  Where \eqn{M}{M} is the number of states, and \eqn{\gamma_{c, m}}{gamma_(c, m)} and \eqn{\beta_m}{beta_m} 
  are paramters when a state is \eqn{m}{m} at \eqn{t}{t}. 

  We assume Gaussian distribution for prior of \eqn{\beta}{beta}:
  \deqn{\beta_m \sim \mathcal{N}(b_0,B_0^{-1}),\;\; m = 1, \ldots, M}{beta_m ~ N(b0,B0^(-1)), m = 1,...,M.}

  And:
  \deqn{p_{mm} \sim \mathcal{B}eta(a, b),\;\; m = 1, \ldots, M}{p_mm ~ Beta(a, b), m = 1,...,M.}
  Where \eqn{M}{M} is the number of states.


  Note that when the fitted changepoint model has very few observations in any of states, 
  the marginal likelihood outcome can be ``nan," which indicates that too many breaks are assumed
  given the model and data. 

}

\references{ 
  Jong Hee Park. 2011. ``Changepoint Analysis of Binary and Ordinal Probit Models:
  An Application to Bank Rate Policy Under the Interwar Gold Standard." \emph{Political Analysis}. 19: 188-204.

 Andrew D. Martin, Kevin M. Quinn, and Jong Hee Park. 2011.
  ``MCMCpack: Markov Chain Monte Carlo in R.'',
  \emph{Journal of Statistical Software}. 42(9): 1-21.
  \url{http://www.jstatsoft.org/v42/i09/}.

  Siddhartha Chib. 1998. ``Estimation and comparison of multiple change-point models.''
  \emph{Journal of Econometrics}. 86: 221-241.

}

\examples{
set.seed(1909)
N <- 200
x1 <- rnorm(N, 1, .5);

## set a true break at 100
z1 <- 1 + x1[1:100] + rnorm(100);
z2 <- 1 -0.2*x1[101:200] + rnorm(100);
z <- c(z1,  z2);
y <- z

## generate y
y[z < 1] <- 1;
y[z >= 1 & z < 2] <- 2;
y[z >= 2] <- 3;

## inputs
formula <- y ~ x1 

## fit multiple models with a varying number of breaks
out1 <- MCMCoprobitChange(formula, m=1, 
      	mcmc=1000, burnin=1000, thin=1, tune=c(.5, .5), verbose=1000, 
     	b0=0, B0=10, marginal.likelihood = "Chib95")
out2 <- MCMCoprobitChange(formula, m=2, 
      	mcmc=1000, burnin=1000, thin=1, tune=c(.5, .5, .5), verbose=1000, 
     	b0=0, B0=10, marginal.likelihood = "Chib95")
out3 <- MCMCoprobitChange(formula, m=3, 
      	mcmc=1000, burnin=1000, thin=1, tune=c(.5, .5, .5, .5), verbose=1000, 
     	b0=0, B0=10, marginal.likelihood = "Chib95")

## find the most reasonable one
BayesFactor(out1, out2, out3)

## draw plots using the "right" model
plotState(out1)
plotChangepoint(out1)
}

\keyword{models}

\seealso{\code{\link{plotState}}, \code{\link{plotChangepoint}}}
