\name{MCMCpoissonChange}
\alias{MCMCpoissonChange}

\title{Markov Chain Monte Carlo for a Poisson Regression Changepoint Model}
\description{
  This function generates a sample from the posterior distribution
  of a Poisson regression model with multiple changepoints. The function uses
  the Markov chain Monte Carlo method of Chib (1998).
  The user supplies data and priors, and
  a sample from the posterior distribution is returned as an mcmc
  object, which can be subsequently analyzed with functions
  provided in the coda package.
}

\usage{MCMCpoissonChange(
		formula, data = parent.frame(), m = 1,
           	b0 = 0, B0 = 1, a = NULL, b = NULL, c0 = NA, d0 = NA, 
           	burnin = 1000, mcmc = 1000, thin = 1, verbose = 0, 
           	seed = NA, beta.start = NA, P.start = NA,   
           	marginal.likelihood = c("none", "Chib95"), ...)}

\arguments{
    \item{formula}{Model formula.}

    \item{data}{Data frame.}

    \item{m}{The number of changepoints.}

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

     \item{c0}{\eqn{c_0}{c0} is the shape parameter for Gamma prior on \eqn{\lambda}{lambda}
    (the mean). When there is no covariate, this should be provided by users. No default value is provided.}

    \item{d0}{\eqn{d_0}{d0} is the scale parameter for Gamma prior on \eqn{\lambda}{lambda}
    (the mean). When there is no covariate, this should be provided by users. No default value is provided.}

    \item{burnin}{The number of burn-in iterations for the sampler.}

    \item{mcmc}{The number of MCMC iterations after burn-in.}

    \item{thin}{The thinning interval used in the simulation.  The number of
      MCMC iterations must be divisible by this value.}

    \item{verbose}{A switch which determines whether or not the progress of
    the sampler is printed to the screen.  If \code{verbose} is greater
    than 0, the iteration number and the posterior density samples are printed to the screen every \code{verbose}th iteration.}

    \item{seed}{The seed for the random number generator.  If NA, current R system seed is used.}
	
    \item{beta.start}{The starting values for the beta vector. This can either be a scalar or a column vector with dimension equal to the number of betas. The default value of NA will use draws from the Uniform distribution with the same boundary with the data as the starting value. If this is a scalar, that value will serve as the starting value mean for all of the betas. When there is no covariate, the log value of means should be used.}
	
    \item{P.start}{The starting values for the transition matrix. A user should provide a square matrix with dimension equal to the number of states. By default, draws from the \code{Beta(0.9, 0.1)} are used to construct a proper transition matrix for each raw except the last raw.}	

    \item{marginal.likelihood}{How should the marginal likelihood be
    calculated? Options are: \code{none} in which case the marginal
    likelihood will not be calculated, and
    \code{Chib95} in which case the method of Chib (1995) is used.}

    \item{...}{further arguments to be passed}
}

\value{
   An mcmc object that contains the posterior sample. This
   object can be summarized by functions provided by the coda package.
   The object contains an attribute \code{prob.state} storage matrix that contains the probability of \eqn{state_i}{state_i} for each period, and
   the log-marginal likelihood of the model (\code{logmarglike}).
}

\details{
  \code{MCMCpoissonChange} simulates from the posterior distribution of
  a Poisson regression model with multiple changepoints using the methods of Chib (1998) and Fruhwirth-Schnatter and Wagner (2006).
  The details of the model are discussed in Park (2010). 

  The model takes the following form:
  \deqn{y_t \sim \mathcal{P}oisson(\mu_t)}{y_t ~ Poisson(mu_t)}
  \deqn{\mu_t = x_t ' \beta_m,\;\; m = 1, \ldots, M}{mu_t = x_t'beta_m,  m = 1,...,M.}
  Where \eqn{M}{M} is the number of states and \eqn{\beta_m}{beta_m} is paramters when a state is \eqn{m}{m} at \eqn{t}{t}. 

  We assume Gaussian distribution for prior of \eqn{\beta}{beta}:
  \deqn{\beta_m \sim \mathcal{N}(b_0,B_0^{-1}),\;\; m = 1, \ldots, M}{beta_m ~ N(b0,B0^(-1)), m = 1,...,M.}

  And:
  \deqn{p_{mm} \sim \mathcal{B}eta(a, b),\;\; m = 1, \ldots, M}{p_mm ~ Beta(a, b), m = 1,...,M.}
  Where \eqn{M}{M} is the number of states.
  }

\author{Jong Hee Park, \email{jhp@uchicago.edu}, \url{http://home.uchicago.edu/~jhp/}.}

\references{
 Jong Hee Park. 2010. ``Structural Change in the U.S. Presidents' Use of Force Abroad.''
 \emph{American Journal of Political Science} 54: 766-782.

 Sylvia Fruhwirth-Schnatter and Helga Wagner 2006. ``Auxiliary Mixture Sampling for Parameter-driven Models of 
 Time Series of Counts with Applications to State Space Modelling.'' \emph{Biometrika}. 93:827--841.

 Siddhartha Chib. 1998. ``Estimation and comparison of multiple change-point models.''
   \emph{Journal of Econometrics}. 86: 221-241.

Andrew D. Martin, Kevin M. Quinn, and Jong Hee Park. 2011.
 ``MCMCpack: Markov Chain Monte Carlo in R.'',
 \emph{Journal of Statistical Software}. 42(9): 1-21.
 \url{http://www.jstatsoft.org/v42/i09/}.

 Siddhartha Chib. 1995. ``Marginal Likelihood from the Gibbs Output.''
   \emph{Journal of the American Statistical Association}. 90:
   1313-1321.
}

\examples{
    \dontrun{
    set.seed(11119)
    n <- 150
    x1 <- runif(n, 0, 0.5)
    true.beta1 <- c(1,  1)
    true.beta2 <- c(1,  -2)
    true.beta3 <- c(1,  2)    
    
    ## set true two breaks at (50, 100)
    true.s <- rep(1:3, each=n/3)
    mu1 <- exp(1 + x1[true.s==1]*1)
    mu2 <- exp(1 + x1[true.s==2]*-2)
    mu3 <- exp(1 + x1[true.s==3]*2)

    y <- as.ts(c(rpois(n/3, mu1), rpois(n/3, mu2), rpois(n/3, mu3)))
    formula = y ~ x1
    
    ## fit multiple models with a varying number of breaks
    model0 <-  MCMCpoissonChange(formula, m=0, 
    	       mcmc = 1000, burnin = 1000, verbose = 500, 
	       b0 = rep(0, 2), B0 = 5*diag(2), marginal.likelihood = "Chib95")    
    model1 <-  MCMCpoissonChange(formula, m=1, 
    	       mcmc = 1000, burnin = 1000, verbose = 500, 
	       b0 = rep(0, 2), B0 = 5*diag(2), marginal.likelihood = "Chib95")    
    model2 <-  MCMCpoissonChange(formula, m=2,
    	       mcmc = 1000, burnin = 1000, verbose = 500, 
	       b0 = rep(0, 2), B0 = 5*diag(2), marginal.likelihood = "Chib95")    
    model3 <-  MCMCpoissonChange(formula, m=3, 
    	       mcmc = 1000, burnin = 1000, verbose = 500, 
	       b0 = rep(0, 2), B0 = 5*diag(2), marginal.likelihood = "Chib95")    
    model4 <-  MCMCpoissonChange(formula, m=4, 
    	       mcmc = 1000, burnin = 1000, verbose = 500, 
	       b0 = rep(0, 2), B0 = 5*diag(2), marginal.likelihood = "Chib95")    
    model5 <-  MCMCpoissonChange(formula, m=5, 
    	       mcmc = 1000, burnin = 1000, verbose = 500, 
	       b0 = rep(0, 2), B0 = 5*diag(2), marginal.likelihood = "Chib95")    

    ## find the most reasonable one
    print(BayesFactor(model0, model1, model2, model3, model4, model5))

    ## draw plots using the "right" model
    par(mfrow=c(attr(model2, "m") + 1, 1), mai=c(0.4, 0.6, 0.3, 0.05))
    plotState(model2, legend.control = c(1, 0.6))
    plotChangepoint(model2, verbose = TRUE, ylab="Density", start=1, overlay=TRUE)
    
    ## No covariate case
    model2.1 <- MCMCpoissonChange(y ~ 1, m = 2, c0 = 2, d0 = 1,
    	       	mcmc = 1000, burnin = 1000, verbose = 500, 
	       	marginal.likelihood = "Chib95")   
    print(BayesFactor(model2, model2.1)) 
    }			    
}

\keyword{models}

\seealso{\code{\link{MCMCbinaryChange}}, \code{\link{plotState}}, \code{\link{plotChangepoint}}}
