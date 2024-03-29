\name{MCMCbinaryChange}
\alias{MCMCbinaryChange}

\title{Markov Chain Monte Carlo for a Binary Multiple Changepoint Model}
\description{
  This function generates a sample from the posterior distribution
  of a binary model with multiple changepoints. The function uses
  the Markov chain Monte Carlo method of Chib (1998).
  The user supplies data and priors, and
  a sample from the posterior distribution is returned as an mcmc
  object, which can be subsequently analyzed with functions
  provided in the coda package.
}

\usage{MCMCbinaryChange(data,  m = 1, c0 = 1,  d0 = 1,  a = NULL, b = NULL,
            burnin = 10000, mcmc = 10000, thin = 1, verbose = 0,
            seed = NA, phi.start = NA, P.start = NA,
            marginal.likelihood = c("none", "Chib95"), ...)}

\arguments{
    \item{data}{The data.}

    \item{m}{The number of changepoints.}

    \item{c0}{\eqn{c_0}{c0} is the shape1 parameter for Beta prior on \eqn{\phi}{phi}
    (the mean).}

    \item{d0}{\eqn{d_0}{d0} is the shape2 parameter for Beta prior on \eqn{\phi}{phi}
    (the mean).}

    \item{a}{\eqn{a}{a} is the shape1 beta prior for transition probabilities. By default,
    the expected duration is computed and corresponding a and b values are assigned. The expected
    duration is the sample period divided by the number of states.}

    \item{b}{\eqn{b}{b} is the shape2 beta prior for transition probabilities. By default,
    the expected duration is computed and corresponding a and b values are assigned. The expected
    duration is the sample period divided by the number of states.}

    \item{burnin}{The number of burn-in iterations for the sampler.}

    \item{mcmc}{The number of MCMC iterations after burn-in.}

    \item{thin}{The thinning interval used in the simulation.  The number of
      MCMC iterations must be divisible by this value.}

    \item{verbose}{A switch which determines whether or not the progress of
    the sampler is printed to the screen.  If \code{verbose} is greater
    than 0, the iteration number and the posterior density samples are printed to the screen every \code{verbose}th iteration.}

    \item{seed}{The seed for the random number generator.  If NA, current R
    system seed is used.}

    \item{phi.start}{The starting values for the mean. The default value of NA will use draws from the Uniform distribution.}
	
    \item{P.start}{The starting values for the transition matrix. A user should provide a square matrix with dimension equal to the number of states. By default, draws from the \code{Beta(0.9, 0.1)} are used to construct a proper transition matrix for each raw except the last raw.}

    \item{marginal.likelihood}{How should the marginal likelihood be
    calculated? Options are: \code{none} in which case the marginal
    likelihood will not be calculated, and
    \code{Chib95} in which case the method of Chib (1995) is used.}

    \item{...}{further arguments to be passed}
}

\value{
   An mcmc object that contains the posterior sample.  This
   object can be summarized by functions provided by the coda package.
   The object contains an attribute \code{prob.state} storage matrix that contains the probability of \eqn{state_i}{state_i} for each period, and
   the log-marginal likelihood of the model (\code{logmarglike}).
}

\details{
  \code{MCMCbinaryChange} simulates from the posterior distribution of
  a binary model with multiple changepoints.

  The model takes the following form:
  \deqn{Y_t \sim \mathcal{B}ernoulli(\phi_i),\;\; i = 1, \ldots, k}{Y_t ~ Bernoulli(phi_i), i = 1,...,k.}
  Where \eqn{k}{k} is the number of states.

  We assume Beta priors for \eqn{\phi_{i}}{phi_i} and for transition probabilities:
  \deqn{\phi_i \sim \mathcal{B}eta(c_0, d_0)}{phi_i ~ Beta(c0, d0)}
  And:
  \deqn{p_{mm} \sim \mathcal{B}eta{a}{b},\;\; m = 1, \ldots, k}{p_mm ~ Beta(a, b), m = 1,...,M.}
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

 Martyn Plummer, Nicky Best, Kate Cowles, and Karen Vines. 2002.
   \emph{Output Analysis and Diagnostics for MCMC (CODA)}.
   \url{http://www-fis.iarc.fr/coda/}.

  Siddhartha Chib. 1995. ``Marginal Likelihood from the Gibbs Output.''
   \emph{Journal of the American Statistical Association}. 90:
   1313-1321.

}

\examples{
    \dontrun{
    set.seed(19173)
    true.phi<- c(0.5, 0.8, 0.4)
    
    ## two breaks at c(80, 180)  
    y1 <- rbinom(80, 1,  true.phi[1])
    y2 <- rbinom(100, 1, true.phi[2])
    y3 <- rbinom(120, 1, true.phi[3])
    y  <- as.ts(c(y1, y2, y3))

    model0 <- MCMCbinaryChange(y, m=0, c0=2, d0=2, mcmc=1000, burnin=1000, verbose=500, marginal.likelihood = "Chib95")
    model1 <- MCMCbinaryChange(y, m=1, c0=2, d0=2, mcmc=1000, burnin=1000, verbose=500, marginal.likelihood = "Chib95")
    model2 <- MCMCbinaryChange(y, m=2, c0=2, d0=2, mcmc=1000, burnin=1000, verbose=500, marginal.likelihood = "Chib95")
    model3 <- MCMCbinaryChange(y, m=3, c0=2, d0=2, mcmc=1000, burnin=1000, verbose=500, marginal.likelihood = "Chib95")
    model4 <- MCMCbinaryChange(y, m=4, c0=2, d0=2, mcmc=1000, burnin=1000, verbose=500, marginal.likelihood = "Chib95")
    model5 <- MCMCbinaryChange(y, m=5, c0=2, d0=2, mcmc=1000, burnin=1000, verbose=500, marginal.likelihood = "Chib95")

    print(BayesFactor(model0, model1, model2, model3, model4, model5))
    
    ## plot two plots in one screen
    par(mfrow=c(attr(model2, "m") + 1, 1), mai=c(0.4, 0.6, 0.3, 0.05))
    plotState(model2, legend.control = c(1, 0.6))
    plotChangepoint(model2, verbose = TRUE, ylab="Density", start=1, overlay=TRUE)

    }
}

\keyword{models}

\seealso{\code{\link{MCMCpoissonChange}},\code{\link{plotState}}, \code{\link{plotChangepoint}}}
