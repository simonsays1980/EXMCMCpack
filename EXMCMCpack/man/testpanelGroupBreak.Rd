\name{testpanelGroupBreak}

\alias{testpanelGroupBreak}

\title{A Test for the Group-level Break using a Multivariate Linear Regression Model with Breaks}

\description{testpanelGroupBreak fits a multivariate linear regression model with parametric breaks using 
  panel residuals to test the existence of group-level breaks in panel residuals. The details are discussed in Park (2011).}

\usage{
testpanelGroupBreak(subject.id, time.id, resid, m=1, 
           mcmc=1000, burnin=1000, thin=1, verbose=0, 
           b0, B0, c0, d0, a = NULL, b = NULL, 
           seed = NA, marginal.likelihood = c("none", "Chib95"), ...)}

\arguments{

  \item{subject.id}{A numeric vector indicating the group number. It should start from 1.}

  \item{time.id}{A numeric vector indicating the time unit. It should start from 1.}
  
  \item{resid}{A vector of panel residuals}

  \item{m}{The number of changepoints.}

   \item{mcmc}{The number of MCMC iterations after burn-in.}

  \item{burnin}{The number of burn-in iterations for the sampler.}

  \item{thin}{The thinning interval used in the simulation.  The number of
      MCMC iterations must be divisible by this value.}

  \item{verbose}{A switch which determines whether or not the progress of
    the sampler is printed to the screen.  If \code{verbose} is greater
    than 0, the iteration number and the posterior density samples are printed to the screen every \code{verbose}th iteration.}

    \item{b0}{The prior mean of the residual mean.}

    \item{B0}{The prior precision of the residual variance}

   \item{c0}{\eqn{c_0/2}{c0/2} is the shape parameter for the inverse
    Gamma prior on \eqn{\sigma^2}{sigma^2}. The amount of information in the inverse Gamma prior
    is something like that from \eqn{c_0}{c0} pseudo-observations.} 
    
    \item{d0}{\eqn{d_0/2}{d0/2} is the scale parameter for the
    inverse Gamma prior on \eqn{\sigma^2}{sigma^2}.}

    \item{a}{\eqn{a}{a} is the shape1 beta prior for transition probabilities. By default,
    the expected duration is computed and corresponding a and b values are assigned. The expected
    duration is the sample period divided by the number of states.}

    \item{b}{\eqn{b}{b} is the shape2 beta prior for transition probabilities. By default,
    the expected duration is computed and corresponding a and b values are assigned. The expected
    duration is the sample period divided by the number of states.}

    \item{seed}{The seed for the random number generator.  If NA, current R system seed is used.}

  \item{marginal.likelihood}{How should the marginal likelihood be
    calculated? Options are: \code{none} in which case the marginal
    likelihood will not be calculated and \code{Chib95} in which case the method of Chib (1995) is used.}
	
  \item{...}{further arguments to be passed}       

}

\details{
   
  \code{testpanelGroupBreak} fits a multivariate linear regression model with parametric breaks using panel residuals 
  to detect the existence of system-level breaks in unobserved factors as discussed in Park (2011). 
 
  The model takes the following form:
  \deqn{e_{i} \sim \mathcal{N}(\beta_{m}, \sigma^2_m I)\;\; m = 1, \ldots, M}{epsilon_i ~ N(beta_m, sigma^2_m I_{k_i}), m = 1,..., M.} 
 
  We assume standard, semi-conjugate priors:
  \deqn{\beta \sim \mathcal{N}(b0, B0)}{beta ~ N(b0, B0)}
  And:
  \deqn{\sigma^{-2} \sim \mathcal{G}amma(c_0/2, d_0/2)}{sigma^(-2) ~
    Gamma(c0/2, d0/2)}
  Where \eqn{\beta}{beta} and \eqn{\sigma^{-2}}{sigma^(-2)} are assumed 
  \emph{a priori} independent.  
  
  And:
  \deqn{p_{mm} \sim \mathcal{B}eta(a, b),\;\; m = 1, \ldots, M}{p_mm ~ Beta(a, b), m = 1,...,M.}
  Where \eqn{M}{M} is the number of states.
  	
 }

\value{
   An mcmc object that contains the posterior sample. This
   object can be summarized by functions provided by the coda package.
   The object contains an attribute \code{prob.state} storage matrix that contains the probability of \eqn{state_i}{state_i} 
   for each period, and
   the log-marginal likelihood of the model (\code{logmarglike}).
}

\references{
   Jong Hee Park, 2011. ``A Unified Method for Dynamic and Cross-Sectional Heterogeneity: 
   Introducing Hidden Markov Panel Models." Working Paper.
   
   Siddhartha Chib. 1998. ``Estimation and comparison of multiple change-point models.''
   \emph{Journal of Econometrics}. 86: 221-241.
	  
}

\author{Jong Hee Park, \email{jhp@uchicago.edu}, \url{http://home.uchicago.edu/~jhp/}.}

\examples{
\dontrun{
   ## data generating
  set.seed(1977)
  Q <- 3
  true.beta1   <-  c(1, 1, 1) ; true.beta2   <-  c(1, -1, -1)
  true.sigma2 <-  c(1, 3); true.D1 <- diag(.5, Q); true.D2 <- diag(2.5, Q)
  N=20; T=100; 
  NT <- N*T
  x1 <- rnorm(NT)
  x2 <- runif(NT, 5, 10)
  X <- cbind(1, x1, x2);   W <- X;   y <- rep(NA, NT)

  ## true break numbers are one and at the center
  break.point = rep(T/2, N); break.sigma=c(rep(1, N));
  break.list <- rep(1, N)
  id  <-  rep(1:N, each=NT/N)
  K <-  ncol(X);   
  ruler <- c(1:T)

  ## compute the weight for the break
  W.mat <- matrix(NA, T, N)
  for (i in 1:N){
    W.mat[, i] <- pnorm((ruler-break.point[i])/break.sigma[i])
  }
  Weight <- as.vector(W.mat)

  ## data generating by weighting two means and variances
  j = 1
  for (i in 1:N){
    Xi <- X[j:(j+T-1), ]
    Wi <- W[j:(j+T-1), ]
    true.V1 <- true.sigma2[1]*diag(T) + Wi\%*\%true.D1\%*\%t(Wi)
    true.V2 <- true.sigma2[2]*diag(T) + Wi\%*\%true.D2\%*\%t(Wi)
    true.mean1 <- Xi\%*\%true.beta1
    true.mean2 <- Xi\%*\%true.beta2
    weight <- Weight[j:(j+T-1)]
    y[j:(j+T-1)] <- (1-weight)*true.mean1 + (1-weight)*chol(true.V1)\%*\%rnorm(T) +
      weight*true.mean2 + weight*chol(true.V2)\%*\%rnorm(T) 
    j <- j + T
  }
  ## model fitting
  subject.id <- c(rep(1:N, each=T))
  time.id <- c(rep(1:T, N))


  resid <- rstandard(lm(y ~X-1 + as.factor(subject.id)))
  G <- 100
  out0 <- testpanelGroupBreak(subject.id, time.id, resid, m=0, 
           mcmc=G, burnin=G, thin=1, verbose=G,	   
           b0=0, B0=1/100, c0=2, d0=2, marginal.likelihood = "Chib95")
  out1 <- testpanelGroupBreak(subject.id, time.id, resid, m=1, 
           mcmc=G, burnin=G, thin=1, verbose=G,	   
           b0=0, B0=1/100, c0=2, d0=2, marginal.likelihood = "Chib95")
  out2 <- testpanelGroupBreak(subject.id, time.id, resid, m=2, 
           mcmc=G, burnin=G, thin=1, verbose=G,	   
           b0=0, B0=1/100, c0=2, d0=2, marginal.likelihood = "Chib95")
  out3 <- testpanelGroupBreak(subject.id, time.id, resid, m=3, 
           mcmc=G, burnin=G, thin=1, verbose=G,	   
           b0=0, B0=1/100, c0=2, d0=2, marginal.likelihood = "Chib95")
  
  ## Note that the code is for a hypothesis test of no break in panel residuals. 	     
  ## When breaks exist, the estimated number of break in the mean and variance of panel residuals 
  ## tends to be larger than the number of break in the data generating process.
  ## This is due to the difference in parameter space, not an error of the code.
  BayesFactor(out0, out1, out2, out3)

  ## In order to identify the number of breaks in panel parameters,
  ## use HMMpanelRE() instead. 

}
}
\keyword{models}
