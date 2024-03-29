\name{HMMpanelRE}

\alias{HMMpanelRE}

\title{Markov Chain Monte Carlo for the Hidden Markov Random-effects Model}

\description{HMMpanelRE generates a sample from the posterior
  distribution of the hidden Markov random-effects model discussed in Park (2011).
  The code works for panel data with the same starting point. 
  The sampling of panel parameters is based on Algorithm 2 of Chib and Carlin (1999). This model uses a multivariate
  Normal prior for the fixed effects parameters and varying individual effects, an Inverse-Wishart prior on the 
  random-effects parameters, an Inverse-Gamma prior on the residual error variance, and Beta prior for transition probabilities. 
  The user supplies data and priors, and a sample from the posterior distribution is returned as an
  mcmc object, which can be subsequently analyzed with functions
  provided in the coda package.}

\usage{
HMMpanelRE(subject.id, time.id, y, X, W, m=1,
           mcmc=1000, burnin=1000, thin=1, verbose=0, 
           b0=0, B0=0.001, c0 = 0.001, d0 = 0.001, r0, R0, a = NULL, b = NULL, 
           seed = NA, beta.start = NA, sigma2.start = NA, D.start= NA, P.start = NA, 
           marginal.likelihood = c("none", "Chib95"), ...)}

\arguments{

  \item{subject.id}{A numeric vector indicating the group number. It should start from 1.}

  \item{time.id}{A numeric vector indicating the time unit. It should start from 1.}
  
  \item{y}{The dependent variable}

  \item{X}{The model matrix of the fixed-effects}

  \item{W}{The model matrix of the random-effects. W should be a subset of X.}

  \item{m}{The number of changepoints.}

   \item{mcmc}{The number of MCMC iterations after burn-in.}

  \item{burnin}{The number of burn-in iterations for the sampler.}

  \item{thin}{The thinning interval used in the simulation.  The number of
      MCMC iterations must be divisible by this value.}

  \item{verbose}{A switch which determines whether or not the progress of
    the sampler is printed to the screen.  If \code{verbose} is greater
    than 0, the iteration number and the posterior density samples are printed to the screen every \code{verbose}th iteration.}

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

  \item{r0}{The shape parameter for the Inverse-Wishart prior on variance
    matrix for the random effects. Set r=q for an uninformative prior where q is the number of random effects}
    
  \item{R0}{The scale matrix for the Inverse-Wishart prior on variance matrix
    for the random effects. This must be a square q-dimension
    matrix. Use plausible variance regarding random effects for the
    diagonal of R.}

    \item{a}{\eqn{a}{a} is the shape1 beta prior for transition probabilities. By default,
    the expected duration is computed and corresponding a and b values are assigned. The expected
    duration is the sample period divided by the number of states.}

    \item{b}{\eqn{b}{b} is the shape2 beta prior for transition probabilities. By default,
    the expected duration is computed and corresponding a and b values are assigned. The expected
    duration is the sample period divided by the number of states.}

    \item{seed}{The seed for the random number generator.  If NA, current R system seed is used.}

    \item{beta.start}{The starting values for the beta vector. This can either be a scalar or a column vector with dimension equal to the number of betas. The default value of NA will use draws from the Uniform distribution with the same boundary with the data as the starting value. If this is a scalar, that value will serve as the starting value mean for all of the betas. When there is no covariate, the log value of means should be used.}
	
    \item{sigma2.start}{The starting values for \eqn{\sigma^2}{sigma^2}. This can either be a scalar or a column vector with dimension equal to the number of states.}

    \item{D.start}{The starting values for the beta vector. This can either be a scalar or a column vector with dimension equal to the number of betas. The default value of NA will use draws from the Uniform distribution with the same boundary with the data as the starting value. If this is a scalar, that value will serve as the starting value mean for all of the betas. When there is no covariate, the log value of means should be used.}

    \item{P.start}{The starting values for the transition matrix. A user should provide a square matrix with dimension equal to the number of states. By default, draws from the \code{Beta(0.9, 0.1)} are used to construct a proper transition matrix for each raw except the last raw.}	

  \item{marginal.likelihood}{How should the marginal likelihood be
    calculated? Options are: \code{none} in which case the marginal
    likelihood will not be calculated and \code{Chib95} in which case the method of Chib (1995) is used.}
	
  \item{...}{further arguments to be passed}       

}

\details{
   
  \code{HMMpanelRE} simulates from the random-effect hidden Markov panel model introduced by Park (2011). 
 
 The model takes the following form:
  \deqn{y_i = X_i \beta_m + W_i b_i + \varepsilon_i\;\; m = 1, \ldots, M}{y_i = X_i * beta_m + W_i *
    b_i + epsilon_i, m = 1,..., M.}
  Where each group \eqn{i} have \eqn{k_i} observations.
  Random-effects parameters are assumed to be time-varying at the system level:
  \deqn{b_i \sim \mathcal{N}_q(0, D_m)}{b_i ~ N_q(0, D_m)}
  \deqn{\varepsilon_i \sim \mathcal{N}(0, \sigma^2_m I_{k_i})}{epsilon_i ~ N(0, sigma^2_m I_{k_i})}

  And the errors:
  We assume standard, conjugate priors:
  \deqn{\beta \sim \mathcal{N}_p(b0, B0)}{beta ~ N_p(b0, B0)}
  And:
  \deqn{\sigma^{2} \sim \mathcal{IG}amma(c0/2, d0/2)}{sigma^2 ~ IGamma(c0/2, d0/2)}
  And:
  \deqn{D \sim \mathcal{IW}ishart(r0, R0)}{D ~ IWishart(r0, R0)}
  See Chib and Carlin (1999) for more details.
   
  And:
  \deqn{p_{mm} \sim \mathcal{B}eta(a, b),\;\; m = 1, \ldots, M}{p_mm ~ Beta(a, b), m = 1,...,M.}
  Where \eqn{M}{M} is the number of states.
  
  \emph{NOTE:} We do not provide default parameters for the priors on
   the precision matrix for the random effects. When fitting one of
   these models, it is of utmost importance to choose a prior that
   reflects your prior beliefs about the random effects. Using the
   \code{dwish} and \code{rwish} functions might be useful in choosing
   these values.
	
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
  true.beta1   <-  c(1, 1, 1) ; true.beta2   <-  c(-1, -1, -1)
  true.sigma2 <-  c(2, 5); true.D1 <- diag(.5, Q); true.D2 <- diag(2.5, Q)
  N=30; T=100; 
  NT <- N*T
  x1 <- runif(NT, 1, 2)
  x2 <- runif(NT, 1, 2)
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

  ## model fitting
  G <- 100 
  b0  <- rep(0, K) ; B0  <- solve(diag(100, K))
  c0  <- 2; d0  <- 2
  r0  <- 5; R0  <- diag(c(1, 0.1, 0.1))
  subject.id <- c(rep(1:N, each=T))
  time.id <- c(rep(1:T, N))
  out1 <- HMMpanelRE(subject.id, time.id, y, X, W, m=1, 
                     mcmc=G, burnin=G, thin=1, verbose=G, 
                     b0=b0, B0=B0, c0=c0, d0=d0, r0=r0, R0=R0) 
  
  ## latent state changes   
  plotState(out1)

  ## print mcmc output 	     
  summary(out1)



}
}
\keyword{models}
