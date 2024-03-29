\name{HMMpanelFE}

\alias{HMMpanelFE}

\title{Markov Chain Monte Carlo for the Hidden Markov Fixed-effects Model}

\description{HMMpanelFE generates a sample from the posterior
  distribution of the fixed-effects model with varying individual effects model discussed in Park (2011).
  The code works for both balanced and unbalanced panel data as long as there is no missing data in the middle of each group. 
  This model uses a multivariate
  Normal prior for the fixed effects parameters and varying individual effects, an Inverse-Gamma
  prior on the residual error variance, and Beta prior for transition probabilities. The user supplies data and
  priors, and a sample from the posterior distribution is returned as an
  mcmc object, which can be subsequently analyzed with functions
  provided in the coda package.}

\usage{HMMpanelFE(subject.id, y, X, m,
           mcmc=1000, burnin=1000, thin=1, verbose=0, 
           b0=0, B0=0.001, c0 = 0.001, d0 = 0.001, delta0=0, Delta0=0.001,
           a = NULL, b = NULL, seed = NA, ...)}

\arguments{

  \item{subject.id}{A numeric vector indicating the group number. It should start from 1.}
  
  \item{y}{The response variable.}

  \item{X}{The model matrix excluding the constant.}

  \item{m}{A vector of break numbers for each subject in the panel.}

 
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

     \item{delta0}{The prior mean of \eqn{\alpha}{alpha}.}

    \item{Delta0}{The prior precision of \eqn{\alpha}{alpha}.}

    \item{a}{\eqn{a}{a} is the shape1 beta prior for transition probabilities. By default,
    the expected duration is computed and corresponding a and b values are assigned. The expected
    duration is the sample period divided by the number of states.}

    \item{b}{\eqn{b}{b} is the shape2 beta prior for transition probabilities. By default,
    the expected duration is computed and corresponding a and b values are assigned. The expected
    duration is the sample period divided by the number of states.}

    \item{seed}{The seed for the random number generator.  If NA, current R system seed is used.}
	
   \item{...}{further arguments to be passed}       

}

\details{
  
  \code{HMMpanelFE} simulates from the fixed-effect hidden Markov panel model introduced by Park (2011). 
  The model takes the following form:
  \deqn{y_{it} = x'_{it} \beta + \varepsilon_{it}\;\; m = 1, \ldots, M}{y_it = x'_it * beta + epsilon_it, m = 1,...,M.}
  Unlike conventional fixed-effects models, individual effects and variances are assumed to be time-varying at the subject level:
  \deqn{\varepsilon_{it} \sim \mathcal{N}(\alpha_{im}, \sigma^2_{im})}{epsilon_it ~ N(alpha_im, sigma^2_im)} 

   We assume standard, semi-conjugate priors:
  \deqn{\beta \sim \mathcal{N}(b_0,B_0^{-1})}{beta ~ N(b0,B0^(-1))}
  And:
  \deqn{\sigma^{-2} \sim \mathcal{G}amma(c_0/2, d_0/2)}{sigma^(-2) ~
    Gamma(c0/2, d0/2)}
  And:
  \deqn{\alpha \sim \mathcal{N}(delta_0,Delta_0^{-1})}{alpha ~ N(delta0, Delta0^(-1))}
  \eqn{\beta}{beta}, \eqn{\alpha}{alpha} and \eqn{\sigma^{-2}}{sigma^(-2)} are assumed 
  \emph{a priori} independent.  
  
  And:
  \deqn{p_{mm} \sim \mathcal{B}eta(a, b),\;\; m = 1, \ldots, M}{p_mm ~ Beta(a, b), m = 1,...,M.}
  Where \eqn{M}{M} is the number of states.
  
  OLS estimates are used for starting values.

  }

\author{Jong Hee Park, \email{jhp@uchicago.edu}, \url{http://home.uchicago.edu/~jhp/}.}

\value{
   An mcmc object that contains the posterior sample. This
   object can be summarized by functions provided by the coda package.
   The object contains an attribute \code{sigma} storage matrix that contains time-varying residual variance, an attribute \code{state} storage   
   matrix that contains posterior samples of hidden states, and an attribute \code{delta} storage matrix containing time-varying intercepts.
}

\references{
   Jong Hee Park, 2011. ``A Unified Method for Dynamic and Cross-Sectional Heterogeneity: 
   Introducing Hidden Markov Panel Models." Working Paper.
   
   Siddhartha Chib. 1998. ``Estimation and comparison of multiple change-point models.''
   \emph{Journal of Econometrics}. 86: 221-241.
	  
}

\examples{
\dontrun{
  ## data generating
  set.seed(1974)
  N <- 30
  T <- 80
  NT <- N*T

  ## true parameter values
  true.beta <- c(1, 1)
  true.sigma <- 3
  x1 <- rnorm(NT)
  x2 <- runif(NT, 2, 4)

  ## group-specific breaks 
  break.point = rep(T/2, N); break.sigma=c(rep(1, N));
  break.list <- rep(1, N)

  X <- as.matrix(cbind(x1, x2), NT, );
  y <- rep(NA, NT)
  id  <-  rep(1:N, each=NT/N)
  K <-  ncol(X);   
  true.beta <- as.matrix(true.beta, K, 1)

  ## compute the break probability
  ruler <- c(1:T)
  W.mat <- matrix(NA, T, N)
  for (i in 1:N){
    W.mat[, i] <- pnorm((ruler-break.point[i])/break.sigma[i])
  }
  Weight <- as.vector(W.mat)

  ## draw time-varying individual effects and sample y
  j = 1
  true.sigma.alpha <- 30	   
  true.alpha1 <- true.alpha2 <- rep(NA, N) 
  for (i in 1:N){
    Xi <- X[j:(j+T-1), ]
    true.mean <- Xi  \%*\% true.beta
    weight <- Weight[j:(j+T-1)]
    true.alpha1[i] <- rnorm(1, 0, true.sigma.alpha)
    true.alpha2[i] <- -1*true.alpha1[i]
    y[j:(j+T-1)] <- ((1-weight)*true.mean + (1-weight)*rnorm(T, 0, true.sigma) + (1-weight)*true.alpha1[i]) +
    		     (weight*true.mean + weight*rnorm(T, 0, true.sigma) + weight*true.alpha2[i])
    j <- j + T
  }

  ## extract the standardized residuals from the OLS with fixed-effects
  FEols <- lm(y ~ X + as.factor(id) -1 )
  resid.all <- rstandard(FEols)
  time.id <- rep(1:80, N)

  ## model fitting
  G <- 100
  BF <- testpanelSubjectBreak(subject.id=id, time.id=time.id,
         resid= resid.all, max.break=3, minimum = 10, 
         mcmc=G, burnin = G, thin=1, verbose=G, 
         b0=0, B0=1/100, c0=2, d0=2, Time = time.id)

  ## get the estimated break numbers
  estimated.breaks <- make.breaklist(BF, threshold=3)

  ## model fitting 
  out <- HMMpanelFE(subject.id = id, y, X=X, m =  estimated.breaks,
             mcmc=G, burnin=G, thin=1, verbose=G, 
             b0=0, B0=1/1000, c0=2, d0=2, delta0=0, Delta0=1/1000)

  ## print out the slope estimate
  ## true values are 1 and 1	     
  summary(out)
  
  ## compare them with the result from the constant fixed-effects 
  summary(FEols)
}
}
\keyword{models}
