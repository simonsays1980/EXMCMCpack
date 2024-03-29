\name{MCMCdynamicIRT1d}
\alias{MCMCdynamicIRT1d}
\alias{MCMCdynamicIRT1d_b}
\title{Markov Chain Monte Carlo for Dynamic One Dimensional Item Response Theory Model} 

\description{ This function generates a sample from the posterior
 distribution of a dynamic one dimensional item response theory (IRT)
 model, with Normal random walk priors on the subject abilities (ideal
 points), and multivariate Normal priors on the item parameters. The
 user supplies data and priors, and a sample from the posterior
 distribution is returned as an mcmc object, which can be subsequently
 analyzed with functions provided in the coda package. }

\usage{
MCMCdynamicIRT1d(datamatrix, item.time.map, theta.constraints=list(), 
	         burnin = 1000, mcmc = 20000, thin = 1, verbose = 0, 
                 seed = NA, theta.start = NA, alpha.start = NA, 
                 beta.start = NA, tau2.start = 1, a0 = 0, A0 = 0.1, 
                 b0 = 0, B0 = 0.1, c0 = -1, d0 = -1, e0 = 0, 
                 E0 = 1, store.ability = TRUE, 
                 store.item = TRUE, ...)
}

\arguments{
  \item{datamatrix}{The matrix of data.  Must be 0, 1, or missing values.  
     The rows of \code{datamatrix} correspond to subjects and the
     columns correspond to items.}

  \item{item.time.map}{A vector that relates each item to a time period. 
        Each element of \code{item.time.map} gives the time period of the 
        corresponding column of \code{datamatrix}. It is assumed that the 
        minimum value of \code{item.time.map} is 1.}

  \item{theta.constraints}{ A list specifying possible simple equality
    or inequality constraints on the ability parameters. A typical
    entry in the list has one of three forms: \code{varname=c} which
    will constrain the ability parameter for the subject named
    \code{varname} to be equal to c, \code{varname="+"} which will
    constrain the ability parameter for the subject named \code{varname}
    to be positive, and  \code{varname="-"} which will constrain the
    ability parameter for the subject named \code{varname} to be
    negative. If x is a matrix without row names defaults names of
    ``V1",``V2", ... , etc will be used. See Rivers (2003) for a
    thorough discussion of identification of IRT models. }

  \item{burnin}{ The number of burn-in iterations for the sampler. }

  \item{mcmc}{The number of Gibbs iterations for the sampler. }

  \item{thin}{The thinning interval used in the simulation.  The number of
    Gibbs iterations must be divisible by this value. }

  \item{verbose}{ A switch which determines whether or not the progress of
      the sampler is printed to the screen.   If \code{verbose} is greater
      than 0 then every \code{verbose}th iteration will be printed to the
      screen. }

  \item{seed}{ The seed for the random number generator.  If NA, the Mersenne
    Twister generator is used with default seed 12345; if an integer is 
    passed it is used to seed the Mersenne twister.  The user can also
    pass a list of length two to use the L'Ecuyer random number generator,
    which is suitable for parallel computation.  The first element of the
    list is the L'Ecuyer seed, which is a vector of length six or NA (if NA 
    a default seed of \code{rep(12345,6)} is used).  The second element of 
    list is a positive substream number. See the MCMCpack 
    specification for more details. }

  \item{theta.start}{ The starting values for the subject
    abilities (ideal points). This can either be a scalar or a
    column vector with dimension equal to the number of voters.  
    If this takes a scalar value, then that value will serve as
    the starting value for all of the thetas.  The default value
    of NA will choose the starting values based on an
    eigenvalue-eigenvector decomposition of the aggreement score matrix
    formed from the \code{datamatrix}. }

  \item{alpha.start}{ The starting values for the
    \eqn{\alpha}{alpha} difficulty parameters. This can either be
    a scalar or a column vector with dimension equal to the
    number of items.   If this takes a scalar value, then that
    value will serve as the starting value for all of the alphas.
     The default value of NA will set the starting values based on
     a series of probit regressions that condition on the starting
     values of theta. }

  \item{beta.start}{ The starting values for the
    \eqn{\beta}{beta} discrimination parameters. This can either
    be a scalar or a column vector with dimension equal to the
    number of items.   If this takes a scalar value, then that
    value will serve as the starting value for all of the betas. 
    The default value of NA will set the starting values based on a
    series of probit regressions that condition on the starting values
    of theta. }

  \item{tau2.start}{ The starting values for the evolution variances
  (the variance of the random walk increments for the ability
  parameters / ideal points. Order corresponds to the rows of
  \code{datamatrix}.}

  \item{a0}{ A vector containing the prior mean of each of the
  difficulty parameters \eqn{\alpha}{alpha}. Should have as many
  elements as items / roll calls. Order corresponds to the columns of
  \code{datamatrix}. If a scalar is passed it is assumed that all
  elements of \code{a0} are equal to the scalar. }

  \item{A0}{ A vector containing the prior precision (inverse
  variance) of each of the difficulty parameters \eqn{\alpha}{alpha}.
  Should have as many elements as items / roll calls. Order
  corresponds to the columns of \code{datamatrix}. If a scalar is
  passed it is assumed that all elements of \code{A0} are equal to the
  scalar. }

  \item{b0}{ A vector containing the prior mean of each of the
  discrimination parameters \eqn{\beta}{beta}. Should have as many
  elements as items / roll calls. Order corresponds to the columns of
  \code{datamatrix}. If a scalar is passed it is assumed that all
  elements of \code{b0} are equal to the scalar.  }

  \item{B0}{ A vector containing the prior precision (inverse
  variance) of each of the discrimination parameters
  \eqn{\beta}{beta}. Should have as many elements as items / roll
  calls. Order corresponds to the columns of \code{datamatrix}. If a
  scalar is passed it is assumed that all elements of \code{B0} are
  equal to the scalar. }

  \item{c0}{\eqn{c_{0/2}}{c0/2} is the shape parameter for the inverse
    Gamma prior on \eqn{\tau^2}{tau^2} (the variance of the random
    walk increments). The amount of information in the inverse Gamma
    prior is something like that from \eqn{c_0}{c0}
    pseudo-observations. \code{c0} can be either a vector with an
    element for each subject or a scalar. If \code{c0}
    is negative then \eqn{\tau^2}{tau^2} is not estimated-- the values
    in \code{tau2.start} are used throughout the sampling.}

  \item{d0}{ \eqn{d_{0/2}}{d0/2} is the scale parameter for the inverse
    Gamma prior on \eqn{\tau^2}{tau^2} (the variance of the random
    walk increments). In constructing the inverse Gamma prior,
    \eqn{d_0}{d0} acts like the sum of squared errors from the
    \eqn{c_0}{c0} pseudo-observations. \code{d0} can be either a
    vector with an element for each subject or a scalar. If \code{d0}
    is negative then \eqn{\tau^2}{tau^2} is not estimated-- the values
    in \code{tau2.start} are used throughout the sampling. }

  \item{e0}{ A vector containing the prior mean of the initial ability
  parameter / ideal point for each subject. Should have as many
  elements as subjects. Order corresponds to the rows of \code{datamatrix}. If a
  scalar is passed it is assumed that all elements of \code{e0} are
  equal to the scalar.}

  \item{E0}{ A vector containing the prior variance of the initial
  ability parameter / ideal point for each subject. Should have as
  many elements as subjects. Order corresponds to the rows of
  \code{datamatrix}. If a scalar is passed it is assumed that all
  elements of \code{E0} are equal to the scalar.}

  \item{store.ability}{ A switch that determines whether or not to
    store the ability parameters for posterior analysis. 
    \emph{NOTE}: In situations with many individuals storing the ability
      parameters takes an enormous amount of memory, so
      \code{store.ability} should only be \code{TRUE} if the chain is thinned
      heavily, or for applications with a small number of individuals.
    By default, the item parameters are stored. }

  \item{store.item}{A switch that determines whether or not to
    store the item parameters for posterior analysis. 
    \emph{NOTE}: In situations with many items storing the item
      parameters takes an enormous amount of memory, so
      \code{store.item} should only be \code{FALSE} if the chain is thinned
      heavily, or for applications with a small number of items.
    By default, the item parameters are not stored.}

  \item{\dots}{further arguments to be passed  }

}

\details{
    \code{MCMCdynamicIRT1d} simulates from the posterior distribution using
   the algorithm of Martin and Quinn (2002). The simulation proper is done in
  compiled C++ code to maximize efficiency.  Please consult the
  coda documentation for a comprehensive list of functions that
  can be used to analyze the posterior sample.
  
  The model takes the following form. We assume that each subject has
  an subject ability (ideal point) denoted
  \eqn{\theta_{j,t}}{theta_jt} (where \eqn{j}{j} indexes subjects and
  \eqn{t}{t} indexes time periods) and that each item has a difficulty
  parameter \eqn{\alpha_i}{alpha_i} and discrimination parameter
  \eqn{\beta_i}{beta_i}. The observed choice by subject \eqn{j}{j} on
  item \eqn{i}{i} is the observed data matrix which is \eqn{(I \times
  J)}{(I * J)}. We assume that the choice is dictated by an unobserved
  utility: \deqn{z_{i,j,t} = -\alpha_i + \beta_i \theta_{j,t} +
  \varepsilon_{i,j,t}}{z_ijt = -alpha_i + beta_i*theta_jt + epsilon_ijt}
  Where the disturbances are assumed to be distributed standard Normal. The
  parameters of interest are the subject abilities (ideal points) and
  the item parameters.

  We assume the following priors.  For the subject abilities (ideal points):
  \deqn{\theta_{j,t} \sim \mathcal{N}(\theta_{j,t-1}, \tau^2_j)}{theta_jt ~ N(theta_j(t-1), tau^2)} with \deqn{\theta_{j,0} \sim \mathcal{N}(e0, E0)}{theta_j0 ~ N(e0, E0).}   
  
  The evolution variance has the following prior: \deqn{\tau^2_j \sim \mathcal{IG}(c0/2, d0/2)}{tau^2_j ~ IG(c0/2, d0/2).} 

  For the item parameters in the standard model, the prior is:
  \deqn{\alpha_i \sim \mathcal{N}(a0, A0^{-1})}{alpha_i ~ N(a0, A0^{-1})} and
  \deqn{\beta_i \sim \mathcal{N}(b0, B0^{-1})}{beta_i ~ N(b0, B0^{-1}).}

  The model is identified by the proper priors on the item parameters
  and constraints placed on the ability parameters.
  
  As is the case with all measurement models, make sure that you have plenty
  of free memory, especially when storing the item parameters.  
}

\value{ An mcmc object that contains the posterior sample. This object
can be summarized by functions provided by the coda package. }

\references{ Andrew D. Martin and Kevin M. Quinn. 2002. "Dynamic Ideal
Point Estimation via Markov Chain Monte Carlo for the U.S. Supreme
Court, 1953-1999." \emph{Political Analysis.} 10: 134-153.

Andrew D. Martin, Kevin M. Quinn, and Jong Hee Park. 2011.
 ``MCMCpack: Markov Chain Monte Carlo in R.'',
 \emph{Journal of Statistical Software}. 42(9): 1-21.
 \url{http://www.jstatsoft.org/v42/i09/}.
}

\author{Kevin M. Quinn }

\seealso{\code{\link[coda]{plot.mcmc}},\code{\link[coda]{summary.mcmc}},
\code{\link[MCMCpack]{MCMCirt1d}} }


\examples{
  \dontrun{
	data(Rehnquist)

	## assign starting values
	theta.start <- rep(0, 9)
	theta.start[2] <- -3 ## Stevens
	theta.start[7] <- 2  ## Thomas

	out <- MCMCdynamicIRT1d(t(Rehnquist[,1:9]),
	                        item.time.map=Rehnquist$time,
	                        theta.start=theta.start,
	                        mcmc=50000, burnin=20000, thin=5,
	                        verbose=500, tau2.start=rep(0.1, 9),
	                        e0=0, E0=1,
	                        a0=0, A0=1,
	                        b0=0, B0=1, c0=-1, d0=-1,
	                        store.item=FALSE, 
	                        theta.constraints=list(Stevens="-", Thomas="+"))

	summary(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }

