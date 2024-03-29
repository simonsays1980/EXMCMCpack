\name{BayesFactor}
\alias{BayesFactor}
\alias{is.BayesFactor}

\title{Create an object of class BayesFactor from MCMCpack output}
\description{
This function creates an object of class \code{BayesFactor} from MCMCpack
output.

}
\usage{
BayesFactor(...)
is.BayesFactor(BF)
}

\arguments{
  \item{...}{MCMCpack output objects. These have to be of class
    \code{mcmc} and have a \code{logmarglike} attribute. In what
    follows, we let \code{M} denote the total number of models to be
    compared.} 

  \item{BF}{An object to be checked for membership in class
    \code{BayesFactor}.} 
}

\value{
  An object of class \code{BayesFactor}. A \code{BayesFactor} object has
  four attributes. They are: \code{BF.mat} an \eqn{M \times M}{M by M}
  matrix in which element \eqn{i,j}{i,j} contains the Bayes factor for
  model \eqn{i}{i} relative to model \eqn{j}{j}; \code{BF.log.mat} an
  \eqn{M \times M}{M by M} matrix in which element \eqn{i,j}{i,j}
  contains the natural log of the Bayes factor for model \eqn{i}{i}
  relative to model \eqn{j}{j}; \code{BF.logmarglike} an \eqn{M}{M}
  vector containing the log marginal likelihoods for models 1 through
  \eqn{M}{M}; and \code{BF.call} an \eqn{M}{M} element list containing
  the calls used to fit  models 1 through \eqn{M}{M}.   
}

\examples{
\dontrun{
data(birthwt)

model1 <- MCMCregress(bwt~age+lwt+as.factor(race) + smoke + ht,
                     data=birthwt, b0=c(2700, 0, 0, -500, -500,
                                        -500, -500),
                     B0=c(1e-6, .01, .01, 1.6e-5, 1.6e-5, 1.6e-5,
                          1.6e-5), c0=10, d0=4500000,
                     marginal.likelihood="Chib95", mcmc=10000)
 
model2 <- MCMCregress(bwt~age+lwt+as.factor(race) + smoke,
                     data=birthwt, b0=c(2700, 0, 0, -500, -500,
                                        -500),
                     B0=c(1e-6, .01, .01, 1.6e-5, 1.6e-5, 1.6e-5),
                     c0=10, d0=4500000,
                     marginal.likelihood="Chib95", mcmc=10000)

model3 <- MCMCregress(bwt~as.factor(race) + smoke + ht,
                     data=birthwt, b0=c(2700, -500, -500,
                                        -500, -500),
                     B0=c(1e-6, 1.6e-5, 1.6e-5, 1.6e-5,
                          1.6e-5), c0=10, d0=4500000,
                     marginal.likelihood="Chib95", mcmc=10000)

BF <- BayesFactor(model1, model2, model3)
print(BF)

}
}

\concept{Bayes factor}
\concept{model comparison}

\seealso{\code{\link{MCMCregress}}}
\keyword{models}
