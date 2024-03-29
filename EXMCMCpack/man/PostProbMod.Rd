\name{PostProbMod}
\alias{PostProbMod}


\title{Calculate Posterior Probability of Model}
\description{
This function takes an object of class \code{BayesFactor} and calculates
the posterior probability that each model under study is correct given
that one of the models under study is correct.

}
\usage{
PostProbMod(BF, prior.probs=1)
}

\arguments{
  \item{BF}{An object of class \code{BayesFactor}.}

  \item{prior.probs}{The prior probabilities that each model is
    correct. Can be either a scalar or array. Must be positive. If the
    sum of the prior probabilities is not equal to 1 prior.probs will be
    normalized so that it does sum to unity.}
}

\value{
  An array holding the posterior probabilities that each model under
  study is correct given that one of the models under study is correct.
}

\examples{
\dontrun{
data(birthwt)

post1 <- MCMCregress(bwt~age+lwt+as.factor(race) + smoke + ht,
                     data=birthwt, b0=c(2700, 0, 0, -500, -500,
                                        -500, -500),
                     B0=c(1e-6, .01, .01, 1.6e-5, 1.6e-5, 1.6e-5,
                          1.6e-5), c0=10, d0=4500000,
                     marginal.likelihood="Chib95", mcmc=10000)
 
post2 <- MCMCregress(bwt~age+lwt+as.factor(race) + smoke,
                     data=birthwt, b0=c(2700, 0, 0, -500, -500,
                                        -500),
                     B0=c(1e-6, .01, .01, 1.6e-5, 1.6e-5, 1.6e-5),
                     c0=10, d0=4500000,
                     marginal.likelihood="Chib95", mcmc=10000)

post3 <- MCMCregress(bwt~as.factor(race) + smoke + ht,
                     data=birthwt, b0=c(2700, -500, -500,
                                        -500, -500),
                     B0=c(1e-6, 1.6e-5, 1.6e-5, 1.6e-5,
                          1.6e-5), c0=10, d0=4500000,
                     marginal.likelihood="Chib95", mcmc=10000)

BF <- BayesFactor(post1, post2, post3)
mod.probs <- PostProbMod(BF)
print(mod.probs)
}
}

\concept{Bayes factor}
\concept{model comparison}

\seealso{\code{\link{MCMCregress}}}
\keyword{models}

