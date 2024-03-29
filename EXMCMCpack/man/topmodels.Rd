\name{topmodels}
\alias{topmodels}

\title{Shows an ordered list of the most frequently visited models sampled during quantile regression stochastic search variable selection (QR-SSVS).}

\description{Given output from quantile regression stochastic search variable selection, this function returns a table of the 'best' models together with their associated empirical posterior probability.}

\usage{
topmodels(qrssvs, nmodels=5, abbreviate=FALSE, minlength=3)
}

\arguments{

  \item{qrssvs}{An object of class \code{qrssvs}. Typically this will be the \code{gamma} component of the list returned by \code{SSVSquantreg}.}

  \item{nmodels}{The number of models to tabulate.}

  \item{abbreviate}{Logical: should the names of the predictors be abbreviated?}

  \item{minlength}{If \code{abbreviate} is set to \code{TRUE}, the minimum length of the abbreviations. }
}

\value{A table with the models and their associated posterior probability. The models are arranged in descending order of probability.}

\author{Craig Reed}

\examples{
\dontrun{
set.seed(1)
epsilon<-rnorm(100)
set.seed(2)
x<-matrix(rnorm(1000),100,10)
y<-x[,1]+x[,10]+epsilon
qrssvs<-SSVSquantreg(y~x)
topmodels(qrssvs$gamma)
}
}

\keyword{models}

\seealso{
\code{\link[MCMCpack]{SSVSquantreg}}}
