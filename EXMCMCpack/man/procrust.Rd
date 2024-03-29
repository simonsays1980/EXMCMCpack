\name{procrustes}
\alias{procrustes}


\title{Procrustes Transformation}
\description{
This function performs a Procrustes transformation on a matrix \code{X} to
minimize the squared distance between \code{X} and another matrix \code{Xstar}.
}
\usage{
procrustes(X, Xstar, translation=FALSE, dilation=FALSE)
}

\arguments{
  \item{X}{The matrix to be transformed.}

  \item{Xstar}{The target matrix.}
  
  \item{translation}{logical value indicating whether \code{X} should be
    translated.}

  \item{dilation}{logical value indicating whether \code{X} should be
    dilated.}   
}

\value{
  A list containing: \code{X.new} the matrix that is the Procrustes
  transformed version of \code{X}, \code{R} the rotation matrix,
  \code{tt} the translation vector, and \code{s} the scale factor. 
}

\details{
  \code{R}, \code{tt}, and \code{s} are chosen so that:

  \deqn{s X R + 1 tt' \approx X^*}{s X R + 1 tt' approximately Xstar}

  \code{X.new} is given by:

  \deqn{X_{new} = s X R + 1 tt'}{X.new = s X R + 1 tt'}
  }


\references{
  Borg and Groenen. 1997. \emph{Modern Multidimensional Scaling}. New
  York: Springer. pp. 340-342.
  }


\seealso{\code{\link{MCMCirtKd}}}
\keyword{manip}

