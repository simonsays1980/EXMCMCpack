\name{make.breaklist}
\alias{make.breaklist}
\title{Vector of break numbers}
\description{This function generates a vector of break numbers using the output of \code{testpanelSubjectBreak}. 
The function performs a pairwise comparison of models using Bayes Factors. }

\usage{
   make.breaklist(BF, threshold=3)
}

\arguments{

\item{BF}{output of \code{testpanelSubjectBreak}.}

\item{threshold}{The Bayes Factor threshold to pick the best model. 
		 If a Bayes factor of two models is smaller than \code{threshold},  
		 the model with a smaller number of break is chosen to avoid the over-identification problem.
		 Users can change threshold into any positive number. 
		 The default value of 3 is chosen as it indicates the existence of "substantial evidence" in favor of the model in the numerator
		 according to Jeffreys' scale.}

}
\references{
   Jong Hee Park, 2011. ``A Unified Method for Dynamic and Cross-Sectional Heterogeneity: 
   Introducing Hidden Markov Panel Models." Working Paper.

   Harold Jeffreys, 1961. The Theory of Probability. Oxford University Press.    
}

\seealso{\code{\link{testpanelSubjectBreak}}}
