\name{Rehnquist}
\alias{Rehnquist}
\title{
  U.S. Supreme Court Vote Matrix, Rehnquist Court (1994-2004)
}
\description{
  This dataframe contains a matrix of votes cast by U.S. Supreme
  Court justices by all cases in the 1994-2004 terms.
}
\usage{
data(SupremeCourt)
}
\format{
  The dataframe has contains data for justices Rehnquist, Stevens,
  O'Connor, Scalia, Kennedy, Souter, Thomas, Ginsburg, and Breyer
  for the 1994-2004 terms of the U.S. Supreme Court.  The dataframe
  also contains the term of the case, and a time variable that counts
  from term 1 to 11.  The votes are coded liberal (1)
  and conservative (0) using the protocol of Spaeth (2003). 
  The unit of analysis
  is the case citation (ANALU=0).  We are concerned with formally
  decided cases issued with written opinions, after full oral
  argument and cases decided by an equally divided vote
  (DECTYPE=1,5,6,7).
}

\source{
  Harold J. Spaeth. 2005. \emph{Original United States Supreme Court Database: 
    1953-2004 Terms.} \url{http://www.as.uky.edu/polisci/ulmerproject/sctdata.htm}.
}

\keyword{datasets}
