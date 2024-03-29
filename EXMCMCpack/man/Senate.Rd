\name{Senate}
\alias{Senate}

\title{
  106th U.S. Senate Roll Call Vote Matrix
}
\description{
  This dataframe contains a matrix of votes cast by U.S. Senators
  in the 106th Congress.
}
\usage{
data(Senate)
}
\format{
  The dataframe contains roll call data for all Senators in the 106th 
  Senate.  The first column (id) is the ICPSR member ID number, the
  second column (statecode) is the ICPSR state code, the third column (party)
  is the member's state name, and the fourth column (member) is the member's
  name.  This is followed by all roll call votes (including unanimous ones)
  in the 106th.   Nay votes are coded 0, yea votes are coded 1,
  and NAs are missing votes.
}

\source{
  Keith Poole. 2005. \emph{106th Roll Call Vote Data}.
  \url{http://voteview.uh.edu/}.
}

\keyword{datasets}
