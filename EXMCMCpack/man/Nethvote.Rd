\name{Nethvote}
\alias{Nethvote}
\docType{data}
\title{Dutch Voting Behavior in 1989}
\description{
Dutch Voting Behavior in 1989.
}
\usage{data(Nethvote)}
\format{
  A data frame with 1754 observations and 11 variables from the 1989
  Dutch Parliamentary Election Study
  (Anker and Oppenhuis, 1993). Each observation is a survey respondent.
  These data are a subset of one of five multiply imputed datasets used in
  Quinn and Martin (2002). For more information see Quinn and Martin
  (2002). 
  \describe{
    \item{vote}{A factor giving the self-reported vote choice of each
      respondent. The levels are CDA (Christen Democratisch Appel), D66
      (Democraten 66), Pvda (Partij van de Arbeid), and VVD (Volkspartij
      voor Vrijheid en Democratie).}
    \item{distD66}{A numeric variable giving the squared ideological distance
      between the respondent and the D66. Larger values indicate ideological
      dissimilarity between the respondent and the party.}
    \item{distPvdA}{A numeric variable giving the squared ideological distance
      between the respondent and the PvdA. Larger values indicate ideological
      dissimilarity between the respondent and the party.}
    \item{distVVD}{A numeric variable giving the squared ideological distance
      between the respondent and the VVD. Larger values indicate ideological
      dissimilarity between the respondent and the party.}
    \item{distCDA}{A numeric variable giving the squared ideological distance
      between the respondent and the CDA. Larger values indicate ideological
      dissimilarity between the respondent and the party.}
    \item{relig}{An indicator variable equal to 0 if the respondent is
      not religious and 1 if the respondent is religious.}
    \item{class}{Social class of respondent. 0 is the lowest social class, 4 is
      the highest social class.}
    \item{income}{Income of respondent. 0 is lowest and 6 is highest.}
    \item{educ}{Education of respondent. 0 is lowest and 4 is highest.}
    \item{age}{Age category of respondent. 0 is lowest and 12 is
      highest.}
    \item{urban}{Indicator variable equal to 0 if the respondent is not
      a resident of an urban area and 1 if the respondent is a resident
      of an urban area.}
  }
}
\source{
  H. Anker and E.V. Oppenhuis. 1993. ``Dutch Parliamentary Election
  Study.'' (computer file). Dutch Electoral Research Foundation and
  Netherlands Central Bureau of Statistics, Amsterdam.
}
\references{
  Kevin M. Quinn and Andrew D. Martin. 2002. ``An Integrated Computational Model
  of Multiparty Electoral Competition.'' \emph{Statistical Science}. 17:
  405-419.   
}
\keyword{datasets}
