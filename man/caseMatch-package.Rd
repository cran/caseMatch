\name{caseMatch-package}
\alias{caseMatch-package}
\alias{caseMatch}
\docType{package}
\title{
  A package for using matching to select cases from a quantitative data set for further qualitative analysis.
}
\description{
  This package uses statistical matching to identify "most similar" cases in a quantitative data set for subsequent qualitative analysis.  Unlike existing matching packages, this package intended to meet some specific needs of analysts using matching for case studies.
}
\details{
  Use the \code{case.match} function.
}
\author{
Maintainer: Rich Nielsen <rnielsen@mit.edu>
}
\references{
Nielsen, Richard A. 2016. "Case Selection Via Matching," Sociological Methods and Research.  http://journals.sagepub.com/doi/abs/10.1177/0049124114547054
}
% Optionally other standard keywords, one per line,
% from the file KEYWORDS in the R documentation.
\keyword{ package }
\seealso{ \code{\link{case.match}}}
\examples{

data(EU)
mvars <- c("socialist","rgdpc","FHc","FHp","trade")
dropvars <- c("countryname","population")

## In this example, I subset to the first 40 obs. to cut run-time
out <- case.match(data=EU[1:40,], id.var="countryname", leaveout.vars=dropvars,
             distance="mahalanobis", case.N=2, greedy.match="pareto", 
             number.of.matches.to.return=10,
             treatment.var="eu", max.variance=TRUE)
out$cases

}
