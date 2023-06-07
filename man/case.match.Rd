\name{case.match}
\alias{case.match}
\title{Uses matching methods to select cases for qualitative analysis}
\description{
Uses matching methods to select cases for qualitative analysis
}
\usage{
case.match(data, id.var, case.N = 2, distance = "mahalanobis", 
    design.type = "most similar", match.case = NULL, 
    greedy.match="pareto", number.of.matches.to.return = 1, 
    treatment.var = NULL, outcome.var= NULL, leaveout.vars = NULL, 
    max.variance = FALSE,  max.variance.outcome=FALSE,
    variance.tolerance = 0.1, max.spread = FALSE, 
    max.spread.outcome=FALSE, varweights = NULL)
}
\arguments{
  \item{data}{A data frame.}
  \item{id.var}{A string variable that uniquely identifies cases within the data}
  \item{case.N}{The number of cases to choose.  Must be 1 or more.}
  \item{distance}{The distance metric, specified as a string.  Options are "mahalanobis", "euclidean", or "standardized", where "standardized" means that variables are standardized by their standard deviations.}
  \item{design.type}{Should the algorithm pick cases that are most similar or most different?  Specify either "most similar" or "most different" as a string.}
  \item{match.case}{If specified, this is the value of \code{id.var} of a specific case to match.}
  \item{number.of.matches.to.return}{How many matches to return.}
  \item{greedy.match}{Specifies which matches to return. Options are "pareto", "greedy", and "all". "all" keeps all matches. "pareto" matches eliminate 'redundant' matches where both units have better available matches. "greedy" matches keeps only the top matches in the data, but does eliminates best matches for some units since it uses a without replacement algorithm.}
  \item{treatment.var}{The name of the treatment variable, specified as a string.}
  \item{outcome.var}{The name of the outcome variable, specified as a string.}
  \item{leaveout.vars}{A vector of variables to not include in the matching.}
  \item{max.variance}{Should the cases be selected to maximize variance on \code{treatment.var}?  If cases should be in opposite treatment conditions, specify \code{max.variance=TRUE}.}
  \item{max.variance.outcome}{Should the cases be selected to maximize variance on \code{outcome.var}?  If cases should have opposite outcomes, specify \code{max.variance.outcome=TRUE}.}
  \item{variance.tolerance}{The proportion of cases to consider if \code{max.variance} is specified but there are too few cases that maximize the variance of \code{treatment.var}.}
  \item{max.spread}{Should the cases be selected to maximize "spread" on the treatment variable, meaning that cases are selected to be have evenly values from the min of \code{treatment.var} to the max?}
  \item{max.spread.outcome}{Should the cases be selected to maximize "spread" on the outcome variable, meaning that cases are selected to be have evenly values from the min of \code{outcome.var} to the max?}
  \item{varweights}{An optional vector of variable weights.  It must line up with the columns of the data after \code{id.var} and \code{leaveout.vars} are removed.  Optionally, element names can be included for \code{varweights} --- if so, the function checks that the names are identical to (and line up with) the names of the matching variables.  It will throw an error if they do not.}
}
\details{
\code{case.match} uses statistical matching to select cases in a quantitative data set for subsequent qualitative analysis in "most similar" and "most different" research designs.
}
\value{
\code{case.match} returns a named list with the following elements:
  \item{cases}{A table of the matched cases.}
  \item{case.distances}{A list of the distances between matched cases.}
}
\references{
Nielsen, Richard. 2016. "Case Selection via Matching," Sociological Methods and Research, 45 (3): 569-597.  http://www.mit.edu/~rnielsen/Case%20Selection%20via%20Matching.pdf
}
\author{
Rich Nielsen
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{

data(EU)
mvars <- c("socialist","rgdpc","FHc","FHp","trade")
dropvars <- c("countryname","population")

## In this example, I subset to the first 40 obs. to cut run-time
out <- case.match(data=EU[1:40,], id.var="countryname", leaveout.vars=dropvars,
             distance="mahalanobis", case.N=2, 
             number.of.matches.to.return=10,
             treatment.var="eu", max.variance=TRUE)
out$cases

\dontrun{
## All cases:
## Find the best matches of EU to non-EU countries
out <- case.match(data=EU, id.var="countryname", leaveout.vars=dropvars,
             distance="mahalanobis", case.N=2, 
             number.of.matches.to.return=10,
             treatment.var="eu", max.variance=TRUE)
out$cases

## Find the best matches while downweighting political variables
myvarweights <- c(1,1,.1,.1,.1)
names(myvarweights) <- c("rgdpc","trade","FHp","FHc","socialist")
myvarweights
(case.match(data=EU, id.var="countryname", leaveout.vars=dropvars,
             distance="mahalanobis", case.N=2, 
             number.of.matches.to.return=10, treatment.var="eu",
             max.variance=TRUE,varweights=myvarweights))$cases

## Find the best non-EU matches for Germany
tabGer <- case.match(data=EU, match.case="German Federal Republic", 
             id.var="countryname",leaveout.vars=dropvars,
             distance="mahalanobis", case.N=2, 
             number.of.matches.to.return=10,max.variance=TRUE,
             treatment.var="eu")
}



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~case selection }
\keyword{ ~matching }
