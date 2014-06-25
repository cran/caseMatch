\name{EU}
\alias{EU}
\docType{data}
\title{
Cross-national data for 189 countries.
}
\description{
A cross-national data set including economic and political variables for 189 countries, averaged from 1980-1992.
}
\usage{data(EU)}
\format{
  A data frame with 185 observations on the following 13 variables.
  \describe{
    \item{\code{countryname}}{The name of the country}
    \item{\code{population}}{Country population from Gleditsch.}
    \item{\code{rgdpc}}{GDP per capita from Gleditsch.}
    \item{\code{trade}}{Trade from Gleditsch.}
    \item{\code{FHp}}{Freedom House political rights.}
    \item{\code{FHc}}{Freedom House civil rights.}
    \item{\code{socialist}}{An indicator for countries that were socialist during the Cold War.}
    \item{\code{eu}}{An indicator for EU members.}
  }
}
\details{
  A cross-national data set including economic and political variables for 189 countries, averaged from 1980-1992.  Data are collected by Gleditsch and Freedom House.
}
\source{
Gleditsch, Kristian Skrede. (2004) Expanded Trade and GDP Data,
Version 4.0. http://privatewww.essex.ac.uk/~ksg/exptradegdp.html

http://www.freedomhouse.org/report-types/freedom-world
}
\references{
Nielsen, Richard A. Forthcoming. "Case Selection Via Matching," Sociological Methods and Research.  http://www.mit.edu/~rnielsen/Case Selection via Matching.pdf

}
\examples{
data(EU)
}
\keyword{datasets}
