\name{pdiff}
\alias{pdiff}
\title{Difference operator for panel data}
\description{Computes the first difference of variables that have a temporal and
  an individual dimension.}
\usage{
pdiff(x,lag=0)
}
\arguments{
  \item{x}{a numeric vector or a numeric matrix for which the
    difference should be computed,}
  \item{lag}{the lag for the first difference.}
}

\value{
  a vector (or matrix) of the same dimensions as \code{x} containing the
  lag values. \code{NA} values are inserted at the correct positions}


\examples{
library(Ecdat)
data(Produc)
Produc=pdata.frame(Produc,state,year)
pmean(Produc$pc)
form <- pmean(log(gsp))~pmean(log(pc))
lm(form,data=Produc)
}
\keyword{arith}
