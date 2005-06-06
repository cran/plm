\name{plag}
\alias{plag}
\title{Lag operator for panel data}
\description{Computes the lag of variables that have a time and an
  individual dimension}
\usage{
plag(x,order=1)
}
\arguments{
  \item{x}{a numeric vector or a numeric matrix for which the
    lag should be computed,}
  \item{order}{the order of the lag to be computed.}
}

\value{
  a vector (or matrix) of the same dimension as \code{x} containing the
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
