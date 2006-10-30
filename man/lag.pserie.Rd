\name{lag.pserie}
\alias{lag}
\alias{lag.pserie}
\title{Lag operator for panel data}
\description{Computes the lag of variables that have a time and an
  individual dimension.}
\usage{
lag.pserie(x,k=1,...)
}
\arguments{
  \item{x}{a vector of class \code{pserie},}
  \item{k}{the order of the lag to be computed,}
  \item{...}{further arguments.}
}

\value{
  a vector of the same length as \code{x} containing the
  lag values. \code{NA} values are inserted at the correct positions.}

\seealso{\code{\link{diff.pserie}}}
\examples{
library(Ecdat)
data(Produc)
Produc=pdata.frame(Produc,"state","year")
lag(Produc$pc)
lag(Produc$pc,k=2)
}
\keyword{arith}
