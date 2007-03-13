\name{diff.pserie}
\alias{diff.pserie}
\alias{diff}
\title{Difference operator for panel data}
\description{Computes the first difference of variables that have a time and
  an individual dimension.}
\usage{
diff.pserie(x,lag=0, ...)
}
\arguments{
  \item{x}{a numeric vector for which the
    difference should be computed,}
  \item{lag}{the lag for the first difference,}
  \item{...}{further arguments.}
}

\value{
  a vector of the same dimension as \code{x} containing the
  lag values. \code{NA} values are inserted at the correct positions.}


\seealso{\code{\link{lag.pserie}}}

\author{Yves Croissant}
\examples{
# Anderson and Hsiao estimator
library(Ecdat)
data(Grunfeld)
pdata.frame(Grunfeld,"firm","year")
form <- diff(inv)~diff(inv,1)+diff(value)+diff(capital)
hs <- plm(form,data=Grunfeld,instruments=~lag(inv,2),endog=~diff(inv,1),model="pooling")
}
\keyword{arith}
