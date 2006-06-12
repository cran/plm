\name{pmean}
\alias{pmean}
\title{Conditional mean of a variable}
\description{Computes the mean of a variable conditionnaly to the
  time or the individual index}
\usage{
pmean(x,cond=NULL,effect="individual",data=NULL)
}
\arguments{
  \item{x}{a numeric vector or a numeric matrix for which the
    conditional means should be computed,}
  \item{cond}{a vector containing the conditional variable,}
  \item{effect}{the mean is conditionnal to the individual 
    or the time index whether effect equals \code{"individual"} or \code{"time"},}
  \item{data}{a \code{pdata.frame}.}
}

\value{
  a vector (or matrix) of the same dimension as \code{x} containing the
  conditional means. An attribute \code{"cm"} contains a vector (or
  matrix) containing the conditional means with one line for each value
  of the conditionnal variable.}

\details{
  If \code{x} is a variable of a \code{pdata.frame}, \code{pmean(x)}
  finds the \code{"data"} attribute of \code{x}, extract the individual
  index from this \code{pdata.frame} and computes the individual
  means. If \code{effect} is \code{"time"}, the time means are
  computed.
  The conditional variable may also be explicitely indicated using the
  \code{cond} (a vector of the same length as \code{x}) or the
  \code{data} (a \code{pdata.frame}) argument. This is the only way to
  use this function if \code{x} is a matrix.
}

\examples{
library(Ecdat)
data(Produc)
Produc=pdata.frame(Produc,state,year)
pmean(Produc$pc)
form <- pmean(log(gsp))~pmean(log(pc))
lm(form,data=Produc)
}
\keyword{arith}
