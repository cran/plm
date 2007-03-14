\name{FE}
\alias{FE}
\alias{FE.default}
\alias{FE.plm}
\alias{FE.plms}
\alias{print.FE}
\alias{summary.FE}
\title{Extract the fixed effects}
\description{
  This function extracts the fixed effects from a \code{plm} or a
  \code{plms} object.
}

\usage{
FE(x)
\method{FE}{plm}(x)
\method{FE}{plms}(x)
\method{print}{FE}(x,digits=5,...)
\method{summary}{FE}(object,...)
}

\arguments{
  \item{x,object}{an object of class \code{"plms"} or \code{"plm"}, an object
    of class \code{FE} for the \code{print} and the \code{summary} method,}
  \item{digits}{digits,}
  \item{...}{further arguments.}
}

\value{
  An object of class \code{FE}. It is a numeric vector containing the
  fixed effects with two attributes :  \code{se} which contains the standard
  errors and \code{intercept} which is the overall intercept.
}

\details{
  The element \code{FE} of a \code{plm} object corresponding to a within
  model is a numerical vector which length is the total number of
  observations. This function returns a vector which length is the number
  of individuals. The summary method prints the effects in deviation
  from  the overall intercept, the standard error and the t--values.
}
\seealso{
  \code{\link{plm}}
}

\author{Yves Croissant}

\examples{
library(Ecdat)
data(Grunfeld)
pdata.frame(Grunfeld,"firm","year")
gi <- plm(inv ~ value + capital, data=Grunfeld)
FE(gi)
FE(gi$within)
summary(FE(gi))
}
\keyword{attribute}
