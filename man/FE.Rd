\name{FE}
\alias{FE}
\alias{FE.default}
\alias{FE.plm}
\alias{FE.plms}
\title{Extract the fixed effects}
\description{
  This function extracts the fixed effects from a \code{plm} or a
  \code{plms} object.
}
\usage{
FE(x)
\method{FE}{plm}(x)
\method{FE}{plms}(x)
}

\arguments{
  \item{x}{an object of class \code{"plms"} or \code{"plm"}}
}
\value{
  A numeric vector containing the fixed effects.
}

\details{
  The element \code{FE} of a \code{plm} object corresponding to a within
  model is a numerical vector which length is the total number of
  observations. This function returns a vector which length is the number
  of individuals.
}
\seealso{\code{\link{plm}}}

\examples{
library(Ecdat)
data(Grunfeld)
pdata.frame(Grunfeld,"firm","year")
gi <- plm(inv ~ value + capital, data=Grunfeld)
FE(gi)
FE(gi$within)
}
\keyword{attribute}
