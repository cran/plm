\name{pFtest}
\alias{pFtest}
\alias{pFtest.default}
\alias{pFtest.plm}
\alias{pFtest.plms}
\title{F test for effects}
\description{
  Test of individuals and/or time effects based on the comparison of the
  \code{within} and the \code{pooling} model.
}
\usage{
pFtest(x, ...)
\method{pFtest}{plm}(x,z,data=NULL, ...)
\method{pFtest}{plms}(x,...)
}

\arguments{
  \item{x}{an object of class \code{"plms"} or \code{"plm"}}
  \item{z}{an object of class \code{"plm"}}
  \item{data}{}
  \item{...}{further arguments}
}
\value{
  An object of class \code{"htest"}.
}

\details{
 The argument of this function is whether a \code{plms} object or two
 \code{plm} objects, the first being a within model, the second a pooling
 model. The effects tested are whether individual, time or twoways effects
 depending on the effects introduced in the model.
}
\seealso{
  \code{\link{plmtest}} for Lagrange multipliers tests of individuals
  and/or time effects.
}
\examples{
library(Ecdat)
data(Grunfeld)
Grunfeld <-pdata.frame(Grunfeld,"firm","year")
gi <- plm(inv ~ value + capital, data=Grunfeld)
gt <- plm(inv ~ value + capital, data=Grunfeld,effect="time")
gd <- plm(inv ~ value + capital, data=Grunfeld,effect="twoways")
pFtest(gi)
pFtest(gi$within,gi$pooling)
pFtest(gt)
pFtest(gd)
}
\keyword{htest}
