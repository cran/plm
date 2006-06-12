\name{phtest}
\alias{phtest}
\alias{phtest.default}
\alias{phtest.plm}
\alias{phtest.plms}
\title{Hausman's test for panel models}
\description{
  Test of consistency for panel models
}
\usage{
phtest(x, ...)
\method{phtest}{plm}(x,x2,data=NULL, ...)
\method{phtest}{plms}(x,...)
}

\arguments{
  \item{x}{an object of class \code{"plms"} or \code{"plm"}}
  \item{x2}{an object of class \code{"plm"}}
  \item{data}{the name of the data}
  \item{...}{further arguments}
}

\value{
  An object of class \code{"htest"}.
}

\details{
  The Hausman test is based on the difference of the vectors of
  coefficients of two different models. If the argument of \code{phtest}
  is an object of class \code{plms}, the Hausman test compares the
  within and the random effects model.
}
\references{
  Hausman, G. (1978), Specification tests in econometrics,
  \emph{Econometrica}, 46, pp.1251--1271.
  }

\examples{
library(Ecdat)
data(Gasoline)
Gasoline <- pdata.frame(Gasoline,country,year)
form <- lgaspcar~lincomep+lrpmg+lcarpcap
gas <-  plm(form,data=Gasoline)
phtest(gas)
phtest(gas$between, gas$random)
}
\keyword{htest}
