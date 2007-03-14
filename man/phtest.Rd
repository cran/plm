\name{phtest}
\alias{phtest}
\alias{phtest.default}
\alias{phtest.panelmodel}
\alias{phtest.plms}
\title{Hausman test for panel models}
\description{
  Specification test for panel models.
}
\usage{
phtest(x, ...)
\method{phtest}{panelmodel}(x,x2,data=NULL, ...)
\method{phtest}{plms}(x,...)
}

\arguments{
  \item{x}{an object of class \code{"plms"} or \code{"panelmodel"},}
  \item{x2}{an object of class \code{"panelmodel"},}
  \item{data}{}
  \item{...}{further arguments.}
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
  Hausman, J.A. (1978), Specification tests in econometrics,
  \emph{Econometrica}, 46, pp.1251--1271.
}

\author{Yves Croissant}
  
\examples{
library(Ecdat)
data(Gasoline)
Gasoline <- pdata.frame(Gasoline,"country","year")
form <- lgaspcar~lincomep+lrpmg+lcarpcap
gas <-  plm(form,data=Gasoline)
phtest(gas)
phtest(gas$between, gas$random)
}
\keyword{htest}
