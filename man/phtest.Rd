\name{phtest}
\alias{phtest}
\alias{phtest.formula}
\alias{phtest.panelmodel}
\title{Hausman Test for Panel Models}
\description{
  Specification test for panel models.
}
\usage{
phtest(x, ...)
\method{phtest}{panelmodel}(x, x2, ...)
\method{phtest}{formula}(x, data, model = c("within","random"),
    effect = "individual", index = NULL, ...)
}

\arguments{
  \item{x}{an object of class  \code{"panelmodel"} or \code{"formula"},}
  \item{x2}{an object of class \code{"panelmodel"},}
  \item{model}{a vector containing the names of two models,}
  \item{effect}{the effects introduced in the model, one of \code{"individual"}, \code{"time"} or \code{"twoways"},}
  \item{data}{a \code{data.frame},}
  \item{index}{the index of the \code{data.frame},}
  \item{\dots}{further arguments.}
}

\value{
  An object of class \code{"htest"}.
}

\details{
  The Hausman test is based on the difference of the vectors of
  coefficients of two different models.
}
\references{
  Hausman, J.A. (1978), Specification tests in econometrics,
  \emph{Econometrica}, \bold{46}, pp.1251--1271.
}

\author{Yves Croissant}
  
\examples{
data("Gasoline",package="Ecdat")
form <- lgaspcar~lincomep+lrpmg+lcarpcap
wi <-  plm(form, data=Gasoline, model="within")
re <-  plm(form, data=Gasoline, model="random")
phtest(wi, re)
phtest(form, data=Gasoline)
}
\keyword{htest}
