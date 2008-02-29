\name{pooltest}
\alias{pooltest}
\alias{pooltest.plm}
\alias{pooltest.formula}
\title{Test of Poolability}
\description{
  A Chow  test for the poolability of the data.
}
\usage{
pooltest(x, ...)
\method{pooltest}{plm}(x, z, ...)
\method{pooltest}{formula}(x, data, effect = "individual", model =
"within", index = NULL, ...)
}

\arguments{
  \item{x}{an object of class  \code{"plm"},}
  \item{z}{an object of class \code{"pvcm"} obtained with
    \code{model="within"},}
  \item{data}{a \code{data.frame},}
  \item{effect}{one of \code{"individual"} or \code{"time"},}

  \item{model}{one of \code{"pooling"} or \code{"within"}; intercepts
  are assumed to be identical in the first case and different in the
  second case.}
  \item{index}{the index of the \code{data.frame},}

  \item{\dots}{further arguments.}
}

\value{
  An object of class \code{"htest"}.
}

\details{
  \code{pooltest} is a \emph{F} test of stability (or Chow
  test) for the coefficients of a panel model.
}

\author{Yves Croissant}


\examples{
data("Gasoline", package="Ecdat")
form <- lgaspcar~lincomep+lrpmg+lcarpcap
gasw <-  plm(form, data=Gasoline, model="within")
gasp <-  plm(form, data=Gasoline, model="pooling")
gasnp <-  pvcm(form, data=Gasoline, model="within")
pooltest(gasw, gasnp)
pooltest(gasp, gasnp)

pooltest(form, data=Gasoline, effect="individual", model="within")
pooltest(form, data=Gasoline, effect="individual", model="pooling")



}
\keyword{htest}
