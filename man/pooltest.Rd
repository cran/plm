\name{pooltest}
\alias{pooltest}
\alias{pooltest.plm}
\alias{pooltest.plms}
\title{Test of poolability}
\description{
  A Chow  test for the poolability of the data.
}
\usage{
pooltest(x,...)
\method{pooltest}{plms}(x,nopool,effect=FALSE,...)
\method{pooltest}{plm}(x,nopool,...)
}

\arguments{
  \item{x}{an object of class \code{"plms"} or \code{"plm"},}
  \item{nopool}{an object of class \code{"pvcm"} obtained with \code{effect="within"},}
  \item{effect}{if \code{x} is a \code{plms} and \code{effect=TRUE} the \code{nopool} model is compared to
    the \code{within} model, if \code{FALSE} to the \code{pooling}
    model.}
  \item{...}{further arguments.}
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
library(Ecdat)
data(Gasoline)
Gasoline <- pdata.frame(Gasoline,"country","year")
form <- lgaspcar~lincomep+lrpmg+lcarpcap
gas <-  plm(form,data=Gasoline)
gasnp <-  pvcm(form,data=Gasoline,model="within")
pooltest(gas,gasnp)
pooltest(gas,gasnp,effect=TRUE)
}
\keyword{htest}
