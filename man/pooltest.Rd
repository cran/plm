\name{pooltest}
\alias{pooltest}
\title{Test of poolability}
\description{
  A Chow  test for the poolability of the data.
}
\usage{
pooltest(plms,nopool=NULL,effect=FALSE)
}

\arguments{
  \item{plms}{an object of class \code{"plms"},}
  \item{nopool}{an object of class \code{"nopool"} (if the \code{plms}
    object does not contain any \code{"nopool"} model,}
  \item{effect}{If \code{TRUE} the \code{nopool} model is compared to
    the \code{within} model, if \code{FALSE} to the \code{pooling} model.}
}

\value{
  An object of class \code{"htest"}.
}

\details{
  \code{pooltest} is a \emph{F} test of stability (or Chow
  test) for the coefficients of a panel model.
}

\examples{
library(Ecdat)
data(Gasoline)
Gasoline <- pdata.frame(Gasoline,"country","year")
form <- lgaspcar~lincomep+lrpmg+lcarpcap
gas <-  plm(form,data=Gasoline,np=TRUE)
pooltest(gas)
pooltest(gas,effect=TRUE)
}
\keyword{htest}
