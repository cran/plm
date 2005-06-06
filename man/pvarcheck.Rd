\name{pvarcheck}
\alias{pvarcheck}
\title{Check wheter variables of a panel exhibit individual and time variations}
\description{This function checks for each variable of a panel if it has an
  individual and a time variation.}
\usage{
pvarcheck(id,time,data)
}
\arguments{
  \item{id}{a vector containing the individual index,}
  \item{time}{a vector containing the time index,}
  \item{data}{a  \code{data.frame}.}
}

\value{
  a list containing the following elements :
  \item{id.var}{ a logical vector with values \code{TRUE} if the
    variable has individual variation, \code{FALSE} otherwise,}
  \item{time.var}{ a logical vector with values \code{TRUE} if the
    variable has time variation, \code{FALSE} otherwise,}
}
\details{
\code{pvarcheck} is called by \code{pdata.frame} and by \code{plm} in
there are some missing values.
}
\seealso{
  \code{\link{pdata.frame}} and \code{\link{plm}}.
}

\examples{
library(Ecdat)
data(Produc)
pvarcheck(Produc$state,Produc$year,Produc)
}
\keyword{attribute}
