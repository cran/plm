\name{pvar}
\alias{pvar}
\alias{pvar.default}
\alias{pvar.data.frame}
\alias{pvar.pdata.frame}
\alias{pvar.matrix}
\alias{print.pvar}
\title{Check whether variables of a panel have individual and time variations}
\description{This function checks for each variable of a panel if it has an
  individual and a time variation.}
\usage{
pvar(x, ...)
\method{pvar}{pdata.frame}(x, ...)
\method{pvar}{data.frame}(x,id,time, ...)
\method{pvar}{matrix}(x,id,time, ...)
%\method{pvar}{default}(x,id,time, ...)
}
\arguments{
  \item{x}{a  \code{data.frame}, a \code{pdata.frame} or a \code{matrix}
    object,}
  \item{id}{a vector containing the individual index,}
  \item{time}{a vector containing the time index,}
  \item{...}{further arguments.}
}

\value{
  an object of class \code{pvar}  containing the following elements :
  \item{id.var}{ a logical vector with \code{TRUE} values if the
    variable has individual variation, \code{FALSE} otherwise,}
  \item{time.var}{ a logical vector with \code{TRUE} values if the
    variable has time variation, \code{FALSE} otherwise,}
}
\details{
\code{pvar} is called by \code{pdata.frame} and by \code{plm} in
there are some missing values.
}
\seealso{
  \code{\link{pdata.frame}} and \code{\link{plm}}.
}

\examples{
library(Ecdat)
data(Produc)
pvar(Produc,"state","year")
}
\keyword{attribute}
