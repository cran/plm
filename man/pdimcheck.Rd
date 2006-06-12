\name{pdimcheck}
\alias{pdimcheck}
\title{Check for the dimensions of the panel}
\description{This function check the number of individuals and time
  observations in the panel and wheter it is balanced or not.}
\usage{
pdimcheck(id,time,data=NULL)
}
\arguments{
  \item{id}{a vector containing the individual index,}
  \item{time}{a vector containing the time index,}
  \item{data}{an optional \code{data.frame}.}
}

\value{
  a list containing the following elements :
  \item{nT}{ a list containing \code{n}, the number of individuals,
    \code{T}, the number of time observations, \code{N} the total number
    of observations,}
  \item{TinT}{ a list containing two vectors : \code{Ti} gives the
    number of observations for each individuals and \code{nt} gives
    the number of individuals observed for each period,}
  \item{balanced}{ a logical value : \code{TRUE} for a balanced panel,
    \code{FALSE} for an unbalanced panel,}
  \item{panel.names}{a list of character vectors : \code{id.names}
    contains the names of each individual and \code{time.names}
    contains the names of each period.}
    }
\details{
\code{pdimcheck} is called by \code{pdata.frame} and by \code{plm} in
there are some missing values.
}
\seealso{
  \code{\link{pdata.frame}} and \code{\link{plm}}
}

\examples{
library(Ecdat)
data(Produc)
pdimcheck(state,year,Produc)
pdimcheck(Produc$state,Produc$year)
}
\keyword{attribute}
