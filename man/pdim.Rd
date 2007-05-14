\name{pdim}
\alias{pdim}
\alias{pdim.pdata.frame}
\alias{pdim.data.frame}
\alias{pdim.default}
\alias{print.pdim}

\title{Check for the dimensions of the panel}
\description{This function checks the number of individuals and time
  observations in the panel and whether it is balanced or not.}
\usage{
pdim(x,...)
\method{pdim}{pdata.frame}(x,...)
\method{pdim}{data.frame}(x,id,time,...)
%\method{pdim}{default}(x,y,...)

}
\arguments{
  \item{x}{a \code{data.frame} or a \code{pdata.frame},}
%  or a vector for the default method}
%  \item{y}{a vector,}
  \item{id}{the individual index,}
  \item{time}{the time index,}
  \item{...}{further arguments.}
}

\value{
  an object of class \code{pdim} containing the following elements :
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
\code{pdim} is called by \code{pdata.frame} and by \code{plm} if
there are some missing values.
}
\seealso{
  \code{\link{pdata.frame}} and \code{\link{plm}}
}

\author{Yves Croissant}

\examples{
library(Ecdat)
data(Produc)
pdim(Produc,"state","year")
}
\keyword{attribute}
