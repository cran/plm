\name{papply}
\alias{papply}
\alias{papply.pserie}
\alias{papply.matrix}
\alias{papply.default}

\title{Conditional operation for panel data}
\description{Apply an operation to a vector or a matrix
  for panel data.}

\usage{
papply(x, ...)
\method{papply}{pserie}(x,func,effect="individual", ...)
\method{papply}{matrix}(x,func,cond, ...)
\method{papply}{default}(x,func,cond, ...)
}

\arguments{
  \item{x}{a numeric vector or matrix,}
  \item{func}{the function to apply to \code{x},}
  \item{cond}{a vector containing the conditional variable,}
  \item{effect}{the effect for the papply transformation (one of
    \code{individual} or \code{time}),}
  \item{...}{further arguments.}
}

\value{
  a vector (or a matrix) of the same dimension as the first argument.
}


\examples{
library(Ecdat)
data(Grunfeld)
pdata.frame(Grunfeld,"firm","year")
papply(Grunfeld$inv,mean)
X <- model.matrix(~value+capital-1,Grunfeld)
papply(X,mean,Grunfeld$firm)
}
\keyword{arith}
