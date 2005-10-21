\name{within}
\alias{within}
\alias{within.pserie}
\alias{within.matrix}
%\alias{within.default}

\title{The within transformation}
\description{Computes the within transformation}

\usage{
within(x, ...)
\method{within}{pserie}(x,effect="individual", ...)
\method{within}{matrix}(x,cond, ...)
%\method{within}{default}(x,cond, ...)
}

\arguments{
  \item{x}{a numeric vector or a matrix,}
  \item{cond}{a vector containing the conditional variable,}
  \item{effect}{the effect for the within transformation (one of
    \code{individual}, \code{time} and \code{twoways}),}
  \item{...}{further arguments.}
}

\value{
  a vector of the same length as \code{x} containing the
  within transformation. \code{NA} values are inserted at the correct positions}

\seealso{\code{\link{between}} and \code{\link{papply}}}

\examples{
library(Ecdat)
data(Grunfeld)
pdata.frame(Grunfeld,"firm","year")
wi <- within(Grunfeld$inv)
w2 <- within(Grunfeld$inv,effect="twoways")
X <- model.matrix(~value+capital-1,Grunfeld)
Xt <- within(X,Grunfeld$firm)
}
\keyword{arith}
