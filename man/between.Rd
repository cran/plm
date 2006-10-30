\name{between}
\alias{Between}
\alias{Between.pserie}
\alias{Between.matrix}
\alias{Between.default}
\alias{between}
\alias{between.pserie}
\alias{between.matrix}
\alias{between.default}

\title{The between transformation}
\description{Computes the between transformation}

\usage{
Between(x, ...)
\method{Between}{pserie}(x,effect="individual", ...)
\method{Between}{matrix}(x,cond, ...)
\method{Between}{default}(x,cond, ...)
between(x, ...)
\method{between}{pserie}(x,effect="individual", ...)
\method{between}{matrix}(x,cond, ...)
\method{between}{default}(x,cond, ...)

}

\arguments{
  \item{x}{a numeric vector or matrix,}
  \item{cond}{a vector containing the conditional variable,}
  \item{effect}{the effect for the Between transformation (one of
    \code{individual} and \code{time}),}
  \item{...}{further arguments.}
}

\value{
  a vector (or a matrix) of the same length as \code{x} containing the
  between transformation (\code{Between}) or of a length equal of the
  number of modalities of the conditional variable (\code{between}).
  \code{NA} values are inserted at the correct positions}

\seealso{
  \code{\link{within}} and \code{\link{papply}}.
  }

\examples{
library(Ecdat)
data(Grunfeld)
pdata.frame(Grunfeld,"firm","year")
Bi <- Between(Grunfeld$inv)
X <- model.matrix(~value+capital-1,Grunfeld)
BX <- Between(X,Grunfeld$firm)
bi <- between(Grunfeld$inv)
bX <- between(X,Grunfeld$firm)
# computes the between estimator with the lm command
blm <- lm(between(inv)~between(value)+between(capital),data=Grunfeld)
}
\keyword{arith}
