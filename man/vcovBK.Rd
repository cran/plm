% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tool_vcovG.R
\name{vcovBK}
\alias{vcovBK}
\alias{vcovBK.plm}
\title{Beck and Katz Robust Covariance Matrix Estimators}
\usage{
vcovBK(x, ...)

\method{vcovBK}{plm}(
  x,
  type = c("HC0", "sss", "HC1", "HC2", "HC3", "HC4"),
  cluster = c("group", "time"),
  diagonal = FALSE,
  ...
)
}
\arguments{
\item{x}{an object of class \code{"plm"},}

\item{\dots}{further arguments.}

\item{type}{the weighting scheme used, one of \code{"HC0"}, \code{"sss"}, \code{"HC1"},
\code{"HC2"}, \code{"HC3"}, \code{"HC4"}, see Details,}

\item{cluster}{one of \code{"group"}, \code{"time"},}

\item{diagonal}{a logical value specifying whether to force
non-diagonal elements to zero,}
}
\value{
An object of class \code{"matrix"} containing the estimate of
the covariance matrix of coefficients.
}
\description{
Unconditional Robust covariance matrix estimators \emph{a la Beck
and Katz} for panel models (a.k.a. Panel Corrected Standard Errors
(PCSE)).
}
\details{
\code{vcovBK} is a function for estimating a robust covariance matrix of
parameters for a panel model according to the
\insertCite{BECK:KATZ:95;textual}{plm} method, a.k.a. Panel
Corrected Standard Errors (PCSE), which uses an unconditional
estimate of the error covariance across time periods (groups)
inside the standard formula for coefficient
covariance. Observations may be clustered either by \code{"group"} to
account for timewise heteroskedasticity and serial correlation or
by \code{"time"} to account for cross-sectional heteroskedasticity and
correlation. It must be borne in mind that the Beck and Katz
formula is based on N- (T-) asymptotics and will not be appropriate
elsewhere.

The \code{diagonal} logical argument can be used, if set to
\code{TRUE}, to force to zero all non-diagonal elements in the
estimated error covariances; this is appropriate if both serial and
cross--sectional correlation are assumed out, and yields a
timewise- (groupwise-) heteroskedasticity--consistent estimator.

Weighting schemes specified by \code{type} are analogous to those in
\code{\link[sandwich:vcovHC]{sandwich::vcovHC()}} in package \CRANpkg{sandwich} and are
justified theoretically (although in the context of the standard
linear model) by \insertCite{MACK:WHIT:85;textual}{plm} and
\insertCite{CRIB:04;textual}{plm} \insertCite{@see @ZEIL:04}{plm}.

The main use of \code{vcovBK} (and the other variance-covariance estimators
provided in the package \code{vcovHC}, \code{vcovNW}, \code{vcovDC}, \code{vcovSCC}) is to pass
it to plm's own functions like \code{summary}, \code{pwaldtest}, and \code{phtest} or
together with testing functions from the \code{lmtest} and \code{car} packages. All of
these typically allow passing the \code{vcov} or \code{vcov.} parameter either as a
matrix or as a function, e.g., for Wald--type testing: argument \code{vcov.} to
\code{coeftest()}, argument \code{vcov} to \code{waldtest()} and other methods in the
\CRANpkg{lmtest} package; and argument \code{vcov.} to
\code{linearHypothesis()} in the \CRANpkg{car} package (see the
examples), see \insertCite{@see also @ZEIL:04}{plm}, 4.1-2, and examples below.
}
\examples{

data("Produc", package="plm")
zz <- plm(log(gsp)~log(pcap)+log(pc)+log(emp)+unemp, data=Produc, model="random")
summary(zz, vcov = vcovBK)
summary(zz, vcov = function(x) vcovBK(x, type="HC1"))

## standard coefficient significance test
library(lmtest)
coeftest(zz)
## robust significance test, cluster by group
## (robust vs. serial correlation), default arguments
coeftest(zz, vcov.=vcovBK)
## idem with parameters, pass vcov as a function argument
coeftest(zz, vcov.=function(x) vcovBK(x, type="HC1"))
## idem, cluster by time period
## (robust vs. cross-sectional correlation)
coeftest(zz, vcov.=function(x) vcovBK(x, type="HC1", cluster="time"))
## idem with parameters, pass vcov as a matrix argument
coeftest(zz, vcov.=vcovBK(zz, type="HC1"))
## joint restriction test
waldtest(zz, update(zz, .~.-log(emp)-unemp), vcov=vcovBK)
\dontrun{
## test of hyp.: 2*log(pc)=log(emp)
library(car)
linearHypothesis(zz, "2*log(pc)=log(emp)", vcov.=vcovBK)
}
}
\references{
\insertRef{BECK:KATZ:95}{plm}

\insertRef{CRIB:04}{plm}

\insertRef{GREE:03}{plm}

\insertRef{MACK:WHIT:85}{plm}

\insertRef{ZEIL:04}{plm}
}
\seealso{
\code{\link[sandwich:vcovHC]{sandwich::vcovHC()}} from the \CRANpkg{sandwich}
package for weighting schemes (\code{type} argument).
}
\author{
Giovanni Millo
}
\keyword{regression}
