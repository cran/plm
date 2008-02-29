\name{pbsytest}
\alias{pbsytest}
\alias{pbsytest.formula}
\title{Bera, Sosa-Escudero and Yoon Locally--Robust Lagrange Mulitiplier
  Tests for Panel Models}

\description{

  Test for residual serial correlation (or individual random effects)
  locally robust vs. individual random effects (serial correlation) for
  panel models and joint test by Baltagi and Li.

}
\usage{
pbsytest(x, ...)
\method{pbsytest}{formula}(x, data, test = c("AR","RE","J"), index =
NULL, ...)
}

\arguments{
  \item{x}{an object of class \code{"formula"},}
  \item{data}{a \code{"data.frame"},}
  \item{test}{a character string indicating which test to perform :
    first-order serial correlation (\code{AR}), random effects
    (\code{RE}) or joint test for either of them (\code{J}).}
  \item{index}{the index of the \code{data.frame},}
  \item{\dots}{further arguments.}
}

\value{
  An object of class \code{"htest"}.
}

\details{

  These Lagrange multiplier tests are robust vs. local misspecification
  of the alternative hypothesis, i.e. they test the null of serially
  uncorrelated residuals against AR(1) residuals in a pooling model,
  allowing for local departures from the assumption of no random
  effects; or they test the null of no random effects allowing for local
  departures from the assumption of no serial correlation in
  residuals. They use only the residuals of the pooled OLS model and
  correct for local misspecification as outlined in Bera and \emph{al.} (2001).

The joint test is due to Baltagi and Li (1991) and is added for
convenience under this same interface.

}
\references{

  Bera, A.K., Sosa-Escudero, W. and  Yoon, M. (2001), Tests for the error
   component model in the presence of local misspecification,
   \emph{Journal of Econometrics}, \bold{101}, pp. 1-23.

   Baltagi, B. and Li, Q. (1991), A joint test for serial correlation and
   random individual effects, \emph{Statistics and Probability Letters},
   \bold{11}, pp. 277-280.
   
}
\seealso{
  
  \code{\link{plmtest}} for individual and/or time random effects tests
  based on a correctly specified model; \code{\link{pbltest}},
  \code{\link{pbgtest}} and \code{\link{pdwtest}} for serial correlation
  tests in random effects models.

}
\examples{
## Example in Bera et al.
data(Grunfeld, package="Ecdat")
## Bera et al. use a subset of the original Grunfeld data,
## so results are slightly different here
## default is AR testing
pbsytest(inv ~ value + capital, data=Grunfeld)
pbsytest(inv ~ value + capital, data=Grunfeld, test="RE")
pbsytest(inv ~ value + capital, data=Grunfeld, test="J")
}
\keyword{htest}
