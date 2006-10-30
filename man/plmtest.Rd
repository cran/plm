\name{plmtest}
\alias{plmtest}
%\alias{plmtest.default}
\alias{plmtest.plm}
\alias{plmtest.plms}
\title{Lagrance mulitplier tests for panel models}
\description{
  Test of individual and/or time effect for panel models.
}
\usage{
plmtest(x, ...)
\method{plmtest}{plm}(x,effect="id",type="bp", ...)
\method{plmtest}{plms}(x,effect="id",type="bp", ...)
%\method{plmtest}{default}(x,n=NULL,T=NULL,balanced=NULL,id=NULL,time=NULL,effect="id",type="bp",data=NULL, ...)
}

\arguments{
  \item{x}{an object of class \code{plm} or \code{plms},}
  \item{effect}{a character string indicating which effects are tested :
    individual effects (\code{"id"}), time effects (\code{"time"}) or both (\code{"twoways"}),}
  \item{type}{a character string indicating the  test to be
    computed : \code{"bp"} for Breush-Pagan (1980), \code{"ghm"} for Gourieroux, Holly and
    Monfort (1982), \code{"honda"} for Honda (1985) or \code{"kw"} for King and Wu (1997),}
%   \item{n}{the number of individuals,}
%   \item{T}{the number of time periods,}
%   \item{balanced}{A logical value indicating wheter the panel is
%     balanced or not,}
%   \item{id}{a vector containing the individual index,}
%   \item{time}{a vector containong the time index,}
%   \item{data}{the name of the data,}
  \item{...}{further arguments.}
}

\value{
  An object of class \code{"htest"}.
}

\details{
  These Lagrange multiplier tests use only the residuals of the pooling
  model. The main argument of this function may be whether a pooling model of class \code{plm} or an object of class \code{plms}.
}
\references{
  Breusch, T.S. and A.R. Pagan (1980), The Lagrange multiplier test and
  its applications to model specification in econometrics, \emph{Review
    of Economic Studies}, 47, pp.239--253.

  Gourieroux, C., A. Holly and A. Monfort (1982), Likelihood ratio test,
  Wald test, and Kuhn--Tucker test in linear models with inequality
  constraints on the regression parameters, \emph{Econometrica}, 50,
  pp.63--80.

  Honda, Y. (1985), Testing the error components model with non--normal
  disturbances, \emph{Review of Economic Studies}, 52, pp.681--690.

  King, M.L. and P.X. Wu (1997), Locally optimal one--sided tests for
  multiparameter hypotheses, \emph{Econometric Reviews}, 33,
  pp.523--529.

}
\seealso{
  \code{\link{pFtest}} for individual and/or time effects tests based on
  the within model.
}
\examples{
library(Ecdat)
data(Grunfeld)
Grunfeld <-pdata.frame(Grunfeld,"firm","year")
g <- plm(inv ~ value + capital, data=Grunfeld)
plmtest(g)
plmtest(g,effect="time")
plmtest(g,type="honda")
plmtest(g,type="ghm",effect="twoways")
plmtest(g,type="kw",effect="twoways")
}
\keyword{htest}
