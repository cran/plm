\name{pgmm}
\alias{pgmm}
\alias{summary.pgmm}
\alias{print.summary.pgmm}
\title{General method of moments estimator for panel data}
\description{
  General method of moments estimator for static or dynamic models with panel data.
}
\usage{
pgmm(formula,data,effect="individual",model="twosteps",
instruments=NULL,inst.transformation="l",lags.endog=1,
first.period=-99,last.period=-1,...)
\method{summary}{pgmm}(object, ...)
\method{print}{summary.pgmm}(x,digits=5,length.line=70, ...)
}

\arguments{
  \item{formula}{a symbolic description for the model to be estimated or
    a character string for a pure auto--regressive model,}
  \item{object,x}{an object of class \code{pgmm},}
  \item{data}{the data, must be an object of class \code{pdata.frame}
    and is mandatory,}
  \item{effect}{the effects introduced in the model, one of
    \code{"individual"} or \code{"twoways"},}
  \item{model}{one of \code{"onestep"} or \code{"twosteps"},}
  \item{instruments}{a one side formula containing the instruments,}
  \item{inst.transformation}{a vector of character strings containing
    \code{"l"} for instruments in level and \code{"d"} for instruments
    in first difference,}
  \item{lags.endog}{an integer between 0 and 2 indicating the number of
    lags of the dependent variable,}
  \item{first.period}{a vector of integer containing the first period
    for which the instrument has to be introduced,}
  \item{last.period}{a vector integer containing the last period for
    which the instrument has to be introduced,}
  \item{digits}{digits,}
  \item{length.line}{the maximum length of the lines in the print output,}
  \item{...}{further arguments.}
}

\value{
  an object of class \code{c("pgmm","panelmodel")}, which has the following elements :
  %  \describe{
  \item{coefficients}{the vector (or the list for fixed effects) of coefficients,}
  \item{residuals}{the vector of residuals,}
  \item{fitted.values}{the vector of fitted.values,}
  \item{vcov}{the covariance matrix of the coefficients,}
  \item{df.residual}{degrees of freedom of the residuals,}
  \item{model}{a data.frame containing the variables used for the
    estimation,}
  \item{omega}{a matrix containing the instruments (each column is an individual),}
  \item{WDWm1}{the inverse of the norm matrix used in the gmm estimation,}
  \item{J}{a 3D table containing the first period, last period and
    number of periods for each instrument and each year,}
  \item{K}{a list containing \code{K} the number of explanatory
    variables, \code{Ky} the number of lags of the dependent variable
    and \code{Kt} the number of time dummies,}
  \item{call}{the call.}
  %  }
  It has \code{print}, \code{summary} and \code{print.summary} methods.
}
\details{
  \code{pgmm} estimate a model for panel data with the general method of
  moments estimator. A dynamic model is specified by fixing the
  \code{lags.endog} argument to 1 or 2. By default, the instruments used
  are the independent variables (and the dependent variable if
  \code{lags.endog>0}). The complete list of instruments can also be
  specified with the \code{instruments} argument. For each instruments,
  the first period used is specified by \code{first.period} , the last
  period by \code{last.period} and the way they are introduced (in level
  or in first difference by \code{inst.transformation}). These three
  arguments may be of length one, in this case, all the instruments are
  treated the same way. In case of a dynamic model, they may be of
  length two, then the first   element is used for the dependent
  variable and the second element is used for the independent
  variables. They finally may be of length equal to the number of
  instruments to specify a specific rule for each instrument.}

\references{
  Arellano, Manuel & Bond, Stephen (1991), Some Tests of Specification
  for Panel Data: Monte Carlo Evidence and an Application to Employment
  Equations, \emph{The Review of Economic Studies}, vol. 58(2), april
  1991, pp.227--297.
}
\seealso{
  \code{\link{plm}} for the estimation of models with instrumental variables.
}

\author{Yves Croissant}

\examples{
library(Ecdat)
data(Snmesp)
pdata.frame(Snmesp,"firm","year")
z <-  pgmm(n~lag(w)+lag(w,2),Snmesp,effect="twoways",model="twosteps",lags.endog=2,
          last.period=c(-2),
          inst.transformation=c("l"),instruments=~n+w)
summary(z)

}
\keyword{regression}
