\name{plm}
\alias{plm}
\alias{plm.default}
\alias{plm.formula}
\alias{print.plm}
\alias{summary.plm}
\alias{print.plms}
\alias{summary.plms}
\alias{print.summary.plm}
\alias{print.summary.plms}
\title{Panel Data Estimators}
\description{
  Estimators for panel data (balanced or unbalanced)
}
\usage{
plm(y, ...)
\method{plm}{formula}(y,instruments=NULL,endog=NULL,data,effect="individual",
theta="swar",trinst="baltagi",model=NULL,np=FALSE,...)
\method{plm}{default}(y,X,W,id,time,pvar,pdim,pmodel, ...)
\method{print}{plm}(x,digits=3, ...)
\method{summary}{plm}(object, ...)
\method{print}{plms}(x,digits=3, ...)
\method{summary}{plms}(object, ...)
\method{print}{summary.plm}(x,digits=3, ...)
\method{print}{summary.plms}(x,digits=3, ...)


}

\arguments{
  \item{y}{a symbolic description for the model to be estimated
    for the \code{formula} method, a numeric vector for the default method,}
  \item{object,x}{an object of class \code{plm} or \code{plms},}
  \item{instruments}{a one side formula containing instrumental
    variables,}
  \item{endog}{a one side formula containing endogenous variable,}
  \item{data}{the data, must be an object of class \code{pdata.frame}
    and is compulsary,}
  \item{effect}{one of \code{"individual"}, \code{"time"} or \code{"twoways"} for a two way estimation,}
  \item{theta}{method of estimation for the variance components in the
    random effect model, one of \code{"swar"}, \code{"amemiya"}, \code{"walhus"} and \code{"nerlove"},}
  \item{trinst}{the instrumental variable transformation : one of
    \code{"baltagi"}, \code{"bvk"}, \code{"ht"},}
  \item{model}{one of \code{"pooling"}, \code{"within"},
    \code{"between"} and \code{"random"} or \code{NULL} : \code{plm}
    returns the model spectified or if \code{NULL }a list containing the
    fout models,}
  \item{W}{a matrix of instrumental variables,}
  \item{X}{a matrix of explanatory variables,}
  \item{id}{the individual index,}
  \item{time}{the time index,}
  \item{pvar}{a list resulting from a call to \code{pvarcheck},}
  \item{pdim}{a list resulting from a call to \code{pdimcheck},}
  \item{pmodel}{a list containing the characteristics of the model to be
  estimated : \code{model}, \code{formula}, \code{effect}, \code{theta},
\code{trinst},}
  \item{np}{a logical value which indicates whether the \code{nopool}
    model has to be estimated or not,}
  \item{digits}{digits,}
  \item{...}{further arguments.}

    
}

\value{
  Wheter :
  
  an object of class \code{"plms"}, which is a list of the
  following models : \code{pooling}, \code{between} (\code{between.id} and
  \code{between.time} if \code{method="twoways"}), \code{within} and
  \code{random} which are all of class \code{"plm"},

  an object of class \code{"plm"} if the argument \code{model} is filled
  or if \code{trinst="ht"}.

  A \code{"plm"} object is a  list of the following elements :
  \code{coefficients}, \code{df.residual}, \code{ssr},
  \code{cov.unscaled} and \code{formula}. It has \code{print}, \code{summary} and
  \code{print.summary} methods which are not unlike \code{lm}'s methods.

  A specific \code{summary} method is provided for objects of class \code{"plms"}, which returns an objects of
  class \code{summary.plms} and prints a table of the coefficients
  of the different models and their standard errors.
}

\details{
  \code{plm} is a general function for the estimation of linear
  panel models. It offers limited support for unbalanced panels and
  estimation of two-ways effects models.

  For random effect models, 4 estimators of the transformation
  parameter are available : "swar","amemiya","walhus" and "nerlove".

  Instrumental variable estimation is obtained using the
  \code{instruments} and/or \code{endog} arguments. If for example, the
  model is y~x1+x2+x3,  x1,x2 are endogenous and z1,z2 are external
  instruments, the model can be estimated with :
  \code{instruments=~x3+z1+z2}, or
  \code{instruments=~z1+z2,endog=~x1+x2}. The four models are estimated by
  instrumental variables if \code{trinstr} equal \code{"bvk"} (Balestra, P. and
  J. Varadharajan--Krishnakumar (1987)) or \code{"baltagi"} (Baltagi
  (1981)). If \code{trinstr="ht"}, the Hausman and Taylor estimator is computed
  and only a random effect model is returned.
}
\references{
  Amemiyia, T. (1971), The estimation of the variances in a
  variance--components model, \emph{International Economic Review}, 12,
  pp.1--13.

  Balestra, P. and J. Varadharajan--Krishnakumar (1987), Full
  information estimations of a system of simultaneous equations with
  error components structure, \emph{Econometric Theory}, 3, pp.223--246.
  
  Baltagi, B.H. (1981), Simultaneous equations with error components,
  \emph{Journal of econometrics}, 17, pp.21--49.
  
  Baltagi, B.H. (2001) \emph{Econometric Analysis of Panel Data}. John
  Wiley and sons. ltd.

  Hausman, J.A. and W.E. Taylor (1981), Panel data and unobservable
  individual effects, \emph{Econometrica}, 49, pp.1377--1398.
  
  Nerlove, M. (1971), Further evidence on the estimation of dynamic
  economic relations from a time--series of cross--sections,
  \emph{Econometrica}, 39, pp.359--382.

  Swamy, P.A.V.B. and S.S. Arora (1972), The exact finite sample
  properties of the estimators of coefficients in the error components
  regression models, \emph{Econometrica}, 40, pp.261--275.

  Wallace, T.D. and A. Hussain (1969), The use of error components
  models in combining cross section with time series data,
  \emph{Econometrica}, 37(1), pp.55--72.
}
\seealso{
  \code{\link{pdata.frame}} for the creation of a \code{pdata.frame}
}
\examples{
library(Ecdat)
data(Produc)
Produc <-pdata.frame(Produc,state,year)
zz <- plm(log(gsp)~log(pcap)+log(pc)+log(emp)+unemp,data=Produc)
summary(zz$random)
}
\keyword{regression}
