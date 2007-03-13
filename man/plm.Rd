\name{plm}
\alias{plm}
\alias{print.panelmodel}
\alias{summary.plm}
\alias{print.plms}
\alias{summary.plms}
\alias{print.summary.plm}
\alias{print.summary.plms}
\title{Panel Data Estimators}
\description{
  Linear models for panel data estimated using the \code{lm} function to
  transformed data.
}
\usage{
plm(formula,data,effect="individual",model=NULL,instruments=NULL,endog=NULL,
random.method="swar",inst.method="bvk", ...)
\method{summary}{plm}(object, ...)
\method{print}{plms}(x,digits=5, ...)
\method{summary}{plms}(object, ...)
\method{print}{summary.plm}(x,digits=5,length.line=70, ...)
\method{print}{summary.plms}(x,digits=5,length.line=70, ...)


}

\arguments{
  \item{formula}{a symbolic description for the model to be estimated,}
  \item{object,x}{an object of class \code{plm} or \code{plms},}
  \item{data}{the data, must be an object of class \code{pdata.frame}
    and is mandatory,}
  \item{effect}{the effects introduced in the model, one of \code{"individual"}, \code{"time"} or \code{"twoways"},}
  \item{model}{one of \code{"pooling"}, \code{"within"},
    \code{"between"}, \code{"random",} \code{"nopool"} and \code{"ht"} : \code{plm}
    returns the model specified or, if \code{NULL}, a list containing 
    four models (\code{"pooling"}, \code{"within"},
    \code{"between"} and  \code{"random"}),}
  \item{instruments}{a one side formula containing instrumental
    variables,}
  \item{endog}{a one side formula containing endogenous variables,}
  \item{random.method}{method of estimation for the variance components in the
    random effect model, one of \code{"swar"} (the default value), \code{"amemiya"}, \code{"walhus"} and \code{"nerlove"},}
  \item{inst.method}{the instrumental variable transformation : one of
    \code{"bvk"} and \code{"baltagi"},}
  \item{digits}{digits,}
  \item{length.line}{the maximum length of the lines in the print output,}
  \item{...}{further arguments.}
}

\value{
  Whether :
  an object of class \code{"plms"}, which is a list of the
  following models : \code{pooling}, \code{between} (\code{between.id} and
  \code{between.time} if \code{method="twoways"}), \code{within} and
  \code{random} which are all of class \code{"plm"},
  
  an object of class \code{c("plm","panelmodel")} if the argument
  \code{model} is filled.

  A \code{"plm"} object  has the following elements :
  
  %  \describe{
  \item{coefficients}{the vector of coefficients,}
  \item{residuals}{the vector of residuals,}
  \item{fitted.values}{the vector of fitted.values,}
  \item{vcov}{the covariance matrix of the coefficients,}
  \item{df.residual}{degrees of freedom of the residuals,}
  \item{model}{a data.frame containing the variables used for the
    estimation,}
  \item{call}{the call,}
  \item{FE}{the fixed effects (only for within models),}
  \item{alpha}{the overall intercept (only for within models),}
  \item{theta}{the parameter of transformation (only for random effect
    models),}
  \item{sigma2}{the variance of the different elements of the error
    (only for random effect models),}
    %  }
    
    It has \code{print}, \code{summary} and \code{print.summary} methods.

    A specific \code{summary} method is provided for objects of class \code{"plms"}, which returns an object of
    class \code{summary.plms} and prints a table of the coefficients
    of the within and random models and their standard errors.
    
}
\details{
  \code{plm} is a general function for the estimation of linear
  panel models. It offers limited support for unbalanced panels and
  estimation of two-ways effects models.

  For random effect models, 4 estimators of the transformation
  parameter are available : \code{swar} (Swamy and Arora),
  \code{amemiya}, \code{walhus} (Walhus and Hussain) and \code{nerlove}.

  Instrumental variables estimation is obtained using the
  \code{instruments} and/or \code{endog} arguments. If for example, the
  model is \code{y~x1+x2+x3},  \code{x1}, \code{x2} are endogenous and
  \code{z1}, \code{z2} are external
  instruments, the model can be estimated with :
  \code{instruments=~x3+z1+z2}, or
  \code{instruments=~z1+z2,endog=~x1+x2}. The four models are estimated 
  using  Balestra and Varadharajan--Krishnakumar's  method if
  \code{inst.method=bvk} or Baltagi's method if \code{inst.method="baltagi"}.

  The Hausman and Taylor estimator is  computed if \code{model="ht"}.

}

\author{Yves Croissant}

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
  \code{\link{pdata.frame}} for the creation of a \code{pdata.frame}.
}
\examples{
library(Ecdat)
data(Produc)
Produc <-pdata.frame(Produc,"state","year")
zz <- plm(log(gsp)~log(pcap)+log(pc)+log(emp)+unemp,data=Produc)
summary(zz$random)
}
\keyword{regression}
