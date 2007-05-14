\name{pvcm}
\alias{pvcm}
\alias{summary.pvcm}
\alias{print.summary.pvcm}
\title{Variable coefficients models for panel data}
\description{
  Estimators for random and fixed effect models with variable coefficients.
}
\usage{
pvcm(formula,data,effect="individual",model,...)
\method{summary}{pvcm}(object, ...)
\method{print}{summary.pvcm}(x,digits=5,length.line=70, ...)

}

\arguments{
  \item{formula}{a symbolic description for the model to be estimated,}
  \item{object,x}{an object of class \code{pvcm},}
  \item{data}{the data, must be an object of class \code{pdata.frame}
    and is mandatory,}
  \item{effect}{the effects introduced in the model, one of
    \code{"individual"} or \code{"time"},}
  \item{model}{one of \code{"within"} or \code{"random"},}
  \item{digits}{digits,}
  \item{length.line}{the maximum length of the lines in the print output,}
  \item{...}{further arguments.}
}

\value{
  an object of class \code{c("pvcm","panelmodel")}, which has the following elements :
  %  \describe{
    \item{coefficients}{the vector (or the list for fixed effects) of coefficients,}
    \item{residuals}{the vector of residuals,}
    \item{fitted.values}{the vector of fitted.values,}
    \item{vcov}{the covariance matrix of the coefficients,}
    \item{df.residual}{degrees of freedom of the residuals,}
    \item{model}{a data.frame containing the variables used for the
      estimation,}
    \item{call}{the call,}
    \item{Delta}{the estimation of the covariance matrix of the
      coefficients (random effect models only),}
    \item{std.error}{the standard errors for all the coefficients for
      each individual (within models only),}
    %  }
    \code{pvcm} objects have \code{print}, \code{summary} and \code{print.summary} methods.
}
\details{
  The \code{pvcm} function enables the estimation of variable
  coefficients models. Time or individual effects are introduced if
  \code{effect} is fixed to \code{"time"} or \code{"individual"}
  (the default value). 
  
  Coefficients are assumed to be fixed if \code{model="within"} and
  random if \code{model="random"}. In the first case, a different
  model is estimated for each individual (or time period). In the second
  case, the \textsc{Swamy} (1970) model is estimated. It is a
  generalized least squares model which use the results of the previous model.
  
}

\author{Yves Croissant}

\references{
  Swamy, P.A.V.B. (1970), Efficient Inference in a Random Coefficient
  Regression Model, \emph{Econometrica}, 38(2), pp.311--323.
}

\examples{
library(Ecdat)
data(Produc)
pdata.frame(Produc,"state","year","pprod")
zw <- pvcm(log(gsp)~log(pcap)+log(pc)+log(emp)+unemp,data=pprod,model="within")
zr <- pvcm(log(gsp)~log(pcap)+log(pc)+log(emp)+unemp,data=pprod,model="random")
}
\keyword{regression}
