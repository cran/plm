\name{pggls}
\alias{pggls}
\alias{summary.pggls}
\alias{print.summary.pggls}
\title{General FGLS Estimators}
\description{
  General FGLS estimators for panel data (balanced or unbalanced)
}
\usage{
pggls(formula,data,effect="individual",model, ...)
\method{summary}{pggls}(object, ...)
\method{print}{summary.pggls}(x,digits=5,length.line=70, ...)
}

\arguments{
  \item{formula}{a symbolic description for the model to be estimated,}
  \item{object,x}{an object of class \code{pggls},}
  \item{data}{the data, must be an object of class \code{pdata.frame}
    and is mandatory,}
  \item{effect}{the effects introduced in the model, one of \code{"individual"} or \code{"time"},}
  \item{model}{one of \code{"within"} or \code{"random"},}
  \item{digits}{digits,}
  \item{length.line}{the maximum length of the lines in the print output,}
  \item{...}{further arguments.}
}

\value{
  an object of class \code{c("pggls","panelmodel")} containing :
  
  %  \describe{
    \item{coefficients}{the vector of coefficients,}
    \item{residuals}{the vector of residuals,}
    \item{fitted.values}{the vector of fitted.values,}
    \item{vcov}{the covariance matrix of the coefficients,}
    \item{df.residual}{degrees of freedom of the residuals,}
    \item{model}{a data.frame containing the variables used for the
    estimation,}
  \item{call}{the call,}
    \item{sigma}{the estimated intragroup (or cross-sectional, if \code{effect="time"}) covariance of errors,}
    %  }
  }
  
\details{
  \code{pggls} is a function for the estimation of linear  panel models
  by general feasible generalized least squares, either with or without
  fixed effects. General FGLS is based on a two-step estimation process:
  first a model is estimated by OLS (\code{random}) or fixed effects
  (\code{within}), then its residuals are used to estimate an error
  covariance matrix for use in a feasible-GLS analysis. This framework
  allows the error covariance structure inside every group (if
  \code{effect="individual"}, else symmetric) of observations to be
  fully unrestricted and is therefore robust against any type of
  intragroup heteroskedasticity and serial correlation. This structure,
  by converse, is assumed identical across groups and thus general FGLS
  estimation is inefficient under groupwise heteroskedasticity. Care
  shall also be taken that this method requires estimation of T(T+1)/2
  variance parameters, thus efficiency requires N  > >T (if \code{effect="individual"}, else the opposite).}

\author{Giovanni Millo}

\references{
  Kiefer, N. M. (1980), Estimation of Fixed Effects Models for Time Series of Cross-Sections with Arbitrary Intertemporal Covariance, \emph{Journal of Econometrics}, \bold{14}, 195--202. 

  Wooldridge J. M. (2003), \emph{Econometric Analysis of Cross Section and Panel Data}, MIT Press
}
\seealso{
  \code{\link{pdata.frame}} for the creation of a \code{pdata.frame}.
}
\examples{
library(Ecdat)
data(Produc)
Produc <-pdata.frame(Produc,"state","year")
zz <- pggls(log(gsp)~log(pcap)+log(pc)+log(emp)+unemp,data=Produc,model="random")
summary(zz)
}
\keyword{regression}
