\name{nopool}
\alias{nopool}
\alias{nopool.default}
\alias{nopool.formula}
\alias{print.nopool}
\alias{summary.nopool}
\alias{print.summary.nopool}
\alias{plot.nopool}
\title{No pooling estimation}
\description{
  Least Squares estimation for each individual or time observation
}
\usage{
nopool(y,...)
\method{nopool}{formula}(y,data=data,effect="individual",...)
\method{nopool}{default}(y,X,cond,...)
\method{print}{nopool}(x,digits=3, ...)
\method{summary}{nopool}(object, ...)
\method{print}{summary.nopool}(x,digits=3, ...)
\method{plot}{nopool}(x,...)
}

\arguments{
  \item{y}{a symbolic description for the model to be estimated for the
    formula method, a numerical vector for the default method,}
  \item{cond}{a vector containing the individual or the time index,}
  \item{object,x}{an object of class \code{nopool},}
  \item{X}{ a matrix of explanatory variables,}
  \item{data}{the data,}
  \item{effect}{one of \code{"individual"} or \code{"time"},
    respectively for an estimation for each individual or time observation,}
  \item{digits}{digits,}
  \item{...}{further arguments.}   
}

\value{
  an object of class \code{"nopool"}, which is a list of  :

  \code{residuals} the vector of residuals,
  
  \code{df.residuals} the degrees of freedom of the regression,

  \code{ssr} the sum of the squared residuals,

  \code{coefficients} a \code{data.frame} containing the coefficients for each
  estimation,
  
  \code{std.error} a \code{data.frame} containing the standard errors of the
  coefficients for each estimation,
  

  The \code{print}  prints the \code{data.frame} containing the
  coefficients, and the \code{summary} method is the
  \code{data.frame}'s method applied to the coefficients.

  The \code{plot} method plots histograms for each coefficient.
}

\details{
  The \code{nopool} model may also be estimated with the \code{plm}
  function if \code{np=TRUE}.
}
\seealso{
  \code{\link{plm}} and \code{pooltest} for poolability tests.
}
\examples{
library(Ecdat)
data(Produc)
Produc <-pdata.frame(Produc,state,year)
zz <- nopool(log(gsp)~log(pcap)+log(pc)+log(emp)+unemp,data=Produc)
zz <- nopool(log(gsp)~log(pcap)+log(pc)+log(emp)+unemp,effect="temp",data=Produc)
summary(zz)
plot(zz)
}
\keyword{regression}
