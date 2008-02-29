\name{pwartest}
\alias{pwartest}
\title{Wooldridge Test for AR(1) Errors in FE Panel Models}
\description{

  Test of serial correlation for (the idiosyncratic component of) the
  errors in fixed-effects panel models.
  
}
\usage{
pwartest(x, data, index = NULL, ...)
}

\arguments{
  \item{x}{an object of class \code{formula},}
  \item{data}{an object of class \code{pdata.frame}}
  \item{index}{the index of the \code{data.frame},}
  \item{...}{further arguments to be passed on to \code{linear.hypothesis} or to \code{pvcovHC}}
}

\value{
  An object of class \code{"htest"}.
}

\details{

  As Wooldridge (2003, 10.5.4) observes, under the null of no serial
  correlation in the errors, the residuals of a FE model must be
  negatively serially correlated, with $cor(\hat{u}_{it},
  \hat{u}_{is})=-1/(T-1)$ for each $t,s$ . He suggests basing a test for
  this null hypothesis on a pooled regression of FE residuals on their
  first lag: $\hat{u}_{i,t}=\alpha + \delta \hat{u}_{i,t-1} +
  \eta_{i,t}$. Rejecting the restriction $\delta = -1/(T-1)$ makes us
  conclude against the original null of no serial correlation.

\code{pwartest} estimates the \code{within} model and retrieves
residuals, then estimates an AR(1) \code{pooling} model on them. The
test statistic is obtained by applying \code{linear.hypothesis()} to the
latter model to test the above restriction on $\delta$, setting the
covariance matrix to \code{pvcovHC} with the option
\code{method="arellano"} to control for serial correlation.

Unlike the \code{\link{pbgtest}} and \code{\link{pdwtest}}, this test does not rely on
T-asymptotics and has therefore good properties in ``short''
panels. Furthermore, it is robust to general heteroskedasticity.

}
\references{
  Wooldridge, J.M. (2002) \emph{Econometric analysis of cross-section and panel data},
    MIT Press, 10.5.4, page 274.
}
\seealso{
    \code{\link{pdwtest}}, \code{\link{pbgtest}}
}
% \examples{
% data("EmplUK", package="plm")
% pwartest(log(emp) ~ log(wage) + log(capital), data=EmplUK)
% }
\keyword{htest}
