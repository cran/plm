\name{sargan}
\alias{sargan}
\title{Hansen--Sargan test of overidentifying restrictions}
\description{
   A test of overidentifying restrictions for GMM models
}
\usage{
sargan(object)
}

\arguments{
  \item{object}{an object of class \code{"pgmm"}.}
}

\value{
  An object of class \code{"htest"}.
}

\details{
  The Hansen--Sargan test calculate the quadratic form of the moment
  restrictions that is  minimized while computing the GMM
  estimator. It follows asymptotically a chi square with a number of
  degrees of freedom equal to the difference of the number of moment conditions and
  the number of coefficients.
}
\references{
  Hansen, L.P. (1982), Large Sample Properties of Generalized Methods of
  Moments Estimators, \emph{Econometrica}, 50, 1029--54.
  
  Sargan, J.D. (1958), The Estimation of Economic Relationships using
  Instrumental Variables, \emph{Econometrica}, 26, pp.393--415.
}

\author{Yves Croissant}
  
\examples{
data(EmplUK)
pdata.frame(EmplUK,"firm","year")
ar <- pgmm(dynformula(log(emp)~log(wage)+log(capital)+log(output),list(2,1,2,2)),EmplUK,
           effect="twoways",model="twosteps",
           gmm.inst=~log(emp),lag.gmm=list(c(2,99)))
sargan(ar)
}
\keyword{htest}
