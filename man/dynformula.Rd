\name{dynformula}
\alias{dynformula}
\title{Dynamic formula}
\description{
  A function to create easily formula with lags
}
\usage{
dynformula(formula,lag.form=NULL,diff.form=NULL,log.form=NULL)
}
\arguments{
  \item{formula}{a formula,}
  \item{lag.form}{a list containing the lag structure of each variable
    in the formula,}
  \item{diff.form}{a list of logical values indicating whether variables 
    should be differenced,}
  \item{log.form}{a list of logical values indicating whether variables should
    be in logarithms.}
}

\value{
  an object of class \code{c("dynformula","formula")}, which is a
  formula with four additional attributes : \code{var}, the names of the variables in
  the formula, \code{lag}, \code{diff} and \code{log} which store the
  information about lags, differences and logs.
}
\details{
  \code{lag.form}, \code{diff.form} and \code{log.form} should be lists of length equal to the total number
  variables. Each  element of these lists is :

  either a vector of length 2 (c(1,4) means
  lags 1,2,3 and 4) or a scalar (3 means lags 0,1,2,3 except for
  the lhs variable for which it is 1,2,3) for \code{lag.form}.

  logical values for \code{diff.form} and \code{log.form}.

  It can also be an incomplete named list (to apply the transformation
  for only some variables) with eventually an unnamed element which then
  is the default value.
  }
  
\author{Yves Croissant}

\examples{
# all variables in log, x1, x2 and x3 laged twice, y laged once and x3 differenced
dynformula(y~x1+x2+x3,lag.form=list(2,y=1),diff.form=list(x3=TRUE),log.form=TRUE)
}
\keyword{regression}
