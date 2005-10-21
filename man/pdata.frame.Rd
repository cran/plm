\name{pdata.frame}
\alias{pdata.frame}
\alias{summary.pdata.frame}
\alias{print.summary.pdata.frame}

\title{Data frame for panel data}
\description{
  a \code{data.frame} with further attributes useful for panel data
}
\usage{
 pdata.frame(x,id,time=NULL,name=NULL)
\method{summary}{pdata.frame}(object,...)
\method{print}{summary.pdata.frame}(x,...)
}
\arguments{
  \item{x}{a \code{data.frame} (a \code{pdata.frame} for the print
    method),}
  \item{object}{a \code{pdata.frame} object,}
  \item{id}{the individual index}
  \item{time}{the time index}
  \item{name}{the name of the \code{pdata.frame} (the name of \code{x} if \code{NULL}),}
  \item{...}{further arguments.}
}

\value{
  This function returns internally an object of class \code{c("pdata.frame","data.frame")}, which is a
  \code{data.frame} with three additional attributes :
  \item{indexes}{a list containing the names of the
    individual (\code{id}) and time (\code{time}) index variables,}
  \item{pdim}{a list containing the dimension of the \code{pdata.frame} (see \code{\link{pdim}}),}
  \item{pvar}{a list of two logical vectors indicating whether each
    variable has variations across individuals (\code{id.var}) and
    accross time periods (\code{time.var}) (see \code{\link{pvar}}).}
}

\details{
The observations are assumed to be sorted by individuals first and by
time observations. \code{id} and \code{time} can be the names of the
indexes contained in the \code{data.frame} if any. \code{time} may be
omitted, in this case a new variable called \code{"time"} is added to
the \code{data.frame}. For balanced panels, the argument \code{id} can
be an integer. In this case, two new variables, \code{"time"} and
\code{"id"} are added to the \code{data.frame}.

\code{pdata.frame} calls \code{pdim} to check for the individual and
time dimensions of the panel and \code{pvar} to check for the
individual and time variation of each variable.

All elements of a \code{pdata.frame} are of class \code{pserie} and have
an attribute \code{data} which is the name of the \code{pdata.frame}.
  
The \code{print} method is inherited from \code{data.frame}. A
\code{summary} method is provided which add to the
\code{summary.data.frame} method informations about the structure of the panel.
}
\seealso{
  \code{\link{plm}}, \code{\link{pdim}}, \code{\link{pvar}}, \code{\link{pread.table}}.
}
\examples{
library(Ecdat)

#There are 595 individuals

data(Wages)
pdata.frame(Wages,595)

# Gasoline contains two variables which are individual and time indexes
# The pdata.frame is called gas
data(Gasoline)
pdata.frame(Gasoline,"country","year",name="gas")
summary(Gasoline)

# Hedonic is an unbalanced panel, townid is the individual index
data(Hedonic)
pdata.frame(Hedonic,"townid")

}
\keyword{attribute}
