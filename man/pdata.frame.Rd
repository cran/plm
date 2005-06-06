\name{pdata.frame}
\alias{pdata.frame}
\alias{niTt}
\alias{indexes}
\alias{balanced}
\alias{nT}
\alias{na.obs}
\alias{summary.panel.frame}
\alias{print.summary.panel.frame}

\title{Data frame for panel data}
\description{
  a \code{data.frame} with further attributes usefull for panel data
}
\usage{
x <- pdata.frame(x,id,time=NULL)
}
\arguments{
  \item{x}{a \code{data.frame}}
  \item{id}{the individual index}
  \item{time}{the time index}
}

\value{
  An object of class \code{c("pdata.frame","data.frame")}, which is a
  \code{data.frame} with an additional attribute \code{panel}. This
  attribute is a list containing :
  \item{indexes}{a list containing the names of the
    individual (\code{id}) and time (\code{time}) index variables,}
  \item{nT}{a list containing the number of individuals (\code{n}), of time
    observations (\code{T}) and of total observations (\code{N}),}
  \item{balanced}{\code{TRUE} for a balanced panel,}
  \item{Tint}{a list of two vectors containing the number of time observation for
    each individual (\code{Ti}) and the number of individual observation for each
    time observation (\code{nt}),}
  \item{panel.names}{a list of two character vectors containing the
    labels of the individuals (\code{id.names}) and of the time
    observations (\code{time.names}),}
  \item{dim.var}{a list of two logical vectors indicating wheter each
    variable presents variations accross individuals (\code{id.var}) and
    accross time periods (\code{time.var}).}
  These attributes may be extracted using synonimus functions.
}

\details{
The observations are assumed to be sorted by individuals first and by
time observations. \code{id} and \code{time} can be the names of the
indexes contained in the \code{data.frame} if any. \code{time} may be
ommited, in this case a new variable called \code{"time"} is added to
the \code{data.frame}. For balanced pannels, the argument \code{id} can
be an integer. In this case, two new variables, \code{"time"} and
\code{"id"} are added to the \code{data.frame}.

\code{pdata.frame} calls \code{pdim} to check for the individual and
time dimensions of the panel and \code{pvarcheck} to check for the
individual and time variation of each variable.

\code{pdata.frame} adds an attribute \code{"data"} to each variable in
the \code{data.frame}. This can be usefull to use simply the special
panel operators in model formulas. This works properly if the result of
the \code{pdata.frame} is stored in an object of the same name as its
\code{data.frame} argument.
  
  
The \code{print} method is inheritaded from \code{data.frame}. A
\code{summary} method is provided which add to the
\code{summary.data.frame} method informations about the structure of the panel.
}
\seealso{
  \code{\link{plm}}
}
\examples{
library(Ecdat)

#There are 595 individuals

data(Wages)
Wages <- pdata.frame(Wages,595)

# Gasoline contains two variables which are individual and time indexes
data(Gasoline)
Gasoline <- pdata.frame(Gasoline,country,year)
summary(Gasoline)

# Hedonic is an unbalanced panel, townid is the individual index
data(Hedonic)
Hedonic=pdata.frame(Hedonic,townid)

}
\keyword{attribute}
