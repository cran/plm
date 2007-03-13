\name{pread.table}
\alias{pread.table}

\title{Read panel data from text file}
\description{
  Create a \code{pdata.frame} directly by reading a text file.
}
\usage{
 pread.table(file,id,time=NULL,name=NULL,...)

}
\arguments{
  \item{file}{a text file,}
  \item{id}{the name of the individual index,}
  \item{time}{the name of the time index,}
  \item{name}{the name of the \code{pdata.frame}. If missing, the name of the
    \code{pdata.frame} is the name of the file, without the path and the
  extension,}
  \item{...}{further arguments passed to \code{read.table}.}
}

\value{
  a \code{pdata.frame} object
}


\author{Yves Croissant}

\seealso{
  \code{\link{pdata.frame}}
}

\keyword{attribute}
