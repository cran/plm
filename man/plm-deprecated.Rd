% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{plm-deprecated}
\alias{plm-deprecated}
\alias{detect_lin_dep}
\alias{pvcovHC}
\alias{plm.data}
\alias{dynformula}
\alias{formula.dynformula}
\alias{print.dynformula}
\title{Deprecated functions of plm}
\usage{
pvcovHC(x, ...)

plm.data(x, indexes = NULL)

dynformula(formula, lag.form = NULL, diff.form = NULL, log.form = NULL)

\method{formula}{dynformula}(x, ...)

\method{print}{dynformula}(x, ...)
}
\arguments{
\item{x}{an object of class \code{"plm"},}

\item{\dots}{further arguments.}

\item{indexes}{a vector (of length one or two) indicating the (individual
and time) indexes (see Details);}

\item{formula}{a formula,}

\item{lag.form}{a list containing the lag structure of each variable in the
formula,}

\item{diff.form}{a vector (or a list) of logical values indicating whether
variables should be differenced,}

\item{log.form}{a vector (or a list) of logical values indicating whether
variables should be in logarithms.}

\item{data}{a \code{data.frame},}
}
\description{
\code{dynformula}, \code{pht}, \code{plm.data}, and \code{pvcovHC} are
deprecated functions which could be removed from \pkg{plm} in a near future.
}
\details{
\code{dynformula} was used to construct a dynamic formula which was the
first argument of \code{pgmm}. \code{pgmm} uses now multi-part formulas.

\code{pht} estimates the Hausman-Taylor model, which can now be estimated
using the more general \code{plm} function.

\code{plm.data} is replaced by \code{pdata.frame}.

\code{pvcovHC} is replaced by \code{vcovHC}.

\code{detect_lin_dep} was renamed to \code{detect.lindep}.
}
