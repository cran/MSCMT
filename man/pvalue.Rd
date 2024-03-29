% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pvalue.r
\name{pvalue}
\alias{pvalue}
\title{P-values for placebo studies}
\usage{
pvalue(
  x,
  what,
  range.pre,
  range.post,
  alternative = c("two.sided", "less", "greater"),
  exclude.ratio = Inf,
  ratio.type = c("rmspe", "mspe")
)
}
\arguments{
\item{x}{An object of class \code{"mscmt"}, usually obtained as
the result of a call to function \code{\link{mscmt}}.}

\item{what}{A character vector. Name of the variable to be considered. If 
missing, the (first) dependent variable will be used.}

\item{range.pre}{A vector of length 2 defining the range of the pre-treatment
period with start and end time given as 
\itemize{
\item annual dates, if the format of start/end time is "dddd", e.g. "2016",
\item quarterly dates, if the format of start/end time is "ddddQd", e.g. 
"2016Q1",
\item monthly dates, if the format of start/end time is "dddd?dd" with "?"
different from "W" (see below), e.g. "2016/03" or "2016-10",
\item weekly dates, if the format of start/end time is "ddddWdd", e.g. 
"2016W23",
\item daily dates, if the format of start/end time is "dddd-dd-dd", e.g. 
"2016-08-18",
}
corresponding to the format of the respective column of the \code{times.dep}
argument of \code{\link{mscmt}}.
If missing, the corresponding column of \code{times.dep} will be used.}

\item{range.post}{A vector of length 2 defining the range of the 
post-treatment period with start and end time given as 
\itemize{
\item annual dates, if the format of start/end time is "dddd", e.g. "2016",
\item quarterly dates, if the format of start/end time is "ddddQd", e.g. 
"2016Q1",
\item monthly dates, if the format of start/end time is "dddd?dd" with "?"
different from "W" (see below), e.g. "2016/03" or "2016-10",
\item weekly dates, if the format of start/end time is "ddddWdd", e.g. 
"2016W23",
\item daily dates, if the format of start/end time is "dddd-dd-dd", e.g. 
"2016-08-18",
}
corresponding to the format of the respective column of the \code{times.dep}
argument of \code{\link{mscmt}}. Will be guessed if missing.}

\item{alternative}{A character string giving the alternative of the test. 
Either \code{"two.sided"} (default), \code{"less"}, or \code{"greater"}.}

\item{exclude.ratio}{A numerical scalar (default: \code{Inf}). Control units
with a pre-treatment-(r)mspe of more than \code{exclude.ratio} times the
pre-treatment-(r)mspe of the treated unit are excluded from the calculations
of the p-value.}

\item{ratio.type}{A character string. Either \code{rmspe} (default) or 
\code{mspe}. 
Selects whether root mean squared errors or mean squared errors are 
calculated.}
}
\value{
A time series containing the p-values for the 
post-treatment periods.
}
\description{
\code{pvalue} calculates p-values for placebo studies.
}
\details{
\code{pvalue} calculates p-values for placebo studies based
on Synthetic Control Methods.
}
\examples{
\dontrun{
## for an example, see the main package vignette:
 vignette("WorkingWithMSCMT",package="MSCMT")
}
}
