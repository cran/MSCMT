#' P-values for placebo studies
#'
#' \code{pvalue} calculates p-values for placebo studies.
#'
#' \code{pvalue} calculates p-values for placebo studies based
#' on Synthetic Control Methods.
#'
#' @param x An object of class \code{"mscmt"}, usually obtained as
#' the result of a call to function \code{\link{mscmt}}.
#' @param what A character vector. Name of the variable to be considered. If 
#' missing, the (first) dependent variable will be used.
#' @param range.pre A vector of length 2 defining the range of the pre-treatment
#' period with start and end time given as 
#' \itemize{
#' \item annual dates, if the format of start/end time is "dddd", e.g. "2016",
#' \item quarterly dates, if the format of start/end time is "ddddQd", e.g. 
#' "2016Q1",
#' \item monthly dates, if the format of start/end time is "dddd?dd", e.g. 
#' "2016/03" or "2016-10",
#' }
#' corresponding to the format of the respective column of the \code{times.dep}
#' argument of \code{\link{mscmt}}.
#' If missing, the corresponding column of \code{times.dep} will be used.
#' @param range.post A vector of length 2 defining the range of the 
#' post-treatment period with start and end time given as 
#' \itemize{
#' \item annual dates, if the format of start/end time is "dddd", e.g. "2016",
#' \item quarterly dates, if the format of start/end time is "ddddQd", e.g. 
#' "2016Q1",
#' \item monthly dates, if the format of start/end time is "dddd?dd", e.g. 
#' "2016/03" or "2016-10",
#' }
#' corresponding to the format of the respective column of the \code{times.dep}
#' argument of \code{\link{mscmt}}. Will be guessed if missing.
#' @param alternative A character string giving the alternative of the test. 
#' Either \code{"two.sided"} (default), \code{"less"}, or \code{"greater"}.
#' @param exclude.ratio A numerical scalar (default: \code{Inf}). Control units
#' with a pre-treatment-(r)mspe of more than \code{exclude.ratio} times the
#' pre-treatment-(r)mspe of the treated unit are excluded from the calculations
#' of the p-value.
#' @param ratio.type A character string. Either \code{rmspe} (default) or 
#' \code{mspe}. 
#' Selects whether root mean squared errors or mean squared errors are 
#' calculated.
#' @return A time series containing the p-values for the 
#' post-treatment periods.
#' @export pvalue
#' @examples
#' \dontrun{
#' ## for an example, see the main package vignette:
#'  vignette("WorkingWithMSCMT",package="MSCMT")
#' }
pvalue <- function(x,what,range.pre,range.post,
                   alternative = c("two.sided", "less", "greater"),
                   exclude.ratio = Inf, ratio.type = c("rmspe","mspe")) {
  ratio.type <- match.arg(ratio.type)
  if (is.null(x$placebo)) stop("results of placebo study are missing")
  alternative <- match.arg(alternative)  
  if (missing(what)) what  <- x[[1]]$dependent[1]
  gaps <- x$placebo[[what]]$gaps
  if (missing(range.pre))  
    if (what %in% x[[1]]$dependent) range.pre <- x[[1]]$times.dep[,what] else
                                    stop("range.pre is missing")
  gaps.pre  <- AQMwindow(gaps,range.pre)
  gaps.post <- if (missing(range.post)) AQMtail(gaps,range.pre[2]) else
                                        AQMwindow(gaps,range.post)
  rmspe.pre <- apply(gaps.pre,2,function(x) mean(x^2))
  if (ratio.type=="rmspe") rmspe.pre <- sqrt(rmspe.pre)
  excl <- (rmspe.pre > exclude.ratio*rmspe.pre[1])
  gaps.post <- gaps.post[,!excl,drop=FALSE]
  res  <- apply(gaps.post,1,switch(alternative,
                  two.sided = function(x) mean(abs(x[1])<=abs(x)),
                  less      = function(x) mean(x[1]>=x),
                  greater   = function(x) mean(x[1]<=x)))
  ts(res,start=start(gaps.post),frequency=frequency(gaps.post))                  
}
