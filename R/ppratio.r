#' Post-pre-(r)mspe-ratios for placebo studies
#'
#' \code{ppratio} calculates post-to-pre-(r)mspe-ratios for placebo studies.
#'
#' \code{ppratio} calculates post-to-pre-(r)mspe-ratios for placebo studies based
#' on Synthetic Control Methods.
#'
#' @param x An object of class \code{"mscmt"}, usually obtained as
#' the result of a call to function \code{\link{mscmt}}.
#' @param what A character vector. Name of the variable to be considered. If 
#' missing, the (first) dependent variable will be used.
##' @param range.pre A vector of length 2 defining the range of the pre-treatment
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
#' @param type A character string. Either \code{rmspe} (default) or \code{mspe}.
#' Selects whether root mean squared errors or mean squared errors are 
#' calculated.
#' @param return.all A logical scalar. If \code{FALSE} (default), only the
#' (named) vector of post-pre-(r)mspe-ratios is returned, if \code{TRUE},
#' a three-column matrix with pre- and post-treatment (r)mspe's as well as the 
#' post-pre-ratios will be returned.
#' @return If \code{return.all} is \code{FALSE}, a (named) vector of 
#' post-pre-(r)mspe-ratios. If \code{return.all} is \code{TRUE}, a matrix with
#' three columns containing the pre-treatment (r)mspe, the post-treatment 
#' (r)mspe, and the post-pre-ratio.
#' @export ppratio
ppratio <- function(x,what,range.pre,range.post,type=c("rmspe","mspe"),
                    return.all=FALSE) {
  type <- match.arg(type)                  
  if (is.null(x$placebo)) stop("results of placebo study are missing")
  if (missing(what)) what  <- x[[1]]$dependent[1]
  gaps <- x$placebo[[what]]$gaps
  if (missing(range.pre))  
    if (what %in% x[[1]]$dependent) range.pre <- x[[1]]$times.dep[,what] else
                                    stop("range.pre is missing")
  gaps.pre   <- AQMwindow(gaps,range.pre)
  gaps.post <- if (missing(range.post)) AQMtail(gaps,range.pre[2]) else
                                        AQMwindow(gaps,range.post)
  rmspe.pre  <- apply(gaps.pre,2,function(x) mean(x^2))
  rmspe.post <- apply(gaps.post,2,function(x) mean(x^2))
  if (type=="rmspe") {
    rmspe.pre  <- sqrt(rmspe.pre)
    rmspe.post <- sqrt(rmspe.post)
  }
  if (return.all) 
    cbind(pre=rmspe.pre,post=rmspe.post,ratio=rmspe.post/rmspe.pre) else
    rmspe.post/rmspe.pre
}
