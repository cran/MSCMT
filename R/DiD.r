#' Difference-in-difference estimator based on SCM
#'
#' \code{did} calculates difference-in-difference estimators based on SCM.
#'
#' \code{did} calculates difference-in-difference estimators with corresponding
#' p-values (if results of a placebo study are present) based on the Synthetic 
#' Control Method.
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
#' @param exclude.ratio A numerical scalar (default: \code{Inf}). When 
#' calculating the p-value, control units with an average pre-treatment gap
#' of more then \code{exclude.ratio} times the average pre-treatment gap of
#' the treated unit are excluded from the analysis.
#' @return A list with components \code{effect.size}, \code{average.pre} and
#' \code{everage.post}. If \code{x} contains the results of a placebo study,
#' three components \code{p.value}, \code{rank}, and \code{excluded} (with the 
#' names of the excluded units) are included additionally.
#' @export did
#' @examples
#' \dontrun{
#' ## for an example, see the main package vignette:
#'  vignette("WorkingWithMSCMT",package="MSCMT")
#' }
did <- function(x,what,range.pre,range.post,
                alternative = c("two.sided", "less", "greater"),
                exclude.ratio = Inf) {
  alternative <- match.arg(alternative)  
  if (has.placebo <- !is.null(x$placebo)) {
    if (missing(what)) what  <- x[[1]]$dependent[1]
    gaps <- x$placebo[[what]]$gaps
    if (missing(range.pre))  
      if (what %in% x[[1]]$dependent) range.pre <- x[[1]]$times.dep[,what] else
                                      stop("range.pre is missing")
  } else {
    if (missing(what)) what  <- x$dependent[1]
    gaps <- x$gaps[[what]]
    if (missing(range.pre))  
      if (what %in% x$dependent) range.pre <- x$times.dep[,what] else
                                 stop("range.pre is missing")
  }
  gaps.pre  <- AQMwindow(gaps,range.pre)
  gaps.post <- if (missing(range.post)) AQMtail(gaps,range.pre[2]) else
                                        AQMwindow(gaps,range.post)
  if (has.placebo) {
    D.pre     <- apply(gaps.pre,2,function(x) mean(x))
    D.post    <- apply(gaps.post,2,function(x) mean(x))
    excl      <- (abs(D.pre) > exclude.ratio*abs(D.pre[1]))
    effect    <- as.numeric(D.post - D.pre)[!excl]
    p.value   <- switch(alternative,
                   two.sided = mean(abs(effect[1])<=abs(effect)),
                   less      = mean(effect[1]>=effect),
                   greater   = mean(effect[1]<=effect))
    list(effect.size=effect[1],average.pre=as.numeric(D.pre[1]),
         average.post=as.numeric(D.post[1]),p.value=p.value,
         rank=sum(effect[1]>=effect),excluded=names(excl)[excl])
  } else list(effect.size=as.numeric(mean(gaps.post) - mean(gaps.pre)),
              average.pre=as.numeric(mean(gaps.pre)),
              average.post=as.numeric(mean(gaps.post)))
}
