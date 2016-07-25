#' Plotting Results of MSCMT
#'
#' \code{plot.mscmt} plots results of \code{mscmt}.
#'
#' A unified basic plot function for gaps plots, comparison of treated and 
#' synthetic values, as well as plots for placebo studies.
#' Consider using \code{\link{ggplot.mscmt}} instead, which is the preferred
#' plot method and has more functionality than \code{\link{plot.mscmt}}.
#'
#' @param x An object of class \code{"mscmt"}, usually obtained as
#' the result of a call to function \code{\link{mscmt}}.
#' @param what A character scalar. Name of the variable to be plotted. If 
#' missing, the (first) dependent variable will be used.
#' @param type A character scalar denoting the type of the plot containing 
#' either \code{"gaps"}, \code{"comparison"}, \code{"placebo.gaps"}, 
#' or \code{"placebo.data"}. 
#' Partial matching allowed, defaults to \code{"placebo.gaps"}, if results
#' of a placebo study are present, and to \code{"gaps"}, else.
#' @param treatment.time An optional numerical scalar. If not missing, a 
#' vertical dotted line at the given point in time is included in the plot. 
#' \code{treatment.time} is measured in years, but may as well be a decimal 
#' number to reflect treatment times different from January 1st.
#' @param zero.line A logical scalar. If \code{TRUE} (default), a horizontal
#' dotted line (at zero level) is plotted for \code{"gaps"} and \code{"placebo"}
#' plots.
#' @param ylab Optional label for the y-axis, automatically generated if 
#' missing.
#' @param xlab Optional label for the x-axis, defaults to \code{"Date"}.
#' @param main Optional main title for the plot, automatically generated 
#' if missing.
#' @param sub Optional subtitle for the plot. If missing, the subtitle is 
#' generated automatically for \code{"comparison"} and \code{"gaps"} plots.
#' @param col Optional character vector with length corresponding to the number
#' of units. Contains the colours for the different units, automatically 
#' generated if missing. 
#' @param lty Optional numerical vector with length corresponding to the number
#' of units. Contains the line types for the different units, automatically 
#' generated if missing.
#' @param lwd Optional numerical vector with length corresponding to the number
#' of units. Contains the line widths for the different units, automatically 
#' generated if missing.
#' @param legend A logical scalar. If \code{TRUE} (default), a legend is 
#' included in the plot.
#' @param bw A logical scalar. If \code{FALSE} (default), the automatically 
#' generated colours and line types are optimized for a colour plot, 
#' if \code{TRUE}, the automatic colours and line types are set for a black and
#' white plot.
#' @param ... Further optional parameters for the underlying 
#' \code{\link[graphics]{plot}} function.
#' @return Nothing useful (function is called for its side effects).
#' @importFrom graphics plot abline
#' @method plot mscmt
#' @export
plot.mscmt <- function(x,what,type=c("gaps","comparison","placebo.gaps",
                       "placebo.data"),
                       treatment.time,zero.line=TRUE,ylab,xlab="Date",main,sub,
                       col,lty,lwd,legend=TRUE,bw=FALSE,...) {
  if (missing(type)&&(!is.null(x$placebo))) type <- "placebo.gaps"
  type <- match.arg(type)                      
  if (missing(what)) {
    if (is.null(x$placebo)) what <- x$dependent[1] else
                            what <- x[[1]]$dependent[1]
  }  
  if ((type=="placebo.gaps")||(type=="placebo.data")) {
    if (is.null(x$placebo)) stop("results of placebo study missing")
    if (!(what %in% names(x$placebo)))
      stop(paste("variable",what,"missing in results"))
    unames <- colnames(x$placebo[[what]]$gaps)
    if (missing(col)&&bw) col=rep("black",length(unames))
    if (missing(col)&&(!bw)) col=c("red",rep("black",length(unames)-1))
    if (missing(lty)&&bw) lty=c(1,rep(2,length(unames)-1))
    if (missing(lty)&&(!bw)) lty=rep(1,length(unames))
    if (missing(lwd)) lwd=c(2,rep(1,length(unames)-1))
    if (missing(ylab)) ylab=paste(if (type=="placebo.gaps") "Gaps" else "Data",
                                  "for",what)
    if (missing(main)) main=paste("Placebo plot for",what)
    if (missing(sub)) sub=""
    plot(if (type=="placebo.gaps") x$placebo[[what]]$gaps else 
                                   x$placebo[[what]]$data.treat,
         xlab=xlab,ylab=ylab,lty=lty,col=col,
         plot.type="single",main=main,sub=sub,lwd=lwd,...)
    if (!missing(treatment.time)) abline(v=treatment.time,lty=3)
    if (zero.line) abline(h=0,lty=3)         
    if (legend) legend("topleft",legend=c(unames[1],"other"),col=col[1:2],
                       lty=lty[1:2],lwd=lwd[1:2])
  }
  if (type=="comparison") {
    if (is.null(x$combined)) stop("input is not an individual mscmt result")
    if (!(what %in% names(x$combined)))
      stop(paste("variable",what,"missing in results"))
    if (missing(col)&&bw) col=rep("black",2)
    if (missing(col)&&(!bw)) col=c("black","red")
    if (missing(lty)&&bw) lty=c(1,2)
    if (missing(lty)&&(!bw)) lty=rep(1,2)
    if (missing(lwd)) lwd=c(2,2)
    if (missing(ylab)) ylab=what
    if (missing(main)) main=paste("Comparison of",what)
    if (missing(sub)) sub=paste("treated unit: ",x$treated.unit)
    plot(x$combined[[what]][,1:2],xlab=xlab,ylab=ylab,lty=lty,col=col,
         plot.type="single",main=main,sub=sub,lwd=lwd,...)
    if (!missing(treatment.time)) abline(v=treatment.time,lty=3)
    if (legend) legend("topleft",legend=c("actual data","synthesized data"),
                       col=col[1:2],lty=lty[1:2],lwd=lwd[1:2])
  }
  if (type=="gaps") {
    if (is.null(x$combined)) stop("input is not an individual mscmt result")
    if (!(what %in% names(x$combined)))
      stop(paste("variable",what,"missing in results"))
    if (missing(col)&&bw) col="black"
    if (missing(col)&&(!bw)) col="black"
    if (missing(lty)&&bw) lty=1
    if (missing(lty)&&(!bw)) lty=1
    if (missing(lwd)) lwd=2
    if (missing(ylab)) ylab=what
    if (missing(main)) main=paste("Gaps for",what)
    if (missing(sub)) sub=paste("treated unit: ",x$treated.unit)
    plot(x$combined[[what]][,3],xlab=xlab,ylab=ylab,lty=lty,col=col,
         plot.type="single",main=main,sub=sub,lwd=lwd,...)
    if (!missing(treatment.time)) abline(v=treatment.time,lty=3)
    if (zero.line) abline(h=0,lty=3)         
  }
}
