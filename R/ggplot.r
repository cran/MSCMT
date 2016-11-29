#' Plotting Results of mscmt with ggplot2
#'
#' \code{ggplot.mscmt} plots results of \code{\link{mscmt}} based on 
#' \code{\link[ggplot2]{ggplot}}.
#'
#' A unified plot method for gaps plots, comparison of treated and synthetic
#' values, as well as plots for placebo studies based on 
#' \code{\link[ggplot2]{ggplot}}. \code{\link{ggplot.mscmt}} is the preferred
#' plot method and has more functionality than \code{\link{plot.mscmt}}.
#'
#' @param x An object of class \code{"mscmt"}, usually obtained as
#' the result of a call to function \code{\link{mscmt}}.
#' @param what A character vector. Name(s) of the variables to be plotted. If 
#' missing, the (first) dependent variable will be used.
#' @param type A character scalar denoting the type of the plot containing 
#' either \code{"gaps"}, \code{"comparison"}, \code{"placebo.gaps"}, 
#' \code{"placebo.data"}, or \code{"p.value"}. 
#' Partial matching allowed, defaults to \code{"placebo.gaps"}, if results
#' of a placebo study are present, and to \code{"gaps"}, else.
#' @param treatment.time An optional scalar (numeric, character, or 
#' \code{\link[base]{Date}}) giving the treatment time. 
#' If \code{treatment.time} is numeric, Jan 01 of that year will be used. If 
#' \code{treatment.time} is a character string, it will be converted to a 
#' \code{\link[base]{Date}} and must thus be in an unambiguous format.
#' A vertical dotted line at the given point in time is included in the plot. 
#' @param zero.line A logical scalar. If \code{TRUE} (default), a horizontal
#' dotted line (at zero level) is plotted for \code{"gaps"} and 
#' \code{"placebo.gaps"} plots.
#' @param ylab Optional label for the y-axis, automatically generated if 
#' missing.
#' @param xlab Optional label for the x-axis, defaults to \code{"Date"}.
#' @param main Optional main title for the plot, automatically generated 
#' if missing.
#' @param col Optional character vector with length 1 (for gaps plots) or
#' 2 (for all other plot types). For comparison plots, \code{col} contains the
#' colours for the actual and synthesized data, for placebo.plots (with 
#' \code{full.legend==FALSE}), \code{col} contains the colours for 
#' the treated unit and the control units. Automatically generated 
#' if missing. 
#' @param lty Optional numerical vector with length 1 (for gaps plots) or
#' 2 (for all other plot types). For comparison plots, \code{lty} contains the
#' linetypes for the actual and synthesized data, for placebo.plots (with 
#' \code{full.legend==FALSE}), \code{col} contains the linetypes for 
#' the treated unit and the control units. Automatically generated 
#' if missing. 
#' @param lwd Optional numerical vector with length 1 (for gaps plots) or
#' 2 (for all other plot types). For comparison plots, \code{lty} contains the
#' linewidths for the actual and synthesized data, for placebo.plots (with 
#' \code{full.legend==FALSE}), \code{col} contains the linewidths for 
#' the treated unit and the control units. Automatically generated 
#' if missing. 
#' @param legend A logical scalar. If \code{TRUE} (default), a legend is 
#' included in the plot.
#' @param bw A logical scalar. If \code{FALSE} (default), the automatically 
#' generated colours and line types are optimized for a colour plot, 
#' if \code{TRUE}, the automatic colours and line types are set for a black and
#' white plot.
#' @param date.format A character string giving the format for the tick labels
#' of the x axis as documented in \code{\link[base]{strptime}}. Defaults to
#' \code{"\%b \%y"} or \code{"\%Y"}, depending on the granularity of the data.
#' @param unit.name A character string with the title of the legend 
#' for comparison and placebo plots. Defaults to "Estimation" for comparison 
#' and "Unit" for placebo plots.
#' @param full.legend A logical scalar. If \code{TRUE} (default), a full legend
#' of all units (donors) is constructed. If \code{FALSE}, only the treated and
#' the control units are distinguished.
#' @param include.smooth A logical scalar. If \code{TRUE}, a geometric smoother
#' based on the control units is added to placebo plots. Default: \code{FALSE}.
#' @param include.mean A logical scalar. If \code{TRUE}, the arithmetic mean
#' of all control units is added to placebo plots. Default: \code{FALSE}.
#' @param include.synth A logical scalar. If \code{TRUE}, the synthesized data
#' for the treated unit are added to plots of type \code{"placebo.data"}.
#' Defaults to \code{FALSE}.
#' @param draw.estwindow A logical scalar. If \code{TRUE} (default), the time
#' range containing all optimization periods is shaded in the corresponding
#' plots.
#' @param what.set An optional character string for a convenient selection of 
#' multiple variables. Accepted values are \code{"dependents"}, 
#' \code{"predictors"}, and \code{"all"}, which collects all dependent, all
#' predictor, or all variables of both types, respectively. Overrides parameter
#' \code{what} (if the latter is present).
#' @param limits An optional vector of length 2 giving the range of the plot or
#' \code{NULL}.
#' If \code{limits} is numeric, Jan 01 of the corresponding years will be used.
#' If \code{limits} is of type character, both strings will be converted to 
#' Dates (via \code{\link[base]{as.Date}}) and must thus be in an unambiguous 
#' format.
#' @param alpha Either a numerical scalar, a numerical vector of length 
#' corresponding to the number of units, or the character string \code{"auto"}.
#' If \code{alpha} is a numerical scalar (default with value \code{1}), a fixed
#' value for the alpha channel (transparency) is included for all units in 
#' placebo plots. If \code{alpha} is numeric and has length corresponding to the
#' number of units, these values are assigned as alpha channel to the individual 
#' units. If \code{"auto"}, the alpha channel information is obtained from the
#' w weights of the control units.
#' @param alpha.min A numerical scalar (default: \code{0.1}). If \code{alpha} is
#' set to \code{"auto"}, the individual alpha channel information for control 
#' unit \code{i} is set to \code{alpha.min + (1-alpha.min) * w[i]}.
#' @param exclude.units An optional (default: \code{NULL}) character vector with 
#' names for control units which shall be excluded from placebo plots and
#' p-value calculations.
#' @param exclude.ratio A numeric scalar (default: \code{Inf}). Control units 
#' with a pre-treatment (r)mspe of more than \code{exclude.ratio} times
#' the pre-treatment (r)mspe of the treated unit are excluded from placebo 
#' plots and p-value calculations.
#' @param ratio.type A character string. Either \code{rmspe} (default) or 
#' \code{mspe}. Selects whether root mean squared errors or mean squared errors 
#' are considered for the exclusion of control units (see \code{exclude.ratio}).
#' @param alternative A character string giving the alternative of the test for 
#' plots of type \code{"p.value"}. Either \code{"two.sided"} (default), 
#' \code{"less"}, or \code{"greater"}.
#' @param draw.points A logical scalar. If \code{TRUE} (default), points are
#' added to the line plots to enhance visibility.
#' @param control.name A character string for the naming of the non-treated
#' units in placebo plots. Defaults to \code{"control units"}.
#' @return An object of class \code{\link[ggplot2]{ggplot}}.
#' @importFrom ggplot2 ggplot aes_string geom_line labs scale_x_date geom_hline
#' @importFrom ggplot2 scale_colour_manual scale_linetype_manual 
#' @importFrom ggplot2 scale_size_manual scale_alpha_manual geom_rect 
#' @importFrom ggplot2 geom_vline aes facet_wrap geom_smooth guides
#' @importFrom ggplot2 guide_legend stat_summary geom_point scale_y_continuous
#' @importFrom stats frequency
#' @importFrom utils stack
#' @method ggplot mscmt
#' @export 
ggplot.mscmt <- function(x,what,type=c("gaps","comparison","placebo.gaps",
                                       "placebo.data","p.value"),
                         treatment.time,zero.line=TRUE,ylab,xlab="Date",main,
                         col,lty,lwd,legend=TRUE,bw=FALSE,date.format,
                         unit.name,full.legend=TRUE,include.smooth=FALSE,
                         include.mean=FALSE,include.synth=FALSE,
                         draw.estwindow=TRUE,what.set,limits=NULL,
                         alpha=1,alpha.min=0.1,exclude.units=NULL,
                         exclude.ratio=Inf,ratio.type=c("rmspe","mspe"),
                         alternative=c("two.sided", "less", "greater"),
                         draw.points=TRUE,control.name="control units") {
  ratio.type  <- match.arg(ratio.type)                         
  alternative <- match.arg(alternative) 
  if (!missing(what.set)) 
    what.set <- match.arg(what.set,c("dependents","predictors","all")) 
  if (missing(type)&&(!is.null(x$placebo))) type <- "placebo.gaps"
  type        <- match.arg(type) 
  if (exclude.ratio<1) { 
    warning("exclude.ratio too small, using 1")
    exclude.ratio <- 1
  }  
  fixTimes <- function(a,b) {                                                   # helper function for results of compare()
    resA <- as.vector(a)
    names(resA) <- rep(rownames(a),times=ncol(a))
    resB <- as.vector(b)
    names(resB) <- rep(rownames(b),times=ncol(b))
    res <- rbind(Start=resA,End=resB)
    res[,apply(res,2,function(x) !all(is.na(x))),drop=FALSE]
  }
  if (is.null(x$placebo)) {                                                     # special treatment for placebo studies and comparions
    if (is.null(x$comparison)) tmp  <- x else {
      tmp <- x$comparison$results
      tmp$times.dep  <- fixTimes(tmp$dependent.start,tmp$dependent.end)
      tmp$times.pred <- fixTimes(tmp$predictor.start,tmp$predictor.end)
      tmp$dependent <- rownames(tmp$dependent)
    }
  } else tmp  <- x[[1]]
  if (missing(what)) what <- tmp$dependent[1]                                   # try to guess good defaults...
  if (!missing(what.set)) what <- switch(what.set,
    "dependents" = colnames(tmp$times.dep),
    "predictors" = colnames(tmp$times.pred),
    "all"        = c(colnames(tmp$times.dep),colnames(tmp$times.pred)))         # make this unique??? CHECK!!!
  if (!missing(treatment.time)) if (!inherits(treatment.time,"Date")) {         # convert treatment.time to Date
    if (is.numeric(treatment.time)) 
      treatment.time <- as.Date(paste0(treatment.time,"-01-01")) else
      treatment.time <- as.Date(treatment.time)
  }  
  if (!is.null(limits)) if (!inherits(limits,"Date")) {                         # convert limits to Date
    if (is.numeric(limits)) 
      limits <- as.Date(paste0(limits,"-01-01")) else
      limits <- as.Date(limits)
  }
  estwindow <- list()
  for (wh in what)
    estwindow[[wh]] <- cbind(
      apply(tmp$times.dep[,which(colnames(tmp$times.dep)==wh),drop=FALSE],2,
            AQM2Date),
      apply(tmp$times.pred[,which(colnames(tmp$times.pred)==wh),drop=FALSE],2,
            AQM2Date))
  estwindow <- lapply(estwindow,function(x) 
                        range(as.Date(as.vector(x),origin="1970-01-01")))
  res <- ggplot()                                                               # initialize the ggplot object
  if (!missing(treatment.time))                                                 # add line for treatment time
    res <- res + geom_vline(xintercept=as.numeric(treatment.time),
                            color="grey50")
  if ((type=="placebo.gaps")||(type=="placebo.data")||(type=="p.value")) {      # placebo-based plots
    if (is.null(x$placebo)) stop("results of placebo study missing")
    what.missing <- !(what %in% names(x$placebo))
    if (any(what.missing))
      stop(paste("variable(s)",paste0(what[what.missing],collapse=", "),
                 "missing in results"))
    unames <- colnames(x$placebo[[what[1]]]$gaps)
    mspe   <- sapply(x[-length(x)],function(x) x$loss.v)
    if (ratio.type=="rmspe") mspe <- sqrt(mspe)
    cunits <- mspe/mspe[1]<=exclude.ratio
    unames <- unames[cunits]
    cunits[exclude.units] <- FALSE
    unames <- setdiff(unames,exclude.units)
    ncontr <- length(unames)-1
    if (missing(col)&&bw) col=c("black","grey20")
    if (missing(col)&&(!bw)) col=c("red","grey20")
    if (missing(lty)&&bw) lty=c(1,2)
    if (missing(lty)&&(!bw)) lty=rep(1,2)
    if (missing(lwd)) lwd=c(2,1)
    if (missing(main)) main=if (length(what)==1) what[[1]] else ""
    if (missing(date.format)) 
      date.format <- if (frequency(x$placebo[[what[[1]]]]$gaps)>1) "%b %y" else 
                                                                   "%Y" 
  }  
  if (type=="p.value") {
    draw.estwindow <- FALSE
    dat  <- NULL
    for (wh in what) {
      gaps.post <- AQMtail(x$placebo[[wh]]$gaps[,cunits,drop=FALSE],
                           x[[1]]$times.dep[2,wh])
      pval  <- apply(gaps.post,1,switch(alternative,
                       two.sided = function(x) mean(abs(x[1])<=abs(x)),
                       less      = function(x) mean(x[1]>=x),
                       greater   = function(x) mean(x[1]<=x)))
      tmpdf <- ts2df(ts(pval,start=start(gaps.post),
                        frequency=frequency(gaps.post)))
      das   <- data.frame(date=tmpdf[,1],pvalue=as.numeric(tmpdf[,2]),
                          which.data=wh)
      dat   <- rbind(dat,das)
    }
    dat  <- dat[!is.na(dat$pvalue),]
    if (length(what)>1) res <- res + facet_wrap(~which.data,scales="free")
    if (missing(ylab)) ylab <- "p-value"
    res <- res + geom_point(data=dat,aes_string("date","pvalue")) +
                 scale_y_continuous(limits=c(0,1))
    if (zero.line) res <- res + geom_hline(yintercept=0,colour="grey50")
  }
  if ((type=="placebo.gaps")||(type=="placebo.data")) {
    if (missing(ylab)) ylab=paste(if (type=="placebo.gaps") "Gaps" else "Data",
                                  if (length(what)==1) 
                                    paste0(" for ",what[[1]],collapse=""))
    if (missing(unit.name)) unit.name <- "Unit"
    dat  <- NULL
    dat2 <- NULL
    for (wh in what) {
      tmpdf <- ts2df(if (type=="placebo.gaps") 
        x$placebo[[wh]]$gaps[,cunits,drop=FALSE] else 
        x$placebo[[wh]]$data.treat[,cunits,drop=FALSE])
      das <- cbind(date=tmpdf[,1],stack(tmpdf,select=-1))
      das <- cbind(das,treated=factor(
                                ifelse(das$ind==das$ind[1],"treated","control"),
                                levels=c("treated","control")),
                   start.estwindow=c(estwindow[[wh]][[1]],
                                   rep(NA,nrow(das)-!is.null(estwindow[[wh]]))),
                   end.estwindow=c(estwindow[[wh]][[2]],
                                   rep(NA,nrow(das)-!is.null(estwindow[[wh]]))),
                   which.data=wh)
      das$ind <- factor(das$ind,levels=c(as.character(das$ind[1]),
                                         setdiff(levels(das$ind),
                                                 as.character(das$ind[1]))))
      estwindow[[wh]] <- NULL                                                   
      dat  <- rbind(dat,das)
      dat2 <- rbind(dat2,cbind(ts2df(x$placebo[[wh]]$data.synth)[,1:2],
                               which.data=wh))
    }
    names(dat)[names(dat)=="ind"] <- unit.name
    names(dat2)[2] <- "value"
    dat  <- dat[!is.na(dat$values),]
    dat2 <- dat2[!is.na(dat2$value),]
    if (length(what)>1) res <- res + facet_wrap(~which.data,scales="free")
    if (alpha=="auto") {
      alpha <- rep(alpha.min,ncontr+1)
      names(alpha) <- levels(dat[[unit.name]])
      alpha[dat[[unit.name]][1]] <- 1
      alpha[names(x[[1]]$w)] <- alpha.min+(1-alpha.min)*x[[1]]$w
    } else if (length(alpha)==1) alpha <- rep(alpha,ncontr+1)
    if (full.legend) {
      lvals <- (1:6)[((seq_along(unames)-1)%%5)+1]                           
      res <- res + 
        geom_line(data=dat,aes_string("date","values",colour=unit.name,
                  linetype=unit.name,size=unit.name,alpha=unit.name),
                  na.rm=TRUE) +
        geom_line(data=dat[dat$treated=="treated",],aes_string("date","values",
                  colour=unit.name,linetype=unit.name,size=unit.name,
                  alpha=unit.name),na.rm=TRUE)
      if (draw.points) res <- res + 
        geom_point(data=dat,aes_string("date","values",
                   colour=unit.name,size=unit.name,alpha=unit.name),na.rm=TRUE) +
        geom_point(data=dat[dat$treated=="treated",],
                   aes_string("date","values",colour=unit.name,size=unit.name,
                   alpha=unit.name),na.rm=TRUE)
      res <- res + 
        scale_linetype_manual(values=lvals) +
        scale_size_manual(values=c(lwd[1],rep(lwd[2],length(unames)-1))) +
        scale_alpha_manual(values=alpha) +                          
        if (legend) guides(colour=guide_legend(override.aes=list(alpha=1))) else
                    guides(colour="none",linetype="none",size="none",
                           alpha="none")
    } else {
      res <- res + 
        geom_line(data=dat,aes_string("date","values",colour="treated",
                  linetype="treated",size="treated",group=unit.name,
                  alpha=unit.name),na.rm=TRUE) +
        geom_line(data=dat[dat$treated=="treated",],aes_string("date","values",
                  colour="treated",linetype="treated",size="treated",
                  alpha=unit.name,group=unit.name),na.rm=TRUE)
      if (draw.points) res <- res + 
        geom_point(data=dat,aes_string("date","values",
                   colour="treated",size="treated",group=unit.name,
                   alpha=unit.name),na.rm=TRUE) + 
        geom_point(data=dat[dat$treated=="treated",],
                   aes_string("date","values",colour="treated",size="treated",
                   alpha=unit.name,group=unit.name),na.rm=TRUE)
      res <- res + 
        scale_colour_manual("",values=col,
                            labels=c("treated unit",control.name)) +
        scale_linetype_manual("",values=lty,
                              labels=c("treated unit",control.name)) +
        scale_size_manual("",values=lwd,
                          labels=c("treated unit",control.name)) +
        scale_alpha_manual("",values=alpha) +                          
        if (legend) guides(colour=guide_legend(override.aes=list(alpha=1)),
                           alpha="none") else
                    guides(colour="none",linetype="none",size="none",
                           alpha="none")
    }  
    if (include.smooth) 
      res <- res + geom_smooth(data=dat[dat$treated=="control",],
                               mapping=aes_string("date","values"),
                               na.rm=TRUE) 
    if (include.mean) 
      res <- res + stat_summary(data=dat[dat$treated=="control",],
                                mapping=aes_string("date","values"),
                                fun.y="mean",
                                colour="black",geom="line",size=2,alpha=0.5,
                                na.rm=TRUE)
    if (include.synth&&(type=="placebo.data"))
      res <- res + geom_line(data=dat2,aes_string("date","value"),size=lwd[1],
                             linetype=if (full.legend) lvals[1] else lwd[1],
                             col="grey20",alpha=0.5)
    if (zero.line&&(type=="placebo.gaps")) 
      res <- res + geom_hline(yintercept=0,colour="grey50")
  }
  if ((type=="comparison")||(type=="gaps")) {                                   # non-placebo based plots
    if (is.null(x$combined)) {                                                  # input is a comparison
      if (missing(unit.name)) unit.name <- "Estimation"
      if (is.null(x$comparison)) stop("input is not an individual mscmt result")
      what.missing <- !(what %in% names(x$comparison$variables))
      if (any(what.missing))
        stop(paste("variable(s)",paste0(what[what.missing],collapse=", "),
                   "missing in results"))
      unames <- colnames(x$comparison$variables[[what[1]]]$gaps)
      nunits <- length(unames)
#      if (missing(col)&&bw) col=c("black","grey20")
#      if (missing(col)&&(!bw)) col=c("red","grey20")
      if (missing(col)) col=rep("black",nunits)
      if (missing(lty)) lty=rep(1,nunits)                                       # change this, see below?
      if (missing(lwd)) lwd=rep(1,nunits)
      if (missing(main)) main=if (length(what)==1) what[[1]] else ""
      if (missing(date.format)) date.format <- 
        if (frequency(x$comparison$variables[[what[[1]]]]$gaps)>1) "%b %y" else 
                                                                   "%Y" 
      if (missing(ylab)) ylab=paste(if (type=="gaps") "Gaps" else "Data",
                                    if (length(what)==1) 
                                      paste0(" for ",what[[1]],collapse=""))
      dat  <- NULL
      dat2 <- NULL
      for (wh in what) {
        tmpdf <- ts2df(if (type=="gaps") 
          x$comparison$variables[[wh]]$gaps else 
          x$comparison$variables[[wh]]$data.synth)
        das <- cbind(date=tmpdf[,1],stack(tmpdf,select=-1))
        das <- cbind(das,
                     start.estwindow=c(estwindow[[wh]][[1]],
                                   rep(NA,nrow(das)-!is.null(estwindow[[wh]]))),
                     end.estwindow=c(estwindow[[wh]][[2]],
                                   rep(NA,nrow(das)-!is.null(estwindow[[wh]]))),
                     which.data=wh)
#        das$ind <- factor(das$ind,levels=c(as.character(das$ind[1]),
#                                           setdiff(levels(das$ind),
#                                                   as.character(das$ind[1]))))
        estwindow[[wh]] <- NULL                                                   
        dat  <- rbind(dat,das)
        dat2 <- rbind(dat2,
                      cbind(ts2df(x$comparison$variables[[wh]]$data.treat)[,1:2],
                            which.data=wh))
      }
      names(dat)[names(dat)=="ind"] <- unit.name
      names(dat2)[2] <- "value"
      dat  <- dat[!is.na(dat$values),]
      dat2 <- dat2[!is.na(dat2$value),]
      if (length(what)>1) res <- res + facet_wrap(~which.data,scales="free")
      if (length(alpha)==1) alpha <- rep(alpha,nunits)
      lvals <- (1:6)[((seq_len(nunits)-1)%%5)+1]                                # change this, see above?
      res <- res + 
        geom_line(data=dat,aes_string("date","values",colour=unit.name,
                  linetype=unit.name,size=unit.name,alpha=unit.name),
                  na.rm=TRUE) +
        geom_line(data=dat[dat$treated=="treated",],aes_string("date","values",
                  colour=unit.name,linetype=unit.name,size=unit.name,
                  alpha=unit.name),na.rm=TRUE)
      if (draw.points) res <- res + 
        geom_point(data=dat,aes_string("date","values",colour=unit.name,
                   size=unit.name,shape=unit.name,alpha=unit.name),
                   na.rm=TRUE) +
        geom_point(data=dat[dat$treated=="treated",],
                   aes_string("date","values",colour=unit.name,size=unit.name,
                   shape=unit.name,alpha=unit.name),na.rm=TRUE) 
      res <- res + 
        scale_linetype_manual(values=lvals) +
        scale_size_manual(values=lwd) +
        scale_alpha_manual(values=alpha) +                          
        if (legend) guides(colour=guide_legend(override.aes=list(alpha=1))) else
                    guides(colour="none",linetype="none",size="none",
                           shape="none",alpha="none")    
        if ((type=="comparison")&&
            (length(unique(x$comparison$results$treated.unit))==1))
          res <- res + geom_line(data=dat2,aes_string("date","value"),size=2,
                                 linetype=1,col="black",alpha=0.5)
      if (zero.line&&(type=="gaps")) 
        res <- res + geom_hline(yintercept=0,colour="grey50")
                            
    } else {                                                                    # input is not a comparison
      what.missing <- !(what %in% names(x$combined))
      if (any(what.missing))
        stop(paste("variable(s)",paste0(what[what.missing],collapse=", "),
                   "missing in results"))
      dat <- NULL
      for (wh in what) {
        das <- ts2df(x$combined[[wh]])
        if (type=="comparison") das <- cbind(date=das[,1],stack(das,select=2:3))
        das <- cbind(das,
                    start.estwindow=c(estwindow[[wh]][[1]],
                                   rep(NA,nrow(das)-!is.null(estwindow[[wh]]))),
                    end.estwindow=c(estwindow[[wh]][[2]],
                                   rep(NA,nrow(das)-!is.null(estwindow[[wh]]))),
                    which.data=wh)
        estwindow[[wh]] <- NULL                                                   
        if (type=="comparison") 
          das$ind <- factor(das$ind,levels=c(as.character(das$ind[1]),
                                             setdiff(levels(das$ind),
                                                     as.character(das$ind[1]))))
        dat <- rbind(dat,das)
      }  
      dat  <- if (type=="comparison") dat[!is.na(dat$values),] else 
                                      dat[!is.na(dat$gaps),]
      if (length(what)>1) res <- res + facet_wrap(~which.data,scales="free")
      if (missing(date.format)) 
        date.format <- if (frequency(x$combined[[what[[1]]]])>1) "%b %y" else 
                                                                 "%Y" 
      if (type=="comparison") {
        if (missing(col)&&bw) col=rep("black",2)
        if (missing(col)&&(!bw)) col=c("black","red")
        if (missing(lty)&&bw) lty=c(1,2)
        if (missing(lty)&&(!bw)) lty=rep(1,2)
        if (missing(lwd)) lwd=c(2,2)
        if (missing(ylab)) ylab=if (length(what)==1) what else ""
        if (missing(main)) main=if (length(what)==1) 
                                  paste("Comparison of",what) else ""
        res <- res + 
          geom_line(data=dat,aes_string("date","values",colour="ind",size="ind",
                    linetype="ind"),na.rm=TRUE)
        if (draw.points) res <- res + 
          geom_point(data=dat,aes_string("date","values",
                     colour="ind",size="ind"),na.rm=TRUE)
        res <- res + 
          scale_colour_manual("",values=col,
                              labels=c("actual data","synthesized data")) +
          scale_linetype_manual("",values=lty,
                                labels=c("actual data","synthesized data")) +
          scale_size_manual("",values=lwd,
                            labels=c("actual data","synthesized data")) +
          if (legend) guides(colour=guide_legend(override.aes=list(alpha=1))) else
                      guides(colour="none",linetype="none",size="none")
      }
      if (type=="gaps") {
        if (missing(col)&&bw) col="black"
        if (missing(col)&&(!bw)) col="black"
        if (missing(lty)&&bw) lty=1
        if (missing(lty)&&(!bw)) lty=1
        if (missing(lwd)) lwd=2
        if (missing(ylab)) ylab="gap"
        if (missing(main)) main=if (length(what)==1) 
                                  paste("Gap for",what) else ""
        res <- res + geom_line(data=dat,
                               aes_string("date","gaps"),colour=col[1],
                               size=lwd[1],linetype=lty[1],na.rm=TRUE)
        if (draw.points) res <- res + geom_point(data=dat,
                                aes_string("date","gaps"),colour=col[1],
                                size=lwd[1],na.rm=TRUE) 
        if (zero.line) res <- res + geom_hline(yintercept=0,colour="grey50")
      }  
    }
  }  
  res <- res + labs(title=main,x=xlab,y=ylab) +
               scale_x_date(date_labels=date.format,limits=limits)
  if (draw.estwindow)
    res <- res + geom_rect(data=dat,
                           aes_string(xmin="start.estwindow",
                                      xmax="end.estwindow"),
                           ymin=-Inf,ymax=+Inf,alpha=0.3,fill="grey70",
                           na.rm=TRUE)
  res  
}
