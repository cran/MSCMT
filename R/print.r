#' Printing Results of MSCMT
#'
#' \code{print.mscmt} prints results of \code{mscmt}.
#'
#' A human-readable summary of \code{mscmt}'s results.
#'
#' @param x An object of class \code{"mscmt"}, usually obtained as
#' the result of a call to function \code{\link{mscmt}}.
#' @param ... Further arguments to be passed to or from other methods. 
#' They are ignored in this function.
#' @return Nothing useful (function is called for its side effects).
#' @importFrom utils capture.output
#' @method print mscmt
#' @export
print.mscmt <- function(x,...) {
  ind <- 16
  mwd <- options()$width

  printSingle <- function(x2) {
    has.placebo   <- !is.null(x2$placebo)
    if (has.placebo) x1 <- x2[[1]] else x1 <- x2
    Vmat <- capture.output(print(rbind(x1$v,"----------"=NA,
                                       "pred. loss"=x1$loss.w),
                                 na.print=""))
    paste0(
      "Specification:\n--------------\n\n",
      csl("Model type",
          paste0(if(length(x1$dependent)>1) "M","SCM",
                 if(x1$dataprep.scaled$trafo.v$has.trafo) "T"),
          mwd,ind),
      csl("Treated unit",x1$treated.unit,mwd,ind),
      csl("Control units",x1$control.units,mwd,ind),
      csl("Dependent(s)",align(x1$dependent," with optimization period from ",
                               x1$times.dep[1,]," to ",x1$times.dep[2,]),
          mwd,ind,TRUE),
      csl("Predictors",align(x1$predictor," from ",x1$times.pred[1,]," to ",
                             x1$times.pred[2,],
                             ifelse(x1$agg.fns=="id","",
                                    paste0(", aggregated via '",x1$agg.fns,"'"))
                             ),
          mwd,ind),"\n\n",
      "Results:\n--------\n\n",
      wrap("Result type",switch(x1$solution.type,
        nosunny = paste0("Perfect predictor fit possible, ie. predictor ",
          "weights V are meaningless since only W's with perfect preditor ",
          "fit are considered for optimization of the 'dependent' loss."),
        fixed = paste0("Vector V was fixed, outer optimization bypassed."),
        regression = paste0("Vector V was determined by regression, ",
          "outer optimization bypassed."),
        paste0("Ordinary solution, ie. no perfect preditor fit possible and ",
          "the predictors always impose some restrictions on the optimization ",
          "of the 'dependent' loss.")),mwd,ind),
      csl("Optimal W",align(names(ess(x1$w)),": ",
                            format(100*ess(x1$w),scientific=FALSE),"%"),
          mwd,ind),
      csl("Dependent loss",align(c("MSPE ('loss V')","RMSPE"),": ",
                                 c(x1$loss.v,x1$rmspe)),mwd,ind,TRUE),
      if(length(x1$dependent)>1) paste0(
        csl("Unscaled MSPE",x1$mspe.unscaled,mwd,ind),                               
        csl("Individ. MSPEs",align(x1$dependent,": ",x1$mspes),mwd,ind)                               
      ),  
      wrap("(Optimal) V",switch(x1$solution.type,
        nosunny = paste0("Predictor weights V are meaningless."),
        fixed = paste0("Fixed predictor weights V (not necessarily optimal)."),
        regression = paste0("Regression-based predictor weights V (not ",
                            "necessarily optimal)."),
        paste0(if (x1$single.v) "Single predictor weights V requested. ",
          if (!x1$single.v) "Some optimal " else "The optimal ",
          "weight vector",
          if (!x1$single.v) "s V are:" else " V is:")),mwd,ind),
      if (x1$solution.type!="nosunny") paste0(
        paste0(paste0(paste0(rep(" ",ind),collapse=""),Vmat),collapse="\n"),
        "\n",
        wrap("",paste0("(Predictor weights V are standardized by ",x1$std.v,
                       "(V)=1)"),mwd,ind)
      ),
      if (has.placebo) "\nResults of a placebo study are included.\n"
    )
  }

  is.univariate <- is.null(x$placebo) && is.null(x$combined)
  if (is.univariate) {
    out <- "Collection of univariate SCM(T) optimizations\n\n"
    for (i in seq_along(x)) out <- paste0(out,
      "Dependent variable: ",names(x)[i],"\n",
      paste0(rep("=",20+nchar(names(x)[i])),collapse=""),"\n\n",
      paste0(printSingle(x[[i]])),"\n"
    )
  } else out <- printSingle(x)
  cat(out,"\n")      
}
