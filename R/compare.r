#' Compare MSCMT estimation results
#'
#' \code{compare} collects estimation results from \code{\link{mscmt}} for
#' comparison purposes.
#'
#' \code{compare} collects (potentially many) estimation results from 
#' \code{\link{mscmt}} in a special object of class \code{"mscmt"}, which 
#' includes a component \code{"comparison"} where the different estimation 
#' results are aggregated.
#' This aggregated information is used by the \code{\link{ggplot.mscmt}} and 
#' \code{\link{print.mscmt}} methods to present summaries of the different
#' results.
#'
#' @param ... Objects of class \code{"mscmt"} or (a) list(s) containing objects
#' of class \code{"mscmt"}.
#' @param auto.name.prefix A character string (default: "") internally used to 
#' facilitate automatic naming in nested lists of unnamed estimation results.
#' 
#' @return An object of class \code{"mscmt"}, which itself contains the 
#' individual estimation results as well as a component \code{"comparison"} 
#' with aggregated information.
#' @importFrom stats sd
#' @export compare
compare <- function(...,auto.name.prefix="") {
  arg      <- list(`...`)
  argnames <- names(arg)
  if (is.null(argnames)) argnames <- rep("",length(arg))
  res <- list()
  agg <- NULL
  scalars <- c("loss.v","rmspe","loss.w","solution.type","treated.unit","std.v",
               "model.type")  
  XwN <- function(nam) {
    res        <- rep("x",length(nam))
    names(res) <- nam
    res
  }         
  myjoin <- function(a,b) {
    newrn <- unique(c(rownames(a),rownames(b)))
    res   <- matrix(NA,nrow=length(newrn),ncol=ncol(a)+ncol(b))
    colnames(res) <- c(colnames(a),colnames(b))
    rownames(res) <- newrn
    res[rownames(a),colnames(a)] <- a
    res[rownames(b),colnames(b)] <- b
    res
  }
  add <- function(a,b) {
    if (is.null(a)) return(b)
    if (length(unique(c(a$names,b$names)))!=length(a$names)+length(b$names)) 
      stop("names are not unique")
    res <- a
    res$names <- c(a$names,b$names)                                             # merge names
    for (le in names(res$variables)) {                                          # merge variables
      for (i in 1:3) {
        res$variables[[le]][[i]] <- 
          cbind(res$variables[[le]][[i]],b$variables[[le]][[i]])
        if (length(res$names)>1) colnames(res$variables[[le]][[i]]) <- res$names
      }    
    }
    for (le in scalars)                                                         # merge results
      res$results[[le]] <- c(a$results[[le]],b$results[[le]])
    for (le in setdiff(names(res$results),scalars)) 
      res$results[[le]] <- myjoin(a$results[[le]],b$results[[le]])
    res
  }
  if (("comparison" %in% names(arg))&&(!inherits(arg,"mscmt")))
    stop("'comparison' must not be one of the results' names")
  if (("placebo" %in% names(arg))&&(!inherits(arg,"mscmt")))
    stop("'placebo' must not be one of the results' names")
  for (i in seq_along(arg)) {
    if (inherits(arg[[i]],"mscmt")&&(!is.null(arg[[i]]$comparison))) {          # current argument already is a comparison
      res    <- c(res,arg[[i]][setdiff(names(arg[[i]]),"comparison")])
      agg    <- add(agg,arg[[i]]$comparison)
    } else
    if (((!inherits(arg[[i]],"mscmt"))&&is.list(arg[[i]]))||                    # current argument is a list, but not an object of class "mscmt" -> 'recursion'
       (inherits(arg[[i]],"mscmt")&&(is.null(arg[[i]]$placebo)&&
                                     is.null(arg[[i]]$combined))))  {           # current argument is a result of univariate "mscmt" estimation  -> 'recursion'
      tmpres <- do.call("compare",c(arg[[i]],
                        list(auto.name.prefix=paste0(auto.name.prefix,
                                                     LETTERS[i]))))
      res    <- c(res,tmpres[setdiff(names(tmpres),"comparison")])
      agg    <- add(agg,tmpres$comparison)
    } else
    if (inherits(arg[[i]],"mscmt")&&(is.null(arg[[i]]$comparison))) {           # current is object of class "mscmt", but not a comparison
      if (!is.null(arg[[i]]$placebo)) arg[[i]] <- arg[[i]][[1]]                 # placebo study? -> take first component!
      if (argnames[i]=="") argnames[i] <- 
        paste0("Result",auto.name.prefix,LETTERS[i])
      tmpR <- list(arg[[i]])
      names(tmpR) <- argnames[i]
      res  <- c(res,tmpR)                                                       # add whole 'mscmt' object
      tmp1 <- list(w=arg[[i]]$w,                                                # prepare estimation results
                   v=arg[[i]]$v[,ncol(arg[[i]]$v)],
                   loss.v=arg[[i]]$loss.v,                                     
                   rmspe=arg[[i]]$rmspe,
                   loss.w=arg[[i]]$loss.w[[length(arg[[i]]$loss.w)]],
                   solution.type=arg[[i]]$solution.type,
                   treated.unit=arg[[i]]$treated.unit,
                   control.units=XwN(arg[[i]]$control.units),
                   dependent=XwN(arg[[i]]$dependent),
                   predictor=XwN(arg[[i]]$predictor),
                   dependent.start=arg[[i]]$times.dep[1,,drop=FALSE],
                   dependent.end=arg[[i]]$times.dep[2,,drop=FALSE],
                   predictor.start=arg[[i]]$times.pred[1,,drop=FALSE],
                   predictor.end=arg[[i]]$times.pred[2,,drop=FALSE],
                   agg.pred=XwN(arg[[i]]$agg.pred),
                   agg.fns=arg[[i]]$agg.fns,
                   std.v=arg[[i]]$std.v,
                   model.type=paste0(if(length(arg[[i]]$dependent)>1) "M","SCM",
                     if(arg[[i]]$dataprep.scaled$trafo.v$has.trafo) "T")
                   )
      for (le in setdiff(names(tmp1),scalars)) {
        tmp1[[le]] <- if (is.matrix(tmp1[[le]])) t(tmp1[[le]]) else 
                                                 cbind(tmp1[[le]])
        colnames(tmp1[[le]]) <- argnames[i]
      }
      for (le in scalars) names(tmp1[[le]]) <- argnames[i]
      tmp2 <- vector("list",length(arg[[i]]$combined))                          # prepare combined data
      names(tmp2) <- names(arg[[i]]$combined)
      for (j in seq_along(arg[[i]]$combined)) 
        tmp2[[j]] <- list(data.treat=arg[[i]]$combined[[j]][,1],
                          data.synth=arg[[i]]$combined[[j]][,2],
                          gaps=arg[[i]]$combined[[j]][,3])
      tmp3 <- argnames[i]
      agg  <- add(agg,list(names=tmp3,variables=tmp2,results=tmp1))
    } else                                                                      # no valid input
#    if (!is.list(arg[[i]])&&(!inherits(arg[[i]],"mscmt")))                      
      warning("skipping argument which is neither of class 'list' nor 'mscmt'")
  }  
  final <- c(res,list(comparison=agg))
  class(final) <- "mscmt"
  if (!is.null(final$comparison)) final else NULL                               # do we have a valid result?
}
