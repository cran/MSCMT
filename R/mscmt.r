#' Multivariate SCM Using Time Series
#'
#' \code{mscmt} performs the Multivariate Synthetic Control Method Using Time 
#' Series.
#'
#' \code{mscmt} combines, if necessary, the preparation of the raw data (which 
#' is expected to be in "list" format, possibly after conversion from a 
#' \code{\link[base]{data.frame}}
#' with function \code{\link{listFromLong}}) and the call to the appropriate
#' MSCMT optimization procedures (depending on the input parameters).
#'
#' @param data Typically, a list of matrices with rows corresponding to times 
#' and columns corresponding to units for all relevant features (dependent as
#' well as predictor variables, identified by the list elements' names).
#' This might be the result of converting from a 
#' \code{\link[base]{data.frame}}
#' by using function \code{\link{listFromLong}}. 
#' 
#' For convenience, \code{data} may alternatively be the 
#' result of function \code{\link[Synth]{dataprep}} of package 
#' \code{'Synth'}. In this case, the parameters \code{treatment.identifier},
#' \code{controls.identifier}, \code{times.dep}, \code{times.pred}, 
#' and \code{agg.fns} are ignored, as these input parameters are generated
#' automatically from \code{data}. The parameters \code{univariate}, 
#' \code{alpha}, \code{beta}, and \code{gamma} are ignored by fixing them to 
#' their defaults.
#' Using results of \code{\link[Synth]{dataprep}} is experimental, because
#' the automatic generation of input parameters may fail due to lack of 
#' information contained in results of \code{\link[Synth]{dataprep}}.
#' 
#' @param treatment.identifier A character scalar containing the name of the 
#' treated unit. 
#' Must be contained in the column names of the matrices in \code{data}.
#' @param controls.identifier A character vector containing the names of at 
#' least two control units.
#' Entries must be contained in the column names of the matrices in \code{data}.
#' @param times.dep A matrix with two rows (containing start times in
#' the first and end times in the second row) and one column for each dependent
#' variable, where the column names must exactly match the names of the
#' corresponding dependent variables. 
#' A sequence of dates with the given start and end times of
#' \itemize{
#' \item annual dates, if the format of start/end time is "dddd", e.g. "2016",
#' \item quarterly dates, if the format of start/end time is "ddddQd", e.g. 
#' "2016Q1",
#' \item monthly dates, if the format of start/end time is "dddd?dd", e.g. 
#' "2016/03" or "2016-10",
#' }
#' will be constructed; these dates are looked for in the row names of
#' the respective matrices in \code{data}. 
#' @param times.pred A matrix with two rows (containing start times in
#' the first and end times in the second row) and one column for each predictor
#' variable, where the column names must exactly match the names of the
#' corresponding predictor variables.
#' A sequence of dates with the given start and end times of
#' \itemize{
#' \item annual dates, if the format of start/end time is "dddd", e.g. "2016",
#' \item quarterly dates, if the format of start/end time is "ddddQd", e.g. 
#' "2016Q1",
#' \item monthly dates, if the format of start/end time is "dddd?dd", e.g. 
#' "2016/03" or "2016-10",
#' }
#' will be constructed; these dates are looked for in the row names of
#' the respective matrices in \code{data}. 
#' @param agg.fns Either \code{NULL} (default) or a character vector containing
#' one name of an aggregation function for each predictor variable (i.e., each
#' column of \code{times.pred}). The character string "id" may be used as a
#' "no-op" aggregation. Each aggregation function must accept a numeric vector
#' and return either a numeric scalar ("classical" MSCM) or a numeric vector 
#' (leading to MSCM*T* if length of vector is at least two).
#' @param placebo A logical scalar. If \code{TRUE}, a placebo study is 
#' performed where, apart from the treated unit, each control unit is considered
#' as treated unit in separate optimizations. Defaults to \code{FALSE}. 
#' Depending on the number of control units and the complexity of the problem, 
#' placebo studies may take a long time to finish.
#' @param placebo.with.treated A logical scalar. If \code{TRUE}, the treated
#' unit is included as control unit (for other treated units in placebo 
#' studies). Defaults to \code{FALSE}.
#' @param univariate A logical scalar. If \code{TRUE}, a series of univariate
#' SCMT optimizations is done (instead of one MSCMT optimization) even if
#' there is more than one dependent variable. Defaults to \code{FALSE}.
#' @param univariate.with.dependent A logical scalar. If \code{TRUE} (and if
#' \code{univariate} is also \code{TRUE}), all dependent variables (contained
#' in the column names of \code{times.dep}) apart from the current (real) 
#' dependent variable are included as predictors in the series of univariate
#' SCMT optimizations. Defaults to \code{FALSE}.
#' @param check.global A logical scalar. If \code{TRUE} (default), a check for
#' the feasibility of the unrestricted outer optimum (where actually no 
#' restrictions are imposed by the predictor variables) is made before 
#' starting the actual optimization procedure.
#' @param inner.optim A character scalar containing the name of the optimization
#' method for the inner optimization. Defaults to \code{"wnnlsOpt"}, which
#' (currently) is the only supported implementation, because it outperforms
#' all other inner optimizers we are aware of. 
#' \code{"ipopOpt"}, which uses \code{\link[kernlab]{ipop}}, and 
#' \code{LowRankQPOpt}, which uses \code{\link[LowRankQP]{LowRankQP}} as inner
#' optimizer have experimental support for benchmark purposes.
#' @param inner.opar A list containing further parameters for the inner 
#' optimizer. Defaults to the empty list. (For \code{"wnnlsOpt"}, there are no
#' meaningful further parameters.)
#' @param outer.optim A character vector containing the name(s) of the 
#' optimization method(s) for the outer optimization. Defaults to \code{"DEoptC"}, 
#' which (currently) is the recommended global optimizer. 
#' The optimizers currently supported can be found in the documentation of
#' parameter \code{outer.opar}, where the default control parameters for
#' the various optimizers are listed.
#' If \code{outer.optim} has length greater
#' than 1, one optimization is invoked for each outer optimizer (and, 
#' potentially, each random seed, see below), and the best result is used.
#' @param outer.par A list containing further parameters for the outer 
#' optimization procedure. Defaults to the empty list. Entries in this list may 
#' override the following hard-coded general defaults:
#' \itemize{
#' \item \code{lb=1e-8}, corresponding to the lower bound for the ratio of
#' predictor weights,
#' \item \code{opt.separate=TRUE}, corresponding
#' to an improved outer optimization where each predictor is treated as the 
#' (potentially) most important predictor (i.e. with maximal weight) in 
#' separate optimizations (one for each predictor).
#' }
#' @param outer.opar A list (or a list of lists, if \code{outer.optim} has
#' length greater than 1) containing further parameters for the outer 
#' optimizer(s). Defaults to the empty list. Entries in this list may override
#' the following hard-coded defaults for the individual optimizers, which
#' are quite modest concerning the computing time. 
#' \code{dim} is a variable holding the problem dimension, 
#' typically the number of predictors minus one.
#' \tabular{lll}{
#' \bold{Optimizer}  \tab \bold{Package}     \tab \bold{Default parameters} \cr 
#' \code{DEoptC}     \tab \code{MSCMT}       \tab \code{nG=500}, \code{nP=20*dim}, \code{waitgen=100}, \cr
#'                   \tab                    \tab \code{minimpr=1e-14}, \code{F=0.5}, \code{CR=0.9} \cr
#' \code{cma_es}     \tab \code{cmaes}       \tab \code{maxit=2500} \cr
#' \code{crs}        \tab \code{nloptr}      \tab \code{maxeval=2.5e4}, \code{xtol_rel=1e-14}, \cr
#'                   \tab                    \tab \code{population=20*dim}, \code{algorithm="NLOPT_GN_CRS2_LM"} \cr
#' \code{DEopt}      \tab \code{NMOF}        \tab \code{nG=100}, \code{nP=20*dim} \cr
#' \code{DEoptim}    \tab \code{DEoptim}     \tab \code{nP=20*dim} \cr
#' \code{ga}         \tab \code{GA}          \tab \code{maxiter=50}, \code{monitor=FALSE}, \cr
#'                   \tab                    \tab \code{popSize=20*dim} \cr
#' \code{genoud}     \tab \code{rgenoud}     \tab \code{print.level=0}, \code{max.generations=70}, \cr
#'                   \tab                    \tab \code{solution.tolerance=1e-12}, \code{pop.size=20*dim}, \cr
#'                   \tab                    \tab \code{wait.generations=dim}, \code{boundary.enforcement=2}, \cr
#'                   \tab                    \tab \code{gradient.check=FALSE}, \code{MemoryMatrix=FALSE} \cr
#' \code{GenSA}      \tab \code{GenSA}       \tab \code{max.call=1e7}, \code{max.time=25/dim},  \cr
#'                   \tab                    \tab \code{trace.mat=FALSE} \cr
#' \code{hydroPSO}   \tab \code{hydroPSO}    \tab \code{maxit=300}, \code{reltol=1e-14}, \code{npart=3*dim} \cr
#' \code{isres}      \tab \code{nloptr}      \tab \code{maxeval=2e4}, \code{xtol_rel=1e-14}, \cr
#'                   \tab                    \tab \code{population=20*dim}, \code{algorithm="NLOPT_GN_ISRES"} \cr
#' \code{malschains} \tab \code{Rmalschains} \tab \code{popsize=20*dim}, \code{maxEvals=25000} \cr
#' \code{nlminbOpt}  \tab \code{MSCMT/stats} \tab \code{nrandom=30} \cr
#' \code{optimOpt}   \tab \code{MSCMT/stats} \tab \code{nrandom=25} \cr
#' \code{PSopt}      \tab \code{NMOF}        \tab \code{nG=100}, \code{nP=20*dim} \cr
#' \code{psoptim}    \tab \code{pso}         \tab \code{maxit=700} \cr
#' \code{soma}       \tab \code{soma}        \tab \code{nMigrations=100} 
#' }
#' If \code{outer.opar} is a list of lists, its names must correspond to (a 
#' subset of) the outer optimizers chosen in \code{outer.optim}.
#' @param std.v A character scalar containing one of the function names
#' "sum", "mean", "min", or "max" for the standardization of the predictor 
#' weights (weights are divided by \code{std.v(weights)} before reporting). 
#' Defaults to "sum", partial matching allowed.
#' @param alpha A numerical vector with weights for the dependent variables
#' in an MSCMT optimization or \code{NULL} (default). If not \code{NULL},
#' the length of \code{alpha} must agree with the number of dependent
#' variables, \code{NULL} is equivalent to weight 1 for all dependent 
#' variables.
#' @param beta Either \code{NULL} (default), a numerical vector, or a list.
#' If \code{beta} is a numerical vector or a list, its length must agree
#' with the number of dependent variables. 
#' \itemize{
#' \item If \code{beta} is a numerical vector,
#' the \code{i}th dependent variable is discounted with discount factor 
#' \code{beta[i]} (the observations of the dependent variables must thus be 
#' in chronological order!). 
#' \item If \code{beta} is a list, the components of \code{beta} must be 
#' numerical vectors with lengths corresponding to the numbers of observations 
#' for the individual dependent variables. These observations are then 
#' multiplied with the corresponding component of \code{beta}.
#' }
#' @param gamma Either \code{NULL} (default), a numerical vector, or a list.
#' If \code{gamma} is a numerical vector or a list, its length must agree
#' with the number of predictor variables. 
#' \itemize{
#' \item If \code{gamma} is a numerical vector,
#' the output of \code{agg.fns[i]} applied to the \code{i}th predictor variable
#' is discounted with discount factor \code{gamma[i]} (the output of 
#' \code{agg.fns[i]} must therefore be in chronological order!). 
#' \item If \code{gamma} is a list, the components of \code{gamma} must be 
#' numerical vectors with lengths corresponding to the lengths of the output of 
#' \code{agg.fns} for the individual predictor variables. The output of 
#' \code{agg.fns} is then multiplied with the corresponding component of 
#' \code{gamma}.
#' }
#' @param return.ts A logical scalar. If \code{TRUE} (default), most results are
#' converted to time series.
#' @param single.v A logical scalar. If \code{FALSE} (default), a selection
#' of feasible (optimal!) predictor weight vectors is generated. If \code{TRUE}, 
#' only one set of optimal predictor weights is generated.
#' @param verbose A logical scalar. If \code{TRUE} (default), output is verbose.
#' @param debug A logical scalar. If \code{TRUE}, output is very verbose. 
#' Defaults to \code{FALSE}.
#' @param seed A numerical vector or \code{NULL}. If not \code{NULL}, the
#' random number generator is initialized with the elements of \code{seed} via
#' \code{set.seed(seed)} (see \link[base]{Random}) before
#' calling the optimizer, performing repeated optimizations (and staying with 
#' the best) if \code{seed} has length greater than 1. Defaults to \code{NULL}. 
#' If not \code{NULL}, the seeds \code{int.seed} (default: 53058) and 
#' \code{unif.seed} (default: 812821) for \code{\link[rgenoud]{genoud}} are 
#' also initialized to the corresponding element of \code{seed}, but this can 
#' be overridden with the list elements \code{int.seed} and \code{unif.seed} 
#' of (the corresponding element of) \code{outer.opar}.
#' @param cl \code{NULL} (default) or an object of class \code{cluster}
#' obtained by \code{\link[parallel]{makeCluster}} of package \code{parallel}. 
#' Repeated estimations (see \code{outer.optim} and \code{seed}) and
#' placebo studies will make use of the cluster \code{cl} (if not \code{NULL}).
#' @return An object of class \code{"mscmt"}, which is essentially a list
#' containing the results of the estimation and, if applicable, the placebo
#' study.
#' The most important list elements are 
#' \itemize{
#' \item the weight vector \code{w} for the control units,
#' \item a matrix \code{v} with weight vectors for the predictors in its 
#' columns,
#' \item scalars \code{loss.v} and \code{rmspe} with the dependent loss and its 
#' square root,
#' \item a vector \code{loss.w} with the predictor losses corresponding to the
#' various weight vectors in the columns of \code{v},
#' \item a list of multivariate time series \code{combined} containing, 
#' for each dependent and predictor variable, a multivariate time series 
#' with elements \code{treated} for the actual values of the treated unit,
#' \code{synth} for the synthesized values and \code{gaps} for the differences.
#' }
#' Placebo studies produce a list containing individual results for each 
#' unit (as treated unit), starting with the original treated unit, as well
#' as a list element named \code{placebo} with aggregated results for each
#' dependent and predictor variable.
#' @importFrom stats sd
#' @importFrom parallel clusterExport clusterApplyLB clusterEvalQ
#' @rdname MSCMTfunction
#' @export mscmt
#' @examples
#' \dontrun{
#' ## for examples, see the package vignettes:
#' browseVignettes(package="MSCMT")
#' }
mscmt <- function(data,treatment.identifier=NULL, controls.identifier=NULL, 
                  times.dep=NULL, times.pred=NULL, agg.fns=NULL,
                  placebo=FALSE, placebo.with.treated=FALSE, univariate=FALSE,
                  univariate.with.dependent=FALSE,
                  check.global=TRUE, inner.optim="wnnlsOpt",inner.opar=list(),
                  outer.optim="DEoptC",outer.par=list(),
                  outer.opar=list(), std.v=c("sum","mean","min","max"),
                  alpha=NULL,beta=NULL,gamma=NULL,return.ts=TRUE,single.v=FALSE,
                  verbose=TRUE, debug=FALSE, seed=NULL, cl=NULL) {
  # check input                  
  std.v <- match.arg(std.v)
  
  # is data a dataprep object?
  is.dataprep <- isTRUE(all.equal(names(data),c("X0","X1","Z0","Z1","Y0plot",
                                  "Y1plot","names.and.numbers","tag")))
								  
  if (placebo&&(!is.null(cl))) { 
    n.placebo <- length(controls.identifier)+1
    n.multi   <- length(outer.optim) * (if (is.null(seed)) 1 else length(seed))
    placebo.on.cluster <- (n.placebo>n.multi)
  } else placebo.on.cluster <- FALSE
    
  if (is.dataprep) {
    storage.mode(data$X0) <- storage.mode(data$X1) <- 
      storage.mode(data$Z0) <- storage.mode(data$Z1) <- "double"
    univariate <- FALSE
    unit.names <- data$names.and.numbers$unit.names
    if (is.null(colnames(data$X1))) colnames(data$X1) <- colnames(data$Z1)      # correct potentially missing names
    if (is.null(colnames(data$Z1))) colnames(data$Z1) <- colnames(data$X1)
    names(unit.names) <- as.character(data$names.and.number$unit.numbers)
    treatment.identifier <- as.character(unit.names[colnames(data$X1)])
    controls.identifier  <- as.character(unit.names[colnames(data$X0)])
    if (any(is.na(treatment.identifier))) 
      treatment.identifier <- colnames(data$X1)
    if (any(is.na(controls.identifier))) 
      controls.identifier <- colnames(data$X0)
    times.dep  <- cbind("Y"=range(data$tag$time.optimize.ssr))
    times.pred <- matrix(rep(range(data$tag$time.predictors.prior),
                             times=length(data$tag$predictors)),nrow=2)
    for (i in seq_along(data$tag$special.predictors))
      times.pred <- cbind(times.pred,
                          range(data$tag$special.predictors[[i]][[2]]))
    colnames(times.pred) <- c(rownames(data$X1)[seq_along(data$tag$predictors)],
                              sapply(data$tag$special.predictors,
                                     function(x) x[[1]]))
    agg.fns <- c(rep(data$tag$predictors.op,length(data$tag$predictors)),
                 if (!is.null(data$tag$special.predictors)) 
                   sapply(data$tag$special.predictors,function(x) x[[3]])) 
  }

  # main workhorse
  synthSingle <- function(treatment.identifier, controls.identifier, 
                          times.dep, times.pred, agg.fns) {
    if (!is.matrix(times.dep)) stop("times.dep must be a matrix")
    dependent <- colnames(times.dep)
    predictor <- colnames(times.pred)
    time.optimize.ssr <- vector("list",length(dependent))
    for (i in seq_along(dependent))
      time.optimize.ssr[[i]] <- seqAQM(times.dep[1,i],times.dep[2,i])
  
    if (!is.matrix(times.pred)) stop("times.pred must be a matrix")
    if (is.null(agg.fns)) agg.fns <- rep("id",ncol(times.pred))
    if ((!is.character(agg.fns))||(length(agg.fns)!=ncol(times.pred)))
      stop("agg.fns is not a character vector with ncol(times.pred) elements")
    names.v    <- colnames(times.pred)
    len.v      <- length(names.v)
    special.predictors <- vector("list",len.v)
    for (i in 1:len.v) 
      special.predictors[[i]] <- 
        list(names.v[i],list(seqAQM(times.pred[1,i],times.pred[2,i])),
             if (is.null(agg.fns)) "id" else agg.fns[i])
  
    # prepare the data
    if (is.dataprep) {
      X  <- cbind(data$X0, data$X1)
      colnames(X) <- unit.names[as.character(colnames(X))]
      if (any(is.na(colnames(X)))) 
        colnames(X) <- c(colnames(data$X0),colnames(data$X1))
      X  <- X/apply(X, 1, sd)                                                   # scale X0,X1
      Z0 <- data$Z0
      colnames(Z0) <- unit.names[as.character(colnames(data$X0))]               # correct potentially wrong column names of Z0
      if (any(is.na(colnames(Z0)))) 
        colnames(Z0) <- colnames(data$X0)
      Z1 <- data$Z1
      colnames(Z1) <- unit.names[as.character(colnames(Z1))]
      if (any(is.na(colnames(Z1)))) 
        colnames(Z1) <- colnames(data$X1)
      Z  <- cbind(Z0,Z1)
      Y  <- rbind(cbind(data$Z0,data$Z1),cbind(data$Y0plot,data$Y1plot))
      Y  <- Y[sort(unique(row.names(Y))),]
      colnames(Y) <- colnames(X)
      Z1 <- Z[,treatment.identifier,drop=FALSE]
      Z0 <- Z[,controls.identifier,drop=FALSE]
      X1 <- X[,treatment.identifier,drop=FALSE]
      X0 <- X[,controls.identifier,drop=FALSE]
      dat <- list(X0 = X0, X1 = X1, Z0 = Z0, Z1 = Z1,
                  trafo.v=genTrafo(n.v=ncol(times.pred),
                                   names.v=colnames(times.pred)),
                  Z.scaled=FALSE)
    } else {
      dat <- prepare(data,predictors=NULL, predictors.op="mean", 
                     special.predictors, dependent, treatment.identifier, 
                     controls.identifier, time.predictors.prior=NULL, 
                     time.optimize.ssr, alpha=alpha, beta=beta, gamma=gamma, 
                     scale.Z=(ncol(times.dep)>1))
    }                 
  
    # run the optimization                           
    res <- multiOpt(X0=dat$X0,X1=dat$X1,Z0=dat$Z0,Z1=dat$Z1,trafo.v=dat$trafo.v,
                    check.global=check.global,
                    inner.optim=inner.optim,inner.opar=inner.opar,
                    starting.values=NULL,outer.par=outer.par,
                    outer.optim=outer.optim,outer.opar=outer.opar,
                    std.v=std.v,single.v=single.v,verbose=verbose,debug=debug,
                    seed=seed,cl=if (placebo.on.cluster) NULL else cl)
                    
    w <- blow(res$w,controls.identifier)
  
    if (is.dataprep) {
      ind        <- (w>0)
      data.synth <- list(Y=(Y[,names(w[ind]),drop=FALSE] %*% w[ind])[,1])
      data.treat <- Y[,treatment.identifier,drop=FALSE]
      colnames(data.treat) <- NULL
      data.treat <- list(Y=drop(data.treat))
    } else {
      ind        <- (w>0)
      data.synth <- lapply(data,function(x) 
                                  (x[,names(w[ind]),drop=FALSE] %*% w[ind])[,1])
      data.treat <- lapply(data,function(x) {
                             if (treatment.identifier %in% colnames(x))
                               tmp <- x[,treatment.identifier,drop=FALSE] else {
                                 tmp <- cbind(rep(NA,nrow(x)))
                                 rownames(tmp) <- rownames(x)
                               }  
                             colnames(tmp) <- NULL
                             drop(tmp)
                           })
    }                       
    gaps <- data.treat
    for (i in seq_along(gaps)) gaps[[i]] <- gaps[[i]] - data.synth[[i]]
    
    if (return.ts) {
      data.synth <- lapply(data.synth,AQM2ts)
      data.treat <- lapply(data.treat,AQM2ts)
      gaps       <- lapply(gaps,AQM2ts)
    }
    
    combined <- vector("list",length(gaps))
    names(combined) <- names(gaps)
    for (i in seq_along(gaps)) 
      combined[[i]] <- cbind(treated=data.treat[[i]],synth=data.synth[[i]],
                             gaps=gaps[[i]])
    names(agg.fns)  <- colnames(times.pred)                             
  
    res <- c(res,list(
               dataprep.scaled=dat,data.synth=data.synth,data.treat=data.treat,
               gaps=gaps,combined=combined,treated.unit=treatment.identifier,
               control.units=controls.identifier,dependent=dependent,
               predictor=predictor,agg.fns=agg.fns,agg.pred=rownames(dat$X0),
               times.dep=times.dep,times.pred=times.pred,std.v=std.v))          # insert alpha, beta, gamma, ...?
               
    if (dat$Z.scaled) {
      Zu    <- dat$Z0u - drop(dat$Z1u)
      mspeu <- lossDep(Zu,w)
      Zl    <- dat$Z.len
      Zw    <- as.numeric(Zu %*% w)^2
      Zsel  <- genA(Zl)
      res   <- c(res,list(mspes=as.numeric(t(Zsel)%*%Zw)/Zl),
                          mspe.unscaled=mspeu)
    }  
    res
  }           

  synthPlacebo <- function(treatment.identifier, controls.identifier, 
                           times.dep, times.pred, agg.fns) {
    all.units <- c(treatment.identifier,controls.identifier)
    mySynth <- function(treated) {
      if (verbose) catn("Using ",treated," as treated unit now.")
      res <- synthSingle(treated, 
        if (placebo.with.treated) setdiff(all.units,treated) else 
                                  setdiff(controls.identifier,treated),
        times.dep, times.pred, agg.fns)
      class(res) <- "mscmt"
      res
    }  
    if (verbose) catn("Starting placebo study, ",
      if (placebo.with.treated) "in" else "ex","cluding original treated ",
      "unit",if (placebo.on.cluster) ", on the cluster. Please hold the line",
      ".")
    if (placebo.on.cluster) clusterExport(cl,c("mySynth","treatment.identifier",
      "controls.identifier","times.dep","times.pred","agg.fns","all.units",
      setdiff(unique(agg.fns),"id"),"multiOpt","atomOpt"),envir=environment())
    res <- if (placebo.on.cluster) clusterApplyLB(cl,all.units,mySynth) else
                                   lapply(all.units,mySynth)
                            
    names(res) <- all.units
    combined <- vector("list",length(res[[1]]$combined))
    names(combined) <- names(res[[1]]$combined)
    for (j in names(combined)) {
      combined[[j]] <- list()
      for (i in seq_along(res)) {
        combined[[j]]$data.synth <- cbind(combined[[j]]$data.synth,
                                          res[[i]]$data.synth[[j]])
        combined[[j]]$data.treat <- cbind(combined[[j]]$data.treat,
                                          res[[i]]$data.treat[[j]])
        combined[[j]]$gaps       <- cbind(combined[[j]]$gaps,
                                          res[[i]]$gaps[[j]])
      }
      colnames(combined[[j]]$data.synth)<-colnames(combined[[j]]$data.treat)<-
        colnames(combined[[j]]$gaps) <- names(res)
    }    
    c(res,list(placebo=combined))
  }

  # separate univariate SCM for many dependent variables
  synthUnivariate <- function(treatment.identifier, controls.identifier, 
                              times.dep, times.pred, agg.fns) {
    if (verbose) catn("Starting univariate SCMT for single dependent ",
      "variables, ",if (univariate.with.dependent) "in" else "ex",
      "cluding other dependent variables as predictors.")
    synthFun <- if (placebo) synthPlacebo else synthSingle                               
    dependent <- colnames(times.dep)
    res <- vector("list",length(dependent))
    names(res) <- dependent
    for (i in seq_along(dependent)) {
      if (verbose) catn("Performing SCMT for ",names(res)[i]," now.")
      times.predT <- cbind(times.pred, if (univariate.with.dependent) 
                                         times.dep[,-i,drop=FALSE] else NULL)
      agg.fnsT <- if (is.null(agg.fns)) NULL else c(agg.fns,
        if (univariate.with.dependent) rep("id",length(dependent)-1) else NULL)

      res[[i]] <- synthFun(treatment.identifier, controls.identifier, 
                           times.dep[,i,drop=FALSE], times.predT, agg.fnsT)
      class(res[[i]]) <- "mscmt"
    }
    res
  }
  
  # what shall we do?
  synthFun <- if (univariate) synthUnivariate else 
                 if (placebo) synthPlacebo else synthSingle
  
  # do it!
  res <- synthFun(treatment.identifier, controls.identifier, times.dep, 
                  times.pred, agg.fns)
  class(res) <- "mscmt"

  if (placebo.on.cluster&&verbose) 
    catn("Placebo study on cluster finished.")
  res                  
}
