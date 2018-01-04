## internal functions and objects
## todo (?): @export multiOpt
## @export multiOpt
#' @importFrom parallel clusterApplyLB
multiOpt <- function(X0,X1=0,Z0,Z1=0,check.global=TRUE,
                     trafo.v,info.v,n.v,names.v,inner.optim="wnnlsOpt",
                     inner.opar=list(),starting.values=NULL,
                     outer.optim="DEoptC",outer.par=list(),outer.opar=list(),
                     std.v=c("sum","mean","min","max"),single.v=FALSE,
                     verbose=TRUE,debug=FALSE,seed=NULL,cl=NULL) {
  outer.optim   <- match.arg(outer.optim,c("DEoptC","DEoptim","GenSA","genoud",
                                           "nlminbOpt","JDEoptim","optimOpt",
                                           "DEopt","PSopt","nloptr",
                                           "hydroPSO","cma_es","malschains",
                                           "soma","psoptim","ga","crs","isres",
                                           "fixed","regression","none"),
                                           several.ok=TRUE)
  inner.optim   <- match.arg(inner.optim,c("wnnlsOpt","ipopOpt",
                                           "LowRankQPOpt","benchmarkOpt"))                                           
  do.optimize   <- !(all(outer.optim %in% c("fixed","regression","none")))
  checkPkg <- function(pkg) if(!requireNamespace(pkg,quietly=TRUE))
    stop("package ",pkg," not installed") else TRUE
  getFun <- function(pkg,fun) 
    if(!requireNamespace(pkg,quietly=TRUE))
      stop("package ",pkg," not installed") else 
      get(fun,envir=asNamespace(pkg),mode="function")
  
  solution.type <- "normal"
  std.v <- match.arg(std.v)
  X     <- X0 - drop(X1)                                                        # X0, X1 (Z0, Z1) have to be scaled already!
  Z     <- Z0 - drop(Z1)                                                            
  nW    <- ncol(X)
  cnX   <- colnames(X)
  storage.mode(X) <- "double"
  storage.mode(Z) <- "double"
  comb   <- NULL
  gRmspe <- NA
  X.orig <- X; Z.orig <- Z                                                      # make backup of original X&Z 
  
  if (missing(trafo.v)) {
    if (missing(info.v)) {
      if (missing(n.v))     n.v     <- nrow(X)
      if (missing(names.v)) names.v <- rownames(X)
      if (missing(trafo.v)) trafo.v <- genTrafo(n.v=n.v,names.v=names.v)
    } else {
      if (missing(n.v))     n.v     <- info.v$n.v 
      if (missing(names.v)) names.v <- info.v$names.v 
      if (missing(trafo.v)) trafo.v <- genTrafo(n.v=n.v,names.v=names.v,
                                                len.v=info.v$len.v)
    } 
  } else {
    if (missing(n.v))     n.v     <- trafo.v$n.v
    if (missing(names.v)) names.v <- trafo.v$names.v
  }  
  
  # recombine parameters for inner optimization
  inner.args <- list(debug=FALSE)
  inner.args[names(inner.opar)] <- inner.opar

  # recombine parameters for outer optimization
  outer.args1 <- list(lb=1e-8,opt.separate=TRUE)
  outer.args1[names(outer.par)] <- outer.par
  
  if ((outer.args1$lb<=0)||(outer.args1$lb>=1)) 
    stop("lb must be positive and must be smaller than 1")

  # initialize parameters for outer optimizer (depending on optimizer)
  if (any(outer.optim == "fixed")) {
    if (length(outer.optim)>1) 
      stop("optimizer 'fixed' not supported for multiple optimizations")
    solution.type <- "fixed"  
    if (is.null(outer.opar$v)) {
      warning("element v of outer.opar is missing, using constant weights.")
      outer.opar$v <- rep(1,n.v)
    }  
    v <- cbind("fixed"=outer.opar$v)
  }
  if (any(outer.optim == "regression")) {
    if (length(outer.optim)>1) 
      stop("optimizer 'regression' not supported for multiple optimizations")
    solution.type <- "regression"  
    if (is.null(outer.opar$add.const)) outer.opar$add.const <- TRUE
    if (is.null(outer.opar$solve.tol)) 
      outer.opar$solve.tol <- .Machine$double.eps
    Xr  <- cbind(if (outer.opar$add.const) 1,t(cbind(X1,X0)))
	  Zr  <- t(cbind(Z1,Z0))
	  v   <- cbind("regression" = if (outer.opar$add.const)
                   rowSums(solve(crossprod(Xr),crossprod(Xr,Zr),
                                 tol=outer.opar$solve.tol)[-1,]^2) else 
                   rowSums(solve(crossprod(Xr),crossprod(Xr,Zr),
                                 tol=outer.opar$solve.tol)^2))
  }
  if (any(outer.optim == "none")) {
    if (length(outer.optim)>1) 
      stop("optimizer 'none' not supported for multiple optimizations")
    solution.type <- "artificial"
    if (is.null(outer.opar$w)) 
      stop("outer.opar$w must be provided when outer.optim=='none'")
#    v     <- outer.opar$v                                                           # may be NULL
    v     <- cbind("fixed"=outer.opar$v)
    w     <- blow(outer.opar$w,cnX)
    rmspe <- sqrt(lossDep(Z,w))
    res   <- list(w=w,v=v,loss.v=rmspe^2,rmspe=rmspe,conv=0,single.v=single.v,
                  ncalls.inner=0)
  }

  # initialize variables for benchmarking if applicable
  if (inner.optim=="benchmarkOpt") {                                           
    globals$raise_wnnls <- numeric(2e6)
    globals$raise_ipop  <- numeric(2e6)
    globals$raise_LRQP  <- numeric(2e6)
    globals$vs          <- matrix(0,ncol=n.v,nrow=2e6)
  }

  if (n.v==1) {                                                                 # search space is null space?
    if (verbose) catn("Only one predictor, optimization not required.")
    res <- list(v=cbind("optimizer"=1),conv=0,single.v=TRUE,ncalls.inner=1)
  } else {

    #########################
    # Look for sunny donors #
    #########################
    is.sunny <- logical(nW)
    for (i in seq_len(nW)) is.sunny[i] <- isSunny(e(i,nW),X)
    nS <- sum(is.sunny)

    if (nS<=1&&(!any(outer.optim == "none"))) {                                   # at most one 'sunny' donor, special cases!
      if (nS==0) {                                                                # no 'sunny' donors
        if (verbose) catn("No 'sunny' donors!")
        solution.type <- "nosunny"
        res <- solveNoSunny(X,Z,trafo.v,verbose=verbose)			
      } else {                                                                    # one 'sunny' donor
        if (verbose) catn("Only one 'sunny' donor!")                              # v does not matter and does not need to be optimized
        solution.type <- "onesunny"
        res <- solveOneSunny(X[,is.sunny,drop=FALSE],Z[,is.sunny,drop=FALSE],
                             trafo.v)			
      }
      if ((!do.optimize)&&(!is.null(v))) {                                        # restore custom v if provided 
        v           <- cbind(v)
        colnames(v) <- outer.optim
        res$v       <- v
      }  
    } else {                                                                      # at least two 'sunny' donors, standard case
      if (verbose) catn("Number of 'sunny' donors: ",nS," out of ",nW)

      # restrict optimization to sunny donors
      X <- X[,is.sunny,drop=FALSE]; Z <- Z[,is.sunny,drop=FALSE] 

      ############################
      # Checks for global optima #
      ############################
      if (check.global) {
        res <- checkGlobalOpt(X.orig,Z.orig,trafo.v,outer.args1$lb,               # check for feasibility of 'true' outer optimum
                              single.v=single.v,verbose=verbose,debug=debug)
        gRmspe <- res$rmspe
        if (!is.na(res$conv)) {
          if (verbose) catn("Unrestricted outer optimum (obtained by ignoring all ",
                            "predictors) is FEASIBLE even when respecting the ",
                            "predictors.")
          if (do.optimize) solution.type <- "global"
          do.optimize   <- FALSE
        } else {
          if (verbose) catn("Unrestricted outer optimum (obtained by ignoring all ",
                            "predictors) with RMSPE ",gRmspe,
                            " and MSPE (loss v) ",gRmspe^2," is INFEASIBLE when ",
                            "respecting the predictors.")
        }
        if (do.optimize) {                                                        # check for feasibility of 'sunny' outer optimum 
          res <- checkGlobalOpt(X,Z,trafo.v,outer.args1$lb,single.v=single.v,
                                verbose=verbose,debug=debug)
          if (!is.na(res$conv)) {
            do.optimize <- FALSE
            v <- res$v                                                            # check this !!! XXX changed 2017-04-19
          }
        }  
      } 

      # check availability of packages for inner optimizer
      if ((inner.optim=="ipopOpt")||(inner.optim=="benchmarkOpt"))
        checkPkg("kernlab")
      if ((inner.optim=="LowRankQPOpt")||(inner.optim=="benchmarkOpt"))
        checkPkg("LowRankQP")
                   
      # pre-allocate workspace for inner optimizer wnnlsOpt
      if ((inner.optim=="wnnlsOpt")||(inner.optim=="benchmarkOpt")) {                                           
        ME  <- 1L
        MA  <- sum(trafo.v$len.v)
        N   <- ncol(X)
        N   <- as.integer(nS)
        MDW <- ME+MA
        globals$Ipar  <- as.integer(c(ME=ME,MA=MA,MDW=MDW,N=N)) 
        globals$IWORK <- integer(MDW+N)
        globals$WORK  <- double(MDW+5*N)
        globals$IWORK[1:2] <- c(length(globals$WORK),length(globals$IWORK))
        globals$RNORM <- double(1)
        globals$MODE  <- integer(1)
        globals$X     <- double(N)
      }
                 
      if (do.optimize&&(is.na(gRmspe)||is.na(res$conv))) {                      # outer optimization needed
        if (is.null(seed)) seed <- NA
        comb <- expand.grid(seed=seed,outer.optim=outer.optim,
                            stringsAsFactors=FALSE)
        if ((inner.optim=="benchmarkOpt")&&(nrow(comb)>1)) 
          stop("benchmarkOpt not supported for multiple optimizations")
        if ((length(outer.opar)>0)&&(!isTRUE(names(outer.opar)%in%outer.optim))) 
        {
          if (length(outer.optim)>1) 
            stop("names of outer.opar do not match outer optimizers")
          outer.opar=list(outer.opar)
          names(outer.opar) <- outer.optim
        }
        results <- cases <- vector("list",nrow(comb))
        times <- rep(NA,length(results))
        for (i in seq_along(cases)) {
          cases[[i]] <- list(outer.optim=comb$outer.optim[i],seed=comb$seed[i],
                             outer.opar=outer.opar[[comb$outer.optim[i]]])
        }
        #######################################
        # Do (multiple) outer optimization(s) #
        #######################################
        if (is.null(cl)) {                                                      # without cluster
          if (verbose&&(length(results)>1)) 
            catn("Starting ",length(results)," optimizations with ",
                 length(outer.optim)," different outer optimizers and ",
                 length(seed)," different seeds.")
          for (i in seq_along(results)) {
            if (verbose&&(length(results)>1)) 
              catn("Optimization ",i," out of ",length(results),":")
            results[[i]] <- 
              atomOpt(cases[[i]],X,Z,trafo.v,single.v,inner.optim,
                      inner.args,outer.args1,starting.values,verbose,
                      debug)
          }            
        } else {                                                                # with cluster
          if (verbose) 
            catn("Distributing ",length(results)," optimizations with ",
                 length(outer.optim)," different outer optimizers and ",
                 length(seed)," different seeds to the cluster. ",
                 "Please hold the line!")
  #          clusterEvalQ(cl,library(MSCMT))
  #          clusterExport(cl,c("atomOpt","cases","X","Z","trafo.v","single.v",
  #                             "inner.optim","inner.args","outer.args1",
  #                             "starting.values","verbose","debug"),
  #                        envir=environment())
          results <- clusterApplyLB(cl,cases,atomOpt,X,Z,trafo.v,single.v,
                                    inner.optim,inner.args,outer.args1,
                                    starting.values,verbose,debug)
          if (verbose) catn("Optimization on cluster finished!")
        }              
        times   <- sapply(results,function(x) x$user.self)                      # collect benchmark information
        rmspes  <- sapply(results,function(x) x$rmspe)
        loss.v  <- sapply(results,function(x) x$loss.v)
        n.inner <- sapply(results,function(x) x$ncalls.inner)
        all.w   <- sapply(results,function(x) x$w)
        all.v   <- sapply(results,function(x) x$v[,ncol(x$v)])
        dim(times) <- dim(rmspes) <- dim(loss.v) <- dim(n.inner) <- 
          c(length(seed),length(outer.optim))
        rownames(times) <- rownames(rmspes) <- rownames(loss.v) <- 
          rownames(n.inner) <- seed  
        colnames(times) <- colnames(rmspes) <- colnames(loss.v) <- 
          colnames(n.inner) <- outer.optim
        colnames(all.w) <- colnames(all.v) <- 
          paste0(rep(outer.optim,each=length(seed)),".",
                 rep(seed,times=length(outer.optim)))
        tmploss <- loss.v
        tmploss[tmploss==0] <- Inf
        bestrun <- which.min(tmploss)
        if (verbose&&(length(cases)>1)) {
          catn("Best optimization run(s):")
          for (i in bestrun) catn("Run ",i," with optimizer ",
            cases[[i]]$outer.optim," and seed ",cases[[i]]$seed,".")
          if (isTRUE(any(is.infinite(tmploss)))) {
          catn("Optimization run(s) with problems:")
          for (i in which(is.infinite(tmploss))) catn("Run ",i," with optimizer ",
            cases[[i]]$outer.optim," and seed ",cases[[i]]$seed,".")
          }
        }
        res     <- results[[bestrun[1]]]               
      } else {                                                                    # outer optimization not needed
        if (any(outer.optim == "none")) {
          if (is.null(v)) {
            w     <- blow(outer.opar$w,colnames(X))
            rmspe <- sqrt(lossDep(Z,w))
            if (exists_v(w,X,Z,trafo.v,outer.args1$lb)) {
              v <- if (single.v) single_v(w,X,Z,trafo.v,outer.args1$lb) else 
                                 all_v(w,X,Z,trafo.v,outer.args1$lb)
            } else v <- NA
          } else {
            w     <- blow(outer.opar$w,colnames(X.orig))
            rmspe <- sqrt(lossDep(Z.orig,w))
            v     <- cbind("fixed"=v)
          }
          res <- list(w=w,v=v,loss.v=rmspe^2,rmspe=rmspe,conv=0,
                      single.v=single.v,ncalls.inner=0)
        } else if (solution.type!="global") {
          fn.min.par  <- c(list(X=X,Z=Z,trafo=trafo.v$trafo),inner.args)   
          if (!is.matrix(v)) {                                                    # check this !!! XXX changed 2017-04-19
            v <- cbind(v)
            colnames(v) <- outer.optim
          }  
          v <- v[,apply(v,2,function(x) !any(is.na(x))),drop=FALSE]
          rownames(v) <- names.v
          mspe <- do.call(inner.optim,args=c(list(v=as.double(drop(v))),
                                             fn.min.par))
          w    <- do.call(inner.optim,args=c(list(v=as.double(drop(v))),
                                             fn.min.par,return.w=TRUE))
          names(w) <- colnames(X)
          res  <- list(w=w, v=v, loss.v=mspe, rmspe=sqrt(mspe), conv=0,
                       single.v=TRUE,ncalls.inner=2)
        }               
      } 
    }
  }
  
  res$v <- apply(res$v,2,function(x) x/do.call(std.v,list(x)))
  if (!is.matrix(res$v)) {                                                      # apply oversimplifies if length(v)==1
    res$v <- t(as.matrix(res$v))
    rownames(res$v) <- trafo.v$names.v
  }  

  ws    <- apply(res$v,2,function(x) wnnlsExt(x,X.orig,NULL,trafo.v,
                                              return.w=TRUE))
  ws    <- sanitize_w(ws)
  loss  <- apply(ws,2,function(w) lossDep(Z.orig,w))
  idx   <- which.min(loss)[1]
  res$loss.v   <- loss[idx]
  res$rmspe    <- sqrt(res$loss.v)
  rownames(ws) <- cnX
  res$w        <- ws[,idx]
  res$loss.w   <- apply(res$v,2,function(v) lossPred(X.orig,res$w,v,trafo.v))
  res$trafo.v     <- trafo.v
  res$outer.rmspe <- gRmspe
  
  if (verbose) {
    catf("Final rmspe: ",res$rmspe,", mspe (loss v): ",res$rmspe^2,"\n",                           
         "Optimal weights:")
    print(ess(res$w))
  }
  if (debug) {
    print(res$v)
    print(res$loss.w)
    print(rbind(loss=loss,
                increase=pmin((loss-min(loss))/min(loss),loss-min(loss))),
          digits=16)
  }
  if (verbose||debug) cat("\n")

  if (inner.optim=="benchmarkOpt") {                                           
    c(res,list(solution.type=solution.type,
      raise_wnnls=globals$raise_wnnls[1:globals$NRUNS],
      raise_ipop=globals$raise_ipop[1:globals$NRUNS],
      raise_LRQP=globals$raise_LRQP[1:globals$NRUNS],
      bench_vs=globals$vs[1:globals$NRUNS,],
      bench_X=X,bench_Z=Z))
  } else {
    if ((!is.null(comb))&&(nrow(comb)>1))
      c(res,list(solution.type=solution.type,
                 multi.opt=list(rmspe=rmspes,loss.v=loss.v,time=times,
                                n.inner=n.inner,w=all.w,v=all.v))) else
      c(res,list(solution.type=solution.type))
  }  
}

ggplotBench <- function(results,rmspe=TRUE,time.lim,loss.lim,xlab="",ylab="",
                        title="Benchmark for different optimizers",
                        what=c("Time","Loss")) {
  what <- match.arg(what,several.ok=TRUE)                        
  if (!inherits(results,"mscmt")) 
    stop("results is not an object of class mscmt")
  if(!requireNamespace("reshape",quietly=TRUE))
    stop("package reshape not installed")
  if (is.null(results$placebo)) {                      
    if (is.null(results$multi.opt))
      stop("mscmt results with multiple optimization needed")
    if (missing(time.lim)) time.lim <- range(results$multi.opt$time)
    if (missing(loss.lim)) loss.lim <- range(if (rmspe) results$multi.opt$rmspe 
                                                   else results$multi.opt$loss.v)
    Time <- cbind(reshape::melt(results$multi.opt$time),typ="Time",
                  min=min(time.lim),max=max(time.lim))
    Loss <- cbind(reshape::melt(if (rmspe) results$multi.opt$rmspe else 
                                           results$multi.opt$loss.v),
                  typ=if (rmspe) "RMSPE" else "MSPE",
                  min=min(loss.lim),max=max(loss.lim))
    All <- rbind(if ("Time" %in% what) Time, if ("Loss" %in% what) Loss)  
    All$X2 <- factor(All$X2,levels = rev(levels(All$X2)),ordered = TRUE)
    ggplot(All,aes_string("X2","value")) + ggplot2::geom_boxplot() + 
      ggplot2::coord_flip() + ggplot2::facet_grid(. ~ typ, scales="free") + 
      xlab(xlab) + ylab(ylab) + ggplot2::geom_blank(aes(y=min)) + 
      ggplot2::geom_blank(aes(y=max)) + ggplot2::ggtitle(title)
  } else {
    ind <- sapply(results[1:(length(results)-1)],
                  function(x) !is.null(x$multi.opt))
    if (!any(ind))
      stop("no mscmt results with multiple optimization in placebo results")
    if (!(missing(time.lim)&&missing(loss.lim)))
      warning("time.lim and loss.lim ignored for placebo results")  
    All <- NULL
    for (i in which(ind)) {  
      Time <- cbind(reshape::melt(results[[i]]$multi.opt$time),typ="Time",
                    treated=names(results)[i])
      Loss <- cbind(reshape::melt(if (rmspe) results[[i]]$multi.opt$rmspe else 
                                             results[[i]]$multi.opt$loss.v),
                    typ=if (rmspe) "RMSPE" else "MSPE",
                    treated=names(results)[i])
      All <- rbind(All,if ("Time" %in% what) Time, if ("Loss" %in% what) Loss)  
    }
    All$X2 <- factor(All$X2,levels = rev(levels(All$X2)),ordered = TRUE)
    ggplot(All,aes_string("X2","value")) + ggplot2::geom_boxplot() + 
      ggplot2::coord_flip() + ggplot2::facet_grid(treated~typ, scales="free") + 
      xlab(xlab) + ylab(ylab) + ggplot2::ggtitle(title)  
  }    
}
