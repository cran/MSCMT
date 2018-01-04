## internal functions and objects
## todo (?): @export atomOpt
## @export atomOpt
atomOpt <- function(arglist,X,Z,trafo.v,single.v,inner.optim,inner.args,
                    outer.args1,starting.values,verbose=TRUE,debug=FALSE) {
  timing <- system.time({                    
    checkPkg <- function(pkg) if(!requireNamespace(pkg,quietly=TRUE))
      stop("package ",pkg," not installed") else TRUE
    getFun <- function(pkg,fun) 
      if(!requireNamespace(pkg,quietly=TRUE))
        stop("package ",pkg," not installed") else 
        get(fun,envir=asNamespace(pkg),mode="function")
    
    # obtain parameters from arglist
    outer.optim  <- arglist$outer.optim
    outer.opar   <- arglist$outer.opar
    seed         <- arglist$seed
    if (!is.null(seed)) if (is.na(seed)) seed <- NULL
    
    # initialize further parameters
    opt.separate <- outer.args1$opt.separate
    n.v          <- trafo.v$n.v
    true.n.v     <- if (opt.separate) n.v-1 else n.v
    lb           <- log(outer.args1$lb,base=10)
    idx          <- if (opt.separate) seq_len(n.v) else 1
    
    # starting values for outer optimization?                 
    if (!is.null(starting.values)) 
      if (is.na(starting.values)) starting.values <- rep(1/(true.n.v),true.n.v)
    if (!is.null(starting.values)) starting.values <- log(starting.values,base=10)
    
    # initialize parameters for outer optimizer (depending on optimizer)
    if (outer.optim=="DEoptC") {
      optfn   <- "DEoptC"
      outer.args <- list(min=lb, max=0, nG=500, nP=20*true.n.v, waitgen=100, 
                         minimpr=1e-14, F=0.5, CR=0.9, 
                         check.ambiguity = FALSE,
                         width=if (verbose&&interactive()) getOption("width")-1 
                                 else 0)
      outer.args[names(outer.opar)] <- outer.opar
    }
    if (outer.optim=="genoud") {
      optfn <- getFun("rgenoud","genoud")
      outer.args <- list(print.level=0,
                         max.generations=70,
                         solution.tolerance=1e-12,pop.size=20*true.n.v,
                         wait.generations=true.n.v,boundary.enforcement=2,
                         gradient.check=FALSE,MemoryMatrix=FALSE,
                         cluster=FALSE,transform=FALSE)
      outer.args <- c(outer.args,
                      list(Domains=cbind(rep(lb,true.n.v),rep(0,true.n.v)),
                           nvars=true.n.v,starting.values=starting.values))
      if (!is.null(seed)) 
        outer.args <- c(outer.args,list(unif.seed=seed,int.seed=seed))
      outer.args[names(outer.opar)] <- outer.opar
    }
    if (outer.optim=="DEoptim") {
      optfn <- getFun("DEoptim","DEoptim")
      control <- list(trace=FALSE,reltol=1e-14,steptol=500)
      control[names(outer.opar)] <- outer.opar
      outer.args <- list(lower=rep(lb,true.n.v),upper=rep(0,true.n.v),
                         control=control)  
    } 
    if (outer.optim=="JDEoptim") {
      optfn <- getFun("DEoptimR","JDEoptim")
      outer.args <- list(lower=rep(lb,true.n.v),upper=rep(0,true.n.v))  
      outer.args[names(outer.opar)] <- outer.opar
    } 
    if (outer.optim=="hydroPSO") {
      optfn <- getFun("hydroPSO","hydroPSO")
      control <- list(maxit=300,write2disk=debug,verbose=verbose,reltol=1e-14,
	                  npart=3*true.n.v)
      control[names(outer.opar)] <- outer.opar
      outer.args <- list(lower=rep(lb,true.n.v),upper=rep(0,true.n.v),
                         control=control)  
    } 
    if (outer.optim=="DEopt") {
      optfn <- getFun("NMOF","DEopt")
      algo <- list(min=rep(lb,true.n.v), max=rep(0,true.n.v), nG=100, 
	               nP=20*true.n.v,
                   minmaxConstr=TRUE, printBar=verbose, printDetail=verbose)
      algo[names(outer.opar)] <- outer.opar
      outer.args <- list(algo=algo)  
    }
    if (outer.optim=="PSopt") {
      optfn <- getFun("NMOF","PSopt")
      algo <- list(min=rep(lb,true.n.v), max=rep(0,true.n.v), nG=100, 
	               nP=20*true.n.v, 
                   minmaxConstr=TRUE, printBar=verbose, printDetail=verbose)
      algo[names(outer.opar)] <- outer.opar
      outer.args <- list(algo=algo)  
    }
    if (outer.optim=="GenSA"){
      optfn <- getFun("GenSA","GenSA")
      control <- list(max.call=1e7,max.time=25/true.n.v,trace.mat=FALSE)
      control[names(outer.opar)] <- outer.opar
      outer.args <- list(lower=rep(lb,true.n.v),upper=rep(0,true.n.v),
                         control=control)  
    }
    if (outer.optim=="ga") {
      optfn <- getFun("GA","ga")
      outer.args <- list(type="real-valued",min=rep(lb,true.n.v),
                         max=rep(0,true.n.v),maxiter=50,monitor=FALSE,
                         popSize=20*true.n.v)
      outer.args[names(outer.opar)] <- outer.opar
    }
    if (outer.optim=="soma") {
      optfn <- getFun("soma","soma")
      opts <- list(nMigrations=100)
      opts[names(outer.opar)] <- outer.opar
      outer.args <- list(bounds=list(min=rep(lb,true.n.v),max=rep(0,true.n.v)),
                         options=opts)
    }
    if (outer.optim=="cma_es"){
      optfn <- getFun("cmaes","cma_es")
      control <- list(maxit=2500)
      control[names(outer.opar)] <- outer.opar
      outer.args <- list(par=if(!is.null(starting.values)) starting.values else
                                log(rep(1/true.n.v,true.n.v),base=10),
                         lower=rep(lb,true.n.v),upper=rep(0,true.n.v),
                         control=control)  
    }
    if (outer.optim=="psoptim"){
      optfn <- getFun("pso","psoptim")
      control <- list(maxit=700)
      control[names(outer.opar)] <- outer.opar
      outer.args <- list(par=if(!is.null(starting.values)) starting.values else
                                log(rep(1/true.n.v,true.n.v),base=10),
                         lower=rep(lb,true.n.v),upper=rep(0,true.n.v),
                         control=control)  
    }
    if (outer.optim=="malschains"){
      optfn <- getFun("Rmalschains","malschains")
      control <- list(popsize=20*true.n.v)
      if (!is.null(outer.opar$maxEvals)) {
        maxEvals <- outer.opar$maxEvals
        outer.opar$maxEvals <- NULL
      } else maxEvals <- 25000
      control[names(outer.opar)] <- outer.opar
      outer.args <- list(maxEvals=maxEvals, 
                         lower=rep(lb,true.n.v),upper=rep(0,true.n.v),
                         control=do.call(getFun("Rmalschains",
                                                "malschains.control"),control),
                         verbosity=0)  
    }
    if ((outer.optim=="nloptr")||(outer.optim=="crs")){
      optfn <- getFun("nloptr","nloptr")
      opts <- list(maxeval=2.5e4,xtol_rel=1e-14,population=20*true.n.v,
                   ranseed=if(!is.null(seed)) as.integer(seed) else NULL,
                   algorithm="NLOPT_GN_CRS2_LM")
      opts[names(outer.opar)] <- outer.opar
      outer.args <- list(x0=if(!is.null(starting.values)) starting.values else
                              log(rep(1/true.n.v,true.n.v),base=10),
                         lb=rep(lb,true.n.v),ub=rep(0,true.n.v),
                         opts=opts)  
    }
    if (outer.optim=="isres"){
      optfn <- getFun("nloptr","nloptr")
      opts <- list(maxeval=2e4,xtol_rel=1e-14,population=20*true.n.v,
                   ranseed=if(!is.null(seed)) as.integer(seed) else NULL,
                   algorithm="NLOPT_GN_ISRES")
      opts[names(outer.opar)] <- outer.opar
      outer.args <- list(x0=if(!is.null(starting.values)) starting.values else
                              log(rep(1/true.n.v,true.n.v),base=10),
                         lb=rep(lb,true.n.v),ub=rep(0,true.n.v),
                         opts=opts)  
    }
    if (outer.optim %in% c("nlminbOpt","optimOpt")) {
      checkPkg("stats")
      optfn <- outer.optim
      outer.args <- list(lower=rep(lb,true.n.v),upper=rep(0,true.n.v),
	    nrandom=if (outer.optim=="nlminbOpt") 30 else 25)  
      outer.args[names(outer.opar)] <- outer.opar
    }
  
    # initialize best target function value
    best <- Inf
  
    # count calls to inner optimizer
    globals$NRUNS <- 0
    
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
      MDW <- ME+MA
      globals$Ipar  <- as.integer(c(ME=ME,MA=MA,MDW=MDW,N=N)) 
      globals$IWORK <- integer(MDW+N)
      globals$WORK  <- double(MDW+5*N)
      globals$IWORK[1:2] <- c(length(globals$WORK),length(globals$IWORK))
      globals$RNORM <- double(1)
      globals$MODE  <- integer(1)
      globals$X     <- double(N)
    }
             
    #########################
    # outer optimization    #
    #########################
    if (verbose) catn("Starting optimization via ",outer.optim,
                      if (!is.null(seed)) paste0(", random seed ",seed),".") 
                                            
    if (outer.optim=="DEoptC") {
      if (!is.null(seed)) set.seed(seed) 
      trafo.inner <- function(v) 10^v
      trafo.fn    <- function(v) trafo.v$trafo(10^v)
      fn.min.par <- c(list(X=X,Z=Z,trafo=trafo.fn),inner.args) 
      rgV.optim <- 
        .Call(C_DE,X,Z,as.integer(trafo.v$len),as.integer(outer.args$nP),
              as.integer(outer.args$nG),as.double(outer.args$F),
              as.double(outer.args$CR),as.double(outer.args$min),
              as.double(outer.args$max),as.double(outer.args$minimpr),
              as.integer(outer.args$waitgen),
              as.logical(opt.separate),
              as.logical(outer.args$check.ambiguity),
              as.integer(outer.args$width))
      rmspe <- sqrt(rgV.optim$value)
      opt.v <- as.numeric(rgV.optim$par)
      conv  <- rgV.optim$counts
      globals$NRUNS <- rgV.optim$nruns
      w <- do.call(inner.optim,args=c(list(v=opt.v),fn.min.par,
                   return.w=TRUE))
      names(w) <- colnames(X)
      res  <- list(w=w, v=trafo.inner(opt.v), loss.v=rmspe^2, rmspe=rmspe, 
                   conv=conv)
    } else {
      for (i in idx) {
        if (opt.separate) {
          trafo.inner <- function(v) 
            c(if (i>1) 10^v[1:(i-1)] else NULL,
              1,if (i<n.v) 10^v[i:(n.v-1)] else NULL)
          trafo.fn    <- function(v) 
            trafo.v$trafo(c(if (i>1) 10^v[1:(i-1)] else NULL,
                            1,if (i<n.v) 10^v[i:(n.v-1)] else NULL))
        } else {
          trafo.inner <- function(v) 10^v
          trafo.fn    <- function(v) trafo.v$trafo(10^v)
        }                    
        
        if (!is.null(seed)) set.seed(seed) 
        fn.min.par <- c(list(X=X,Z=Z,trafo=trafo.fn),inner.args) 
        obj.fun <- if (outer.optim=="ga") 
                     function(v) 
                       (-1)*do.call(inner.optim,args=c(list(v=v),fn.min.par))
                   else 
                     function(v) 
                       do.call(inner.optim,args=c(list(v=v),fn.min.par))
        obj.fun.name <- if (outer.optim %in% c("DEopt","PSopt")) "OF" 
                        else if (outer.optim %in% c("nloptr","isres","crs")) 
                                                        "eval_f" 
                        else if (outer.optim=="ga")     "fitness"
                        else if (outer.optim=="soma")   "costFunction"
                        else "fn"                      
        obj.fun.list <- if (outer.optim=="hydroPSO") list("obj.fun") else 
                                                     list(obj.fun)
        names(obj.fun.list) <- obj.fun.name                
  
        # this is the actual call of the outer optimizer
        rgV.optim  <- do.call(optfn,args=c(obj.fun.list,outer.args))
  
        if (outer.optim=="genoud") {
          rmspe <- sqrt(rgV.optim$value)
          opt.v <- as.numeric(rgV.optim$par)
          conv  <- rgV.optim$generations
        }
        if (outer.optim=="DEoptim") {
          rmspe <- sqrt(rgV.optim$optim$bestval)
          opt.v <- as.numeric(rgV.optim$optim$bestmem)
          conv  <- rgV.optim$optim$iter
        }
        if (outer.optim=="JDEoptim") {
          rmspe <- sqrt(rgV.optim$value)
          opt.v <- as.numeric(rgV.optim$par)
          conv  <- rgV.optim$iter
        }
        if (outer.optim %in% c("DEopt","PSopt")) {
          rmspe <- sqrt(rgV.optim$OFvalue)
          opt.v <- as.numeric(rgV.optim$xbest)
          conv  <- algo$nG
        }
        if (outer.optim=="GenSA") {
          rmspe <- sqrt(rgV.optim$value)
          opt.v <- as.numeric(rgV.optim$par)
          conv  <- rgV.optim$counts
        }
        if (outer.optim %in% c("nloptr","isres","crs")) {
          rmspe <- sqrt(rgV.optim$objective)
          opt.v <- as.numeric(rgV.optim$solution)
          conv  <- rgV.optim$iterations
        }
        if (outer.optim=="hydroPSO") {
          rmspe <- sqrt(rgV.optim$value)
          opt.v <- as.numeric(rgV.optim$par)
          conv  <- rgV.optim$counts[2]
        }
        if (outer.optim=="ga") {
          rmspe <- sqrt((-1)*rgV.optim@fitnessValue)
          opt.v <- as.numeric(rgV.optim@solution)
          conv  <- rgV.optim@iter
        }
        if (outer.optim=="soma") {
          rmspe <- sqrt(rgV.optim$cost[rgV.optim$leader])
          opt.v <- as.numeric(rgV.optim$population[,rgV.optim$leader])
          conv  <- rgV.optim$migrations
        }
        if (outer.optim=="cma_es") {
          rmspe <- sqrt(rgV.optim$value)
          opt.v <- as.numeric(rgV.optim$par)
          conv  <- rgV.optim$counts[1]
        }
        if (outer.optim=="malschains") {
          rmspe <- sqrt(rgV.optim$fitness)
          opt.v <- as.numeric(rgV.optim$sol)
          conv  <- rgV.optim$numEvalEA + rgV.optim$numEvalLS
        }
        if (outer.optim=="psoptim") {
          rmspe <- sqrt(rgV.optim$value)
          opt.v <- as.numeric(rgV.optim$par)
          conv  <- rgV.optim$counts[2]
        }
        if (outer.optim=="nlminbOpt") {
          rmspe <- sqrt(rgV.optim$objective)
          opt.v <- as.numeric(rgV.optim$par)
          conv  <- rgV.optim$iterations
        }
        if (outer.optim=="optimOpt") {
          rmspe <- sqrt(rgV.optim$value)
          opt.v <- as.numeric(rgV.optim$par)
          conv  <- rgV.optim$counts[1]
        }
        
        # call the inner optimizer again to obtain w
        w <- do.call(inner.optim,args=c(list(v=opt.v),fn.min.par,
                     return.w=TRUE))
        names(w) <- colnames(X)
  
        if (verbose&&opt.separate)
          catn("Fixing v[",i,"] (as maximum) yields rmspe ",rmspe,
               " and mspe (loss v) ",rmspe^2," after ",conv,
               " generations/iterations.")
        if (rmspe < best) {
          best <- rmspe
          res  <- list(w=w, v=trafo.inner(opt.v), loss.v=rmspe^2, rmspe=rmspe, 
                       conv=conv)
        }
      }
    }
    res <- sanitize_res(res,X,Z,trafo.v,outer.args1$lb,verbose=verbose,
                        debug=debug)
    w <- res$w
    v <- res$v
    if (!check_v(v,w,X,Z,trafo.v,lb=outer.args1$lb,verbose=verbose,
                 debug=debug)) {
      if(verbose) catn("optimizer's v violates optimality constraints!") else 
                  warning("optimizer's v violates optimality constraints")
    }
    if (!is.na(single.v)) {
      v_minlossw  <- loss_v(w,X,Z,trafo.v,outer.args1$lb,verbose=verbose,
                            debug=debug)
      minlossw_ok <- check_v(v_minlossw,w,X,Z,trafo.v,outer.args1$lb,
                             verbose=verbose,debug=debug)
      v_maxorder  <- single_v(w,X,Z,trafo.v,outer.args1$lb,verbose=verbose,
                              debug=debug)
      maxorder_ok <- check_v(v_maxorder,w,X,Z,trafo.v,outer.args1$lb,
                             verbose=verbose,debug=debug)
      if (minlossw_ok) {
        v_minlossw <- bestifnear(cbind(v_minlossw,v),X,Z,trafo.v)
        if (maxorder_ok) 
          v_maxorder <- bestifnear(cbind(v_maxorder,v_minlossw),X,Z,trafo.v) else
          v_maxorder <- bestifnear(cbind(v,v_minlossw),X,Z,trafo.v)
      } else {
        v_minlossw <- v
        if (maxorder_ok) 
          v_maxorder <- bestifnear(cbind(v_maxorder,v),X,Z,trafo.v) else
          v_maxorder <- v
      }
      res$v <- if (single.v) cbind("max.order"=v_maxorder) else
                             cbind(if (debug) "optimizer"=v,
                                   "min.loss.w"=v_minlossw,
                                   "max.order"=v_maxorder)
    } else res$v <- cbind("optimizer"=res$v)
    res$single.v <- single.v
    rownames(res$v) <- trafo.v$names.v
  
    if (verbose) catn("Optimization finished (",globals$NRUNS," calls to ",
                      "inner optimizer), rmspe: ",res$rmspe,", mspe: ",
                      res$loss.v,".")
  })[["user.self"]]  
  c(res,list(ncalls.inner=globals$NRUNS,user.self=timing))
}
