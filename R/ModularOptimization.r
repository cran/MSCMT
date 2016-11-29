## internal functions and objects
## but: export modOptim to be able to use X0,X1,Z0,Z1 representation?
## rather not, because improveSynth already has this capability...
## todo (?): @export modOptim

globals         <- new.env()

## @export modOptim
modOptim <- function(X0,X1=0,Z0,Z1=0,check.global=TRUE,
                     trafo.v,info.v,n.v,names.v,inner.optim="wnnlsOpt",
                     inner.opar=list(),starting.values=NULL,
                     outer.optim="DEoptC",outer.par=list(),outer.opar=list(),
                     std.v=c("sum","mean","min","max"),single.v=FALSE,
                     verbose=TRUE,debug=FALSE,seed=NULL) {
  outer.optim   <- match.arg(outer.optim,c("DEoptC","DEoptim","GenSA","genoud",
                                           "nlminbOpt","JDEoptim","optimOpt",
                                           "optimxOpt","DEopt","PSopt","nloptr",
                                           "hydroPSO","cma_es","malschains",
                                           "soma","psoptim","ga",
                                           "fixed","regression","none"))
  inner.optim   <- match.arg(inner.optim,"wnnlsOpt")                                           
  do.optimize   <- !(outer.optim %in% c("fixed","regression","none"))
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
  opt.separate <- outer.args1$opt.separate

  true.n.v <- if (opt.separate) n.v-1 else n.v
  
  # starting values for outer optimization?                 
  if (!is.null(starting.values)) 
    if (is.na(starting.values)) starting.values <- rep(1/(true.n.v),true.n.v)
  
  if ((outer.args1$lb<=0)||(outer.args1$lb>=1)) 
    stop("lb must be positive and must be smaller than 1")
  lb  <- log(outer.args1$lb,base=10)
  idx <- if (opt.separate) seq_len(n.v) else 1
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
                       max.generations=1000,
                       solution.tolerance=1e-12,pop.size=1000,
                       wait.generations=10,boundary.enforcement=2,
                       gradient.check=FALSE,MemoryMatrix=FALSE,
                       cluster=FALSE,transform=FALSE)
    outer.args <- c(outer.args,
                    list(Domains=cbind(rep(lb,true.n.v),rep(0,true.n.v)),
                         nvars=true.n.v,starting.values=starting.values))
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
    control <- list(maxit=1000)
    control[names(outer.opar)] <- outer.opar
    outer.args <- list(lower=rep(lb,true.n.v),upper=rep(0,true.n.v),
                       control=control)  
  } 
  if (outer.optim=="DEopt") {
    optfn <- getFun("NMOF","DEopt")
    algo <- list(min=rep(lb,true.n.v), max=rep(0,true.n.v), nG=500, 
                 minmaxConstr=TRUE, printBar=verbose, printDetail=verbose)
    algo[names(outer.opar)] <- outer.opar
    outer.args <- list(algo=algo)  
  }
  if (outer.optim=="PSopt") {
    optfn <- getFun("NMOF","PSopt")
    algo <- list(min=rep(lb,true.n.v), max=rep(0,true.n.v), nG=500, 
                 minmaxConstr=TRUE, printBar=verbose, printDetail=verbose)
    algo[names(outer.opar)] <- outer.opar
    outer.args <- list(algo=algo)  
  }
  if (outer.optim=="GenSA"){
    optfn <- getFun("GenSA","GenSA")
    control <- list(max.call=1e7)
    control[names(outer.opar)] <- outer.opar
    outer.args <- list(lower=rep(lb,true.n.v),upper=rep(0,true.n.v),
                       control=control)  
  }
  if (outer.optim=="ga") {
    optfn <- getFun("GA","ga")
    outer.args <- list(type="real-valued",min=rep(lb,true.n.v),
                       max=rep(0,true.n.v),maxiter=500,monitor=FALSE,
                       popSize=100)
    outer.args[names(outer.opar)] <- outer.opar
  }
  if (outer.optim=="soma") {
    optfn <- getFun("soma","soma")
    opts <- list(nMigrations=50)
    opts[names(outer.opar)] <- outer.opar
    outer.args <- list(bounds=list(min=rep(lb,true.n.v),max=rep(0,true.n.v)),
                       options=opts)
  }
  if (outer.optim=="cma_es"){
    optfn <- getFun("cmaes","cma_es")
    control <- list(maxit=1000)
    control[names(outer.opar)] <- outer.opar
    outer.args <- list(par=if(!is.null(starting.values)) starting.values else
                              log(rep(1/true.n.v,true.n.v),base=10),
                       lower=rep(lb,true.n.v),upper=rep(0,true.n.v),
                       control=control)  
  }
  if (outer.optim=="psoptim"){
    optfn <- getFun("pso","psoptim")
    control <- list(maxit=1000)
    control[names(outer.opar)] <- outer.opar
    outer.args <- list(par=if(!is.null(starting.values)) starting.values else
                              log(rep(1/true.n.v,true.n.v),base=10),
                       lower=rep(lb,true.n.v),upper=rep(0,true.n.v),
                       control=control)  
  }
  if (outer.optim=="malschains"){
    optfn <- getFun("Rmalschains","malschains")
    control <- list(popsize=100)
    if (!is.null(outer.opar$maxEvals)) {
      maxEvals <- outer.opar$maxEvals
      outer.opar$maxEvals <- NULL
    } else maxEvals <- 10000
    control[names(outer.opar)] <- outer.opar
    outer.args <- list(maxEvals=maxEvals, 
                       lower=rep(lb,true.n.v),upper=rep(0,true.n.v),
                       control=do.call(getFun("Rmalschains",
                                              "malschains.control"),control),
                       verbosity=0)  
  }
  if (outer.optim=="nloptr"){
    optfn <- getFun("nloptr","nloptr")
    opts <- list(maxeval=1e6,xtol_rel=1e-8,population=200,
                 ranseed=if(!is.null(seed)) as.integer(seed) else NULL,
                 algorithm="NLOPT_GN_CRS2_LM")
    opts[names(outer.opar)] <- outer.opar
    outer.args <- list(x0=if(!is.null(starting.values)) starting.values else
                            log(rep(1/true.n.v,true.n.v),base=10),
                       lb=rep(lb,true.n.v),ub=rep(0,true.n.v),
                       opts=opts)  
  }
  if (outer.optim=="optimxOpt") {
    checkPkg("optimx")
    optfn <- outer.optim
    outer.args <- list(lower=rep(lb,true.n.v),upper=rep(0,true.n.v))  
    outer.args[names(outer.opar)] <- outer.opar
  }
  if (outer.optim %in% c("nlminbOpt","optimOpt")) {
    checkPkg("stats")
    optfn <- outer.optim
    outer.args <- list(lower=rep(lb,true.n.v),upper=rep(0,true.n.v))  
    outer.args[names(outer.opar)] <- outer.opar
  }
  if (outer.optim == "fixed") {
    solution.type <- "fixed"  
    if (is.null(outer.opar$v)) {
      warning("element v of outer.opar is missing, using constant weights.")
      outer.opar$v <- rep(1,nrow(X))
    }  
    v <- cbind("fixed"=outer.opar$v)
  }
  if (outer.optim == "regression") {
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
  if (outer.optim == "none") {
    solution.type <- "artificial"
    if (is.null(outer.opar$w)) 
      stop("outer.opar$w must be provided when outer.optim=='none'")
  }

  # initialize best target function value
  best <- Inf

  # count calls to inner optimizer
  globals$NRUNS <- 0
  
  #########################
  # Look for sunny donors #
  #########################
  is.sunny <- logical(nW)
  for (i in seq_len(nW)) is.sunny[i] <- isSunny(e(i,nW),X)
  nS <- sum(is.sunny)
  if (nS>0) { 
    if (verbose) catn("Number of 'sunny' donors: ",nS," out of ",nW)
    
      if (nS==1) { # only one sunny donor?!
        # if there is only one sunny donor, 
        #   v does not matter and does not need to be optimized
        #   w is always given by 100% for the only sunny donor
        w <- 1
        names(w) <- cnX[is.sunny]
        res  <- list(w=w) # , v=opt.v, loss.v=rmspe^2, rmspe=rmspe, conv=conv)
        res$w       <- blow(res$w,cnX)
        return(c(res,list(solution.type="none",ncalls.inner=0)))
      }
    if (nrow(X)==1) 
    {
      # warning("Something's strange! Dimension of V equal to one, but more than one sunny donor! Returning NA.")
      # Hier koennte man auch aussen optimieren!
      w <- wnnlsGetGlobalOpt(Z[,is.sunny])$w
      return(list(w=blow(w,cnX)))
    }


  ############################
  # Check for global optimum #
  ############################
    res <- if (check.global&&do.optimize) 
             checkGlobalOpt(X,Z,trafo.v,outer.args1$lb,single.v=single.v) else 
             NULL 
    if (check.global&&do.optimize&&(!is.na(res$conv))) {
      if (verbose) catn("Global optimum is feasible.")
      solution.type <- "global"
    } else {
      if (verbose&&check.global&&do.optimize) 
        catn("Infeasible(!) global optimum has rmspe ",res$rmspe,
             " and mspe (loss v) ",res$rmspe^2,".")
             
             
      # make backup of original X&Z and restrict optimization to sunny donors
      X.orig <- X; Z.orig <- Z
      X <- X[,is.sunny,drop=FALSE]; Z <- Z[,is.sunny,drop=FALSE] 
      
      # pre-allocate workspace for inner optimizer wnnlsOpt
      if (inner.optim=="wnnlsOpt") {                                           
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
           
      if (do.optimize) {
  #########################
  # outer optimization    #
  #########################
        if (verbose) catn("Starting optimization via ",outer.optim,".") 
                                                
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
          res  <- list(w=w, v=opt.v, loss.v=rmspe^2, rmspe=rmspe, conv=conv)
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
                            else if (outer.optim=="nloptr") "eval_f" 
                            else if (outer.optim=="ga")     "fitness"
                            else if (outer.optim=="soma")   "costFunction"
                            else "fn"                      
            obj.fun.list <- if (outer.optim=="hydroPSO") list("obj.fun") else 
                                                         list(obj.fun)
            names(obj.fun.list) <- obj.fun.name                
    
            # this is the actual call of the outer optimizer:
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
            if (outer.optim=="nloptr") {
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
            if (outer.optim=="optimxOpt") {
              rmspe <- sqrt(rgV.optim$value)
              opt.v <- as.numeric(rgV.optim$par)
              conv  <- 0
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
              res  <- list(w=w, v=opt.v, loss.v=rmspe^2, rmspe=rmspe, conv=conv)
            }
          }
  	    }
    	  w <- res$w
    	  if (exists_v(w,X,trafo.v,outer.args1$lb)) {
          res$v <- if (single.v) 
                     single_v(w,X,trafo.v,outer.args1$lb) else
                     all_v(w,X,trafo.v,outer.args1$lb)
          res$single.v <- single.v
    	  } else { 
          warning("Optimal V could not be verified, using 'backup' V") 
          res$v           <- cbind("backup"=trafo.inner(opt.v))
          rownames(res$v) <- trafo.v$names.v
          res$single.v    <- NA
        } 
  
        if (verbose) catn("Optimization finished (with ",globals$NRUNS,
                          " calls to inner optimizer).")
      } else {
        if (outer.optim=="none") {
          if (is.null(outer.opar$v)) {
            w     <- blow(outer.opar$w,colnames(X))
            rmspe <- sqrt(lossDep(Z,w))
            if (exists_v(w,X,trafo.v,outer.args1$lb)) {
              v <- if (single.v) single_v(w,X,trafo.v,outer.args1$lb) else 
                                 all_v(w,X,trafo.v,outer.args1$lb)
            } else v <- NA
          } else {
            w     <- blow(outer.opar$w,colnames(X.orig))
            rmspe <- sqrt(lossDep(Z.orig,w))
            v     <- cbind("fixed"=outer.opar$v)
          }
          res <- list(w=w,v=v,loss.v=rmspe^2,rmspe=rmspe,conv=0,
                      single.v=single.v)
        } else {
          fn.min.par  <- c(list(X=X,Z=Z,trafo=trafo.v$trafo),inner.args)   
          rownames(v) <- rownames(X)
          mspe <- do.call(inner.optim,args=c(list(v=as.double(drop(v))),
                                             fn.min.par))
          w    <- do.call(inner.optim,args=c(list(v=as.double(drop(v))),
                                             fn.min.par,return.w=TRUE))
          names(w) <- colnames(X)
          v <- cbind(v)
          colnames(v) <- outer.optim
          res  <- list(w=w, v=v, loss.v=mspe, rmspe=sqrt(mspe), conv=0,
                       single.v=TRUE)
        }               
      } 
    }  
  } else {
    if (verbose) catn("No 'sunny' donors!")
    solution.type <- "nosunny"
    res <- solveNoSunny(X,Z,trafo.v,verbose=verbose)			
    if (!do.optimize) {
      v <- cbind(v)
      colnames(v) <- "outer.optim"
      res$v <- v
    }  
  } 

  res$v       <- apply(res$v,2,function(x) x/do.call(std.v,list(x)))
  if (!is.matrix(res$v)) {                                                      # apply oversimplifies if length(v)==1
    res$v <- t(as.matrix(res$v))
    rownames(res$v) <- trafo.v$names.v
  }  
  res$loss.w  <- apply(res$v,2,function(v) lossPred(X,res$w,v,trafo.v))
  res$trafo.v <- trafo.v
  res$w       <- blow(res$w,cnX)
  
  if (verbose) {
    catf("Final rmspe: ",res$rmspe,", mspe (loss v): ",res$rmspe^2,"\n",                           
         "Optimal weights:")
    print(ess(res$w))
#    if (!all(is.na(res$v))) {                                                   # v==NA if there are no sunny donors, no need to report then
#      catf(if (do.optimize) "Optimal ","v's",
#           if (single.v) " (single v requested)",
#           ":")
#      print(rbind(res$v,"----------"=NA,"pred. loss"=res$loss.w),na.print="")
#    }  
    cat0("\n")
  }	
  c(res,list(solution.type=solution.type,ncalls.inner=globals$NRUNS))
}

## inner optimization methods
## currently, only wnnlsOpt is supported
wnnlsOpt <- function(v,X,Z,trafo,debug=FALSE,return.w=FALSE,
                     check.ambiguity=FALSE,...) {
  if (any(is.na(v))) return(Inf) else tv  <- trafo(v)                           # needed for nlminb to work. Better drop nlminb-support for speed reasons?
  sol <-.Fortran(C_wnnls,W=.Call(C_prepareW4,X,tv),
                 MDW=globals$Ipar[3],ME=globals$Ipar[1],MA=globals$Ipar[2],
                 N=globals$Ipar[4],L=0L,PRGOPT=1.0,X=globals$X,
                 RNORM=globals$RNORM,MODE=globals$MODE,IWORK=globals$IWORK,
                 WORK=globals$WORK)
  globals$NRUNS <- globals$NRUNS + 1                 
  if ((any(is.infinite(sol$X))) || (sol$MODE>0)) warning("error in wnnls")                
  w   <- if (check.ambiguity) improveZw(Z,X,sol$X) else sol$X
  if (return.w) w else lossDep(Z,w)
}

## outer optimization with nlminb
nlminbOpt <- function(fn,lower,upper,nrandom=1,debug=FALSE,
                      frandom=function(n,lb,ub) stats::runif(n,lb,ub),
                      distance=1,iter.max=1000L,rel.tol=1e-10,trace=FALSE,...) {
  n       <- length(lower)
  best    <- Inf
  res     <- NULL
  control <- list(trace=trace,rel.tol=rel.tol,iter.max=iter.max)
  control[names(list(`...`))] <- list(`...`)
  if (nrandom==0) {
    if (n>16) stop("dimension of v too big for grid of starting values")        # this is quite arbitrary...
    mask <- bitwShiftL(1, 0:(n-1))
  }    
  nruns   <- if (nrandom>0) nrandom else 2^n
  for (i in 1:nruns) {
    s.v <- if (nrandom>0) frandom(n,lower,upper) else 
             ifelse(bitwAnd(i-1,mask),upper-distance,lower+distance)
    tmp <- stats::nlminb(s.v,fn,lower=lower,upper=upper,control=control)
    if (debug) cat("sv ",i," out of ",nruns," has objective ",tmp$objective, 
                   "(best: ",best,") after ",tmp$iterations," Iterations.\n")
    if (tmp$objective<best) {
      res  <- tmp
      best <- res$objective
    }
  }
  res
}

## outer optimization with optim
optimOpt <- function(fn,lower,upper,nrandom=1,debug=FALSE,
                     frandom=function(n,lb,ub) stats::runif(n,lb,ub),
                     distance=1,maxit=1000L,rel.tol=1e-10,trace=FALSE,...) {
  n       <- length(lower)
  best    <- Inf
  res     <- NULL
  control <- list(trace=trace,factr=rel.tol/.Machine$double.eps,maxit=maxit)
  control[names(list(`...`))] <- list(`...`)
  if (nrandom==0) {
    if (n>16) stop("dimension of v too big for grid of starting values")        # this is quite arbitrary...
    mask <- bitwShiftL(1, 0:(n-1))
  }    
  nruns   <- if (nrandom>0) nrandom else 2^n
  for (i in 1:nruns) {
    s.v <- if (nrandom>0) frandom(n,lower,upper) else 
             ifelse(bitwAnd(i-1,mask),upper-distance,lower+distance)
    tmp <- stats::optim(par=s.v,fn=fn,lower=lower,upper=upper,method="L-BFGS-B",
                        control=control)
    if (debug) cat("sv ",i," out of ",nruns," has objective ",tmp$value, 
                   "(best: ",best,") with conv. code ",tmp$convergence,".\n")
    if (tmp$value<best) {
      res  <- tmp
      best <- res$value
    }
  }
  res
}

## outer optimization with optimx
optimxOpt <- function(fn,lower,upper,nrandom=1,debug=FALSE,                     # this is broken ATM
                      frandom=function(n,lb,ub) stats::runif(n,lb,ub),
                      distance=1,itnmax=10000L,rel.tol=1e-10,trace=FALSE,...) {
  n       <- length(lower)
  best    <- Inf
  res     <- NULL
  control <- list(trace=trace,info=trace,
                  factr=rel.tol/.Machine$double.eps,reltol=rel.tol,
                  all.methods=TRUE)
  control[names(list(`...`))] <- list(`...`)
  if (nrandom==0) {
    if (n>16) stop("dimension of v too big for grid of starting values")        # this is quite arbitrary...
    mask <- bitwShiftL(1, 0:(n-1))
  }    
  nruns   <- if (nrandom>0) nrandom else 2^n
  for (i in 1:nruns) {
    s.v <- if (nrandom>0) frandom(n,lower,upper) else 
             ifelse(bitwAnd(i-1,mask),upper-distance,lower+distance)
    res <- optimx::optimx(par=s.v,fn=fn,lower=lower,upper=upper,                #method="L-BFGS-B",
                          itnmax=itnmax,control=control)
    res <- res[order(res$value,decreasing=FALSE),]
    val <- which(colnames(res)=="value") 
    tmp <- list(out.list=res,par=res[1,1:(val-1)],value=res[1,val])
    if (debug) print(tmp)
    if (debug) cat("sv ",i," out of ",nruns," has objective ",tmp$value, 
                   "(best: ",best,") with conv. code ",tmp$convergence,".\n")
    if (tmp$value<best) {
      res  <- tmp
      best <- res$value
    }
  }
  res
}
