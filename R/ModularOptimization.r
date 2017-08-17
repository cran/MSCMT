## internal functions and objects
globals         <- new.env()

timingBench <- function(vs,X,Z,max.iter=1000L,debug=FALSE,sigf=5,margin=5e-04,
                        bound=10,...) {
  ME  <- 1L
  MA  <- ncol(vs)
  n <- N <- ncol(X)
  MDW <- ME+MA
  globals$Ipar  <- as.integer(c(ME=ME,MA=MA,MDW=MDW,N=N)) 
  globals$IWORK <- integer(MDW+N)
  globals$WORK  <- double(MDW+5*N)
  globals$IWORK[1:2] <- c(length(globals$WORK),length(globals$IWORK))
  globals$RNORM <- double(1)
  globals$MODE  <- integer(1)
  globals$X     <- double(N)  
  
  wnnls.time <- system.time(
    for (i in 1:nrow(vs)) sol <- .Fortran(C_wnnls,W=.Call(C_prepareW4,X,vs[i,]),
                 MDW=globals$Ipar[3],ME=globals$Ipar[1],MA=globals$Ipar[2],
                 N=globals$Ipar[4],L=0L,PRGOPT=1.0,X=globals$X,
                 RNORM=globals$RNORM,MODE=globals$MODE,IWORK=globals$IWORK,
                 WORK=globals$WORK))
  ipop.time <- system.time(
    for (i in 1:nrow(vs)) w.ipop <- try(as.numeric(kernlab::ipop(c=rep(0,n),
                           H=fastMpdVM(X,vs[i,]),A=matrix(rep(1,n),nrow=1),b=1,
                           l=rep(0,n),u=rep(1,n),r=0,sigf=sigf,margin=margin,
                           maxiter=max.iter,bound=bound,verb=debug,...)$primal),
                           silent=TRUE))
  LRQP.time <- system.time(
    for (i in 1:nrow(vs)) w <- LowRankQP::LowRankQP(Vmat=t(sqrt(vs[i,])*X), 
                           dvec=rep(0,n), Amat=matrix(rep(1,n),nrow=1), bvec=1, 
                           uvec=rep(1,n),niter=max.iter,verbose=debug,
                           ...)$alpha)
  list(wnnls=wnnls.time,ipop=ipop.time,LRQP=LRQP.time)
}

## inner optimization methods
## currently, only wnnlsOpt is 'officially' supported
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

wnnlsExt <- function(v,X,Z,trafo.v,debug=FALSE,return.w=FALSE,
                     check.ambiguity=FALSE,...) {
  if (any(is.na(v))) return(Inf) else tv  <- trafo.v$trafo(v)                   # needed for nlminb to work. Better drop nlminb-support for speed reasons?
  ME  <- as.integer(1L)
  MA  <- as.integer(sum(trafo.v$len.v))
  N   <- as.integer(ncol(X))
  MDW <- as.integer(ME+MA)
  IWORK <- integer(MDW+N)
  IWORK[1:2] <- as.integer(c(MDW+5*N,MDW+N))
  sol <-.Fortran(C_wnnls,W=.Call(C_prepareW4,X,tv),
                 MDW=MDW,ME=ME,MA=MA,N=N,L=0L,PRGOPT=1.0,X=double(N),
                 RNORM=double(1),MODE=integer(1),IWORK=IWORK,
                 WORK=double(MDW+5*N))
  if ((any(is.infinite(sol$X))) || (sol$MODE>0)) warning("error in wnnls")                
  w   <- if (check.ambiguity) improveZw(Z,X,sol$X) else sol$X
  if (return.w) w else lossDep(Z,w)
}

LowRankQPOpt <- function(v,X,Z,trafo,debug=FALSE,return.w=FALSE,max.iter=1000L,
                         ...) {
  if (any(is.na(v))) return(Inf) else tv  <- trafo(v)                           # needed for nlminb to work. Better drop nlminb-support for speed reasons?
  n <- ncol(X)
  w <- LowRankQP::LowRankQP(Vmat=t(sqrt(tv)*X), dvec=rep(0,n), 
                            Amat=matrix(rep(1,n),nrow=1), bvec=1, uvec=rep(1,n),
                            niter=max.iter,verbose=debug,...)$alpha
  if (return.w) w else lossDep(Z,w)
}

ipopOpt <- function(v,X,Z,trafo,debug=FALSE,return.w=FALSE,max.iter=1000L,
                    sigf=5,margin=5e-04,bound=10,...) {
  if (any(is.na(v))) return(Inf) else tv  <- trafo(v)                           # needed for nlminb to work. Better drop nlminb-support for speed reasons?
  n <- ncol(X)
  w <- try(as.numeric(kernlab::ipop(c=rep(0,n),H=fastMpdVM(X,tv),
                                    A=matrix(rep(1,n),nrow=1),b=1,l=rep(0,n),
                                    u=rep(1,n),r=0,sigf=sigf,margin=margin,
                                    maxiter=max.iter,bound=bound,verb=debug,
                                    ...)$primal),silent=TRUE)
  if (inherits(w,"try-error")) { 
    if (return.w) return(NA) else return(Inf)
  } else {
    if (return.w) w else lossDep(Z,w)
  }  
}

benchmarkOpt <- function(v,X,Z,trafo,debug=FALSE,return.w=FALSE,max.iter=1000L,
                         sigf=5,margin=5e-04,bound=10,...) {
  if (any(is.na(v))) return(Inf) else tv  <- trafo(v)                           # needed for nlminb to work. Better drop nlminb-support for speed reasons?
  globals$NRUNS <- globals$NRUNS + 1                 
  globals$vs[globals$NRUNS,] <- tv
  
  # wnnls
  sol <-.Fortran(C_wnnls,W=.Call(C_prepareW4,X,tv),
                 MDW=globals$Ipar[3],ME=globals$Ipar[1],MA=globals$Ipar[2],
                 N=globals$Ipar[4],L=0L,PRGOPT=1.0,X=globals$X,
                 RNORM=globals$RNORM,MODE=globals$MODE,IWORK=globals$IWORK,
                 WORK=globals$WORK)
  if ((any(is.infinite(sol$X))) || (sol$MODE>0)) warning("error in wnnls")                
  w.wnnls <- sol$X

  # LowRankQP
  n <- ncol(X)
  w.LRQP <- LowRankQP::LowRankQP(Vmat=t(sqrt(tv)*X), dvec=rep(0,n), 
                            Amat=matrix(rep(1,n),nrow=1), bvec=1, uvec=rep(1,n),
                            niter=max.iter,verbose=debug,...)$alpha

  # ipop                      
  w.ipop  <- try(as.numeric(kernlab::ipop(c=rep(0,n),H=fastMpdVM(X,tv),
                                    A=matrix(rep(1,n),nrow=1),b=1,l=rep(0,n),
                                    u=rep(1,n),r=0,sigf=sigf,margin=margin,
                                    maxiter=max.iter,bound=bound,verb=debug,
                                    ...)$primal),silent=TRUE)
  if (inherits(w.ipop,"try-error")) { 
    trouble <- TRUE; w.ipop <- w.wnnls 
  } else trouble <- FALSE
  
  w.wnnls <- abs(w.wnnls)/sum(abs(w.wnnls))
  w.LRQP  <- abs(w.LRQP)/sum(abs(w.LRQP))
  w.ipop  <- abs(w.ipop)/sum(abs(w.ipop))
  
  Ws  <- cbind("wnnls"=w.wnnls,"LRQP"=w.LRQP,"ipop"=w.ipop)
  TFs <- sqrt(diag(t(Ws)%*%t(Z)%*%Z%*%Ws)/nrow(Z))
  LWs <- diag(t(Ws)%*%t(X)%*%diag(tv)%*%X%*%Ws)
  mLW <- min(LWs)
  
  globals$raise_wnnls[globals$NRUNS] <- (LWs[1]-mLW)/mLW
  globals$raise_LRQP[globals$NRUNS]  <- (LWs[2]-mLW)/mLW
  globals$raise_ipop[globals$NRUNS]  <- if (trouble) NA else (LWs[3]-mLW)/mLW
  
  if (return.w) w.wnnls else lossDep(Z,w.wnnls)
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
