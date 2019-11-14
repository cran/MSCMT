## internal functions

M <- function(w,X,trafo.v) {
  M <- (diag(1,ncol(X))-rep(1,ncol(X))%*%t(w)) %*% t(X) %*% 
    (if (nrow(X)>1) diag(drop(X%*%w)) else drop(X%*%w))
  if (trafo.v$has.trafo) M %*% trafo.v$A else M 
}

all_v <- function(w,X,Z,trafo.v,lb=1e-8,verbose=FALSE,debug=FALSE) 
  cbind("min.loss.w" = loss_v(w,X,Z,trafo.v,lb,verbose=verbose,debug=debug),
#        if (debug) "max.min"    = single_v(w,X,Z,trafo.v,lb,check.exists=TRUE,
#                                           verbose=verbose),
        "max.order"  = single_v(w,X,Z,trafo.v,lb,verbose=verbose,debug=debug))

exists_v <- function(w,X,Z,trafo.v,lb=1e-8,all=FALSE,...) 
  check_v(single_v(w,X,Z,trafo.v,lb,check.exists=TRUE,...),
          w,X,Z,trafo.v,lb,...) || 
  if (all) check_v(loss_v(w,X,Z,trafo.v,lb,...),
                   w,X,Z,trafo.v,lb,...) else FALSE

## used in improveSynth
favourite_v <- function(v,w,X,Z,trafo.v,lb=1e-8,...) {
  tmp <- loss_v(w,X,Z,trafo.v,lb,...)
  if (check_v(tmp,w,X,Z,trafo.v,lb,...)) tmp else {
    tmp <- single_v(w,X,Z,trafo.v,lb,check.exists=TRUE,...)
    if (check_v(tmp,w,X,Z,trafo.v,lb,...)) tmp else v
  }
}

check_v <- function(v,w,X,Z,trafo.v,lb=1e-8,tol.lp=1e-12,tol.loss=1e-8,
                    verbose=FALSE,debug=FALSE) {
  if (any(is.na(v))) {
    if (verbose||debug) catn("checking v: v contains NAs!") 
    invisible(FALSE)
  } else if ((min(v/max(v))<lb)||(max(v)>1)) {
    if (verbose||debug) {
      if (min(v/max(v))<lb) catn("smallest predictor weight too small with ",
                                 min(v)," (lb=",lb,")")
      if (max(v)>1)  catn("biggest predictor weight too big with ",max(v))
    }  
    invisible(FALSE)
  } else {
    M0        <- M(w,X,trafo.v)
    K         <- ncol(M0)
    J         <- nrow(M0)
    const.mat <- M0
    const.rhs <- rep(0,nrow(M0))
    const.dir <- ifelse(w==0,">=","==")
    lhs       <- const.mat%*%v
    eqmax     <- max(abs(lhs-const.rhs)[const.dir=="=="])
    ineqmax   <- max(c(0,(lhs-const.rhs)[const.dir=="<="],
                         (const.rhs-lhs)[const.dir==">="]))
    wnew      <- wnnlsExt(v,X,NULL,trafo.v,return.w=TRUE)
    maxwdev   <- max(pmax(abs(wnew-w),ifelse(w>0,abs((wnew-w)/w),0)))
    
    if (!is.null(Z)) {
      loss.old <- lossDep(Z,w)
      loss.new <- lossDep(Z,wnew)
      loss.d <- if (loss.old==0) loss.new else 
                  min((loss.new-loss.old)/loss.old,loss.new-loss.old)
    } else loss.d <- 0

    if (debug) 
      catn("max. eq. residuum: ",eqmax,", max. ineq. residuum: ",ineqmax,
           ", w deviation: ",maxwdev,", loss difference: ",loss.d)

    invisible((max(eqmax,ineqmax)<=tol.lp)&&(loss.d<=tol.loss))
  }
}

calc_lpq <- function(v,w,X,trafo.v) {
  if (any(is.na(v))) NA else {
    M0        <- M(w,X,trafo.v)
    K         <- ncol(M0)
    J         <- nrow(M0)
    const.mat <- M0
    const.rhs <- rep(0,nrow(M0))
    const.dir <- ifelse(w==0,">=","==")
    lhs       <- const.mat%*%v
    eqmax     <- max(abs(lhs-const.rhs)[const.dir=="=="])
    ineqmax   <- max(c(0,(lhs-const.rhs)[const.dir=="<="],
                         (const.rhs-lhs)[const.dir==">="]))
    max(eqmax,ineqmax)
  }
}

# function delivers "min.loss.w"-v
loss_v <- function(w,X,Z,trafo.v,lb=1e-8,verbose=FALSE,debug=FALSE) {
  M0   <- M(w,X,trafo.v)
  K    <- ncol(M0)

  # linear program written in terms of x_k := v_k - lb
  # conditions that x_k <= 1-lb for all k
  const.mat <- diag(K)
  const.rhs <- rep(1-lb,K)
  const.dir <- rep("<=",K)

  # conditions concerning derivatives w.r.t. w 
  const.mat <- rbind(const.mat,M0)
  const.rhs <- c(const.rhs,-lb*rowSums(M0))
  const.dir <- c(const.dir,ifelse(w==0,">=","=="))

  obj <- if (trafo.v$has.trafo) (X%*%w)[,1]^2 %*% trafo.v$A else (X%*%w)[,1]^2
  lp_res <- lprsolve(direction="min",a0=lb*sum(obj),a=obj,d0=lb*K,d=rep(1,K),   # solve lp via various solvers
                 const.mat=const.mat,const.rhs=const.rhs,const.dir=const.dir,
                 verbose=verbose,debug=debug)
  if (!any(is.na(lp_res))) res <- pmax(lb,pmin(1,lp_res[1:K] + lb)) else 
                           res <- lp_res
  names(res) <- trafo.v$names.v
  res <- sanitize_v(res,X,Z,trafo.v,lb,verbose=verbose,debug=debug)
  if (debug&&!check_v(res,w,X,Z,trafo.v,lb,verbose=TRUE,debug=TRUE))
    catn("min.loss.w v probably unreliable: ",calc_lpq(res,w,X,trafo.v))
  res
}

# function delivers "max.order"-v
single_v <- function(w,X,Z,trafo.v,lb=0,check.exists=FALSE,tol_PUFAS=0.5*lb,
                     tol_min=1e-6,tol_lp=1e-13,scale=196L,                      #,scale=4L+8L+64L+256L
                     cv.alpha=0,v.special=integer(),backup.v=NULL,verbose=FALSE,
                     debug=FALSE) {

  M0   <- M(w,X,trafo.v)                                                        # matrix for calculating derivatives w.r.t. donor weights
  K    <- ncol(M0)                                                              # number of predictors
  J    <- nrow(M0)                                                              # number of donors

  # linear program written in terms of v_k
  # conditions concerning derivatives w.r.t. w 
  const.mat <- M0
  const.rhs <- rep(0,J)
  const.dir <- ifelse(w==0,">=","==")
  if (sum(w>0)>1) {
    J <- J - 1
    const.mat <- const.mat[-which.max(w>0),,drop=FALSE]
    const.rhs <- const.rhs[-which.max(w>0)]
    const.dir <- const.dir[-which.max(w>0)]
  }
  
  # conditions that v_k <= 1 for all k
  const.mat <- rbind(const.mat,diag(K)) # 
  const.rhs <- c(const.rhs,rep(1,K))
  const.dir <- c(const.dir,rep("<=",K))

  if ((cv.alpha>0) && (length(v.special)>0)) {                                  # cross-validation?
    obj.special <- rep(0,K)
    obj.special[v.special] <- 1
    
    tmp <- lpsolve(direction="max",objective.in=obj.special,
                   const.mat=rbind(const.mat,diag(K)),
                   const.rhs=c(const.rhs,rep(lb,K)),
                   const.dir=c(const.dir,rep(">=",K)),scale=scale,tol.lp=tol_lp,
                   verbose=verbose,debug=debug)

    dtmp <- lp_diag(obj.special,const.mat=rbind(const.mat,diag(K)),
                    const.rhs=c(const.rhs,rep(lb,K)),
                    const.dir=c(const.dir,rep(">=",K)),tmp)

    if (max(dtmp[2:3])>tol_lp) {
      tmpv <- backup.v[,ncol(backup.v)]
      tmpv <- tmpv/max(tmpv)
      tmp <- sum(tmpv[v.special])
#      if (debug) {
        catn("maximum of special.v unreliable (inf:",max(dtmp[2:3]),"), using ",tmp,
             " determined by backup v.")
        print(dtmp)
#      }
    } else tmp <- tmp$objval    
    
    if (tmp*cv.alpha>lb*length(v.special)) {
      const.mat <- rbind(const.mat,obj.special)                                 # add constraint for cross-validation
      const.rhs <- c(const.rhs,cv.alpha*tmp)
      const.dir <- c(const.dir,">=")
      CV <- 1
    } else CV <- 0
  } else CV <- 0
  
  # add new variable (minimum of v_1,...,v_K)
  const.mat <- rbind(cbind(const.mat,0),cbind(diag(K),-1))
  const.rhs <- c(const.rhs,rep(0,K))
  const.dir <- c(const.dir,rep(">=",K))
  
  fixed      <- c()                                                             # components of v which are fixed
  not.fixed  <- 1:K
  res        <- numeric(K)
  names(res) <- trafo.v$names.v

  solve_fixed <- function(fixed=c(),PUFAS=FALSE) {                              # helper function
    not.fixed <- setdiff(1:K,fixed)
    n.fixed   <- length(fixed)
    n.free    <- K - n.fixed
    if (n.fixed>0) {                                                            # remove currently fixed variables
      eq_change <- J+fixed                                                      # change equations that currently fixed variables do not exceed upper bound 
      const.rhs[eq_change] <- res[fixed]
      const.dir[eq_change] <- "=="
      eq_remove <- J+K+CV+fixed
      const.mat <- const.mat[-c(eq_remove),,drop=FALSE]                         
      const.dir <- const.dir[-c(eq_remove)]                                     
      const.rhs <- const.rhs[-c(eq_remove)]
    }
    obj.in <- c(rep(0,K),1)
    lp_res <- lpsolve(obj.in,const.mat=const.mat,const.rhs=const.rhs,           # try to solve lp with various solvers
                      const.dir=const.dir,direction="max",tol.lp=tol_lp,
                      not.fixed=not.fixed,check.min=TRUE,verbose=verbose,
                      debug=debug)
                  
    sol <- if (lp_res$status!=0) rep(NA,K) else                                 # check whether lp was solved properly
                                 lp_res$solution[-length(lp_res$solution)]
    sol[(sol>0)&(sol<lb)] <- lb
    if (!PUFAS) return(sol) else return(
       list(solution = sol,
            pufas    = if (n.free < 3) TRUE else 
                         isTRUE(pufas(objective.in=c(rep(0,K),1),               # change this to FALSE if sol contains NA ???
                                const.mat=const.mat,const.rhs=const.rhs,
                                const.dir=const.dir,solution=lp_res$solution,
                                scale=scale,tol=tol_PUFAS,verbose=verbose,
                                debug=debug)))) 
  }
  
  maxmin <- solve_fixed()
  if (min(maxmin)>0) maxmin <- pmax(lb,pmin(1,maxmin))
  maxmin <- sanitize_v(maxmin,X,Z,trafo.v,lb,verbose=verbose,debug=debug)
  if (check.exists) return(maxmin)

  for (i in 1:max(1,(K-2))) {                                                   # loop through orders 

    # checking for uniqueness via PUFAS
    lptmp     <- solve_fixed(fixed,PUFAS=TRUE)
    if (any(is.na(lptmp$solution))||any(lptmp$solution<=0)) {
      if (debug) catn("max.order v possibly unreliable!")
      break
    }
    res[1:K] <- lptmp$solution                                                  # update free (non-fixed) components
    res <- sanitize_v(res,X,Z,trafo.v,lb,verbose=verbose,debug=debug)
    if (lptmp$pufas) break                                                      # PUFAS condition fulfilled (or only 2 free variables left)? => DONE!

    ind <- intersect(which(res<=min(res[not.fixed])*(1+tol_min)),not.fixed)     # determine components which - up to some tolerance - produce the current order statistics
    if (length(ind) == 1) {                                                     # if there is only one candidate, we use this one.
      best_ind <- ind
    }	else {                                                                    # if there is more than one candidate:
      best_value <- 0                                                           # initialize best_value
      best_ind   <- NA
      for (index in ind) {                                                      # loop through all candidates to find the one which maximizes the next order statistics
        cur_fixed <- c(fixed,index)
        res_tmp   <- solve_fixed(cur_fixed,PUFAS=FALSE)
        if (!any(is.na(res_tmp)))
          if (min(res_tmp[-c(cur_fixed)])>best_value) {                         # if currently found next order statistics is larger than current best
            best_value <- min(res_tmp[-c(cur_fixed)])                           # update best value ... 
            best_ind   <- index                                                 # ... and index
        }
      }
    }
    if (is.na(best_ind)) break
    fixed     <- c(fixed,best_ind)                                              # update set of fixed indices
    not.fixed <- setdiff(1:K,fixed)
  }
  res <- res/max(res)
  if (!isTRUE(min(res)>0)) res[1:K] <- NA else res <- pmax(lb,pmin(1,res))      # if lb <= v[i], return res (vector of NAs else)
  res <- sanitize_v(res,X,Z,trafo.v,lb,verbose=verbose,debug=debug)
  maxmin.ok   <- check_v(maxmin,w,X,Z,trafo.v,lb,verbose=debug,debug=debug)
  maxorder.ok <- check_v(res,w,X,Z,trafo.v,lb,verbose=debug,debug=debug)
  if (maxorder.ok) {
    if (maxmin.ok) bestifnear(cbind(res,maxmin),X,Z,trafo.v) else res
  } else {
    if (maxmin.ok) maxmin else res
  }
  if (debug&&!maxorder.ok)
    catn("max.order v probably unreliable: ",calc_lpq(res,w,X,trafo.v))
  res
}
