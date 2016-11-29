## internal functions

M <- function(w,X,trafo.v) {
  M <- (diag(1,ncol(X))-rep(1,ncol(X))%*%t(w)) %*% t(X) %*% 
    (if (nrow(X)>1) diag(drop(X%*%w)) else drop(X%*%w))
  if (trafo.v$has.trafo) M %*% trafo.v$A else M 
}

all_v <- function(w,X,trafo.v,lb=1e-8) 
  cbind("min.loss.w" = loss_v(w,X,trafo.v,lb)$loss.w,single_v(w,X,trafo.v,lb))
                   
exists_v <- function(w,X,trafo.v,lb=1e-8,scale=196L,tol=1e-14)
  isTRUE(single_v(w,X,trafo.v,lb,check.exists=TRUE,tol_feasible=tol))

## currently not used
unique_v <- function(w,X,trafo.v,lb=1e-8,tol=1e-12) 
  loss_v(w,X,trafo.v,lb=lb,tol=tol)$unique.v
                   
## used in improveSynth
favourite_v <- function(w,X,trafo.v,lb=1e-8) {
  tmp <- loss_v(w,X,trafo.v,lb)$loss.w
  if (!any(is.na(tmp))) return(tmp)
  drop(single_v(w,X,trafo.v,lb))
}

# function delivers "loss.w"-v as well as whether v in general is unique
loss_v <- function(w,X,trafo.v,lb=1e-8,scale=196L,tol=1e-8) {                   #,scale=4L+8L+64L+256L
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
  lp_res <- lpr(direction="min",a0=lb*sum(obj),a=obj,d0=lb*K,d=rep(1,K),
                const.mat=const.mat,const.rhs=const.rhs,const.dir=const.dir,
                scale=scale)
  res <- lp_res[1:K] + lb
  comp <- res/sum(res)
  if (any(is.na(res))) unique.v <- NA else {
    k <- 1
    unique.v <- TRUE
    while (unique.v && (k<=K)) {
    	tmp <- lpr(direction="max",a0=lb,a=1*(k==1:K),d0=lb*K,d=rep(1,K),
                 const.mat=const.mat,const.rhs=const.rhs,const.dir=const.dir,
                 scale=scale)
      if (any(is.na(tmp))) 
        warning("lp for calculating unique.v could not be solved") else 
        unique.v <- (max(abs((tmp[1:K]+lb)/sum(tmp[1:K]+lb)-comp)) <= tol) 
      k <- k+1
    }
  }  	
  names(res) <- trafo.v$names.v
  list(loss.w=res,unique.v=unique.v)
}

#' @importFrom lpSolve lp
single_v <- function(w,X,trafo.v,lb=0,check.exists=FALSE,tol_feasible=1e-6,
                     tol_PUFAS=0.5*lb,tol_min=1e-6,scale=196L) {                #,scale=4L+8L+64L+256L

	M0   <- M(w,X,trafo.v)                                                        # matrix for calculating derivatives w.r.t. donor weights
	K    <- ncol(M0)                                                              # number of predictors
	J    <- nrow(M0)                                                              # number of donors
	
	tol_lb <- max(1e-14,lb*tol_feasible)
	
	# linear program written in terms of v_k
	const.mat <- const.rhs <- const.dir <- c()
	
	# conditions concerning derivatives w.r.t. w 
	const.mat <- rbind(const.mat,M0)
	const.rhs <- c(const.rhs,rep(0,J))
	const.dir <- c(const.dir,ifelse(w==0,">=","=="))

	# conditions that v_k <= 1 for all k
	const.mat <- rbind(const.mat,diag(K)) # 
	const.rhs <- c(const.rhs,rep(1,K))
	const.dir <- c(const.dir,rep("<=",K))
	
	# add new variable (minimum of v_1,...,v_K)
	const.mat <- cbind(const.mat,0)
	const.mat <- rbind(const.mat,cbind(diag(K),-1))
	const.rhs <- c(const.rhs,rep(0,K))
	const.dir <- c(const.dir,rep(">=",K))
	
	fixed <- c()                                                                  # components of v which are fixed
	res   <- numeric(K)
	names(res) <- trafo.v$names.v

  solve_fixed <- function(fixed,PUFAS=FALSE,return.min=FALSE) {                 # helper function
    n.fixed <- length(fixed)
    n.free  <- K - n.fixed
    if (n.fixed>0) {                                                            # remove currently fixed variables
	    eq_change <- J+fixed                                                      # change equations that currently fixed variables do not exceed upper bound 
	    const.rhs[eq_change] <- res[fixed]
	    const.dir[eq_change] <- "=="
      eq_remove <- J+K+fixed
	    const.mat <- const.mat[-c(eq_remove),,drop=FALSE]                         
      const.dir <- const.dir[-c(eq_remove)]                                     
      const.rhs <- const.rhs[-c(eq_remove)]
    }  
 	  lp_res <- lp(direction="max",objective.in=c(rep(0,K),1),                    # maximize minimum of free (non-fixed) variables
                 const.mat=const.mat,const.rhs=const.rhs,
                 const.dir=const.dir,scale=scale)
    if (return.min) 
      return(if (lp_res$status!=0) NA else lp_res$solution[K+1])                # return minimum if requested
 	  sol <- if (lp_res$status!=0) rep(NA,K) else lp_res$solution[1:K]            # Check whether lp was solved properly
 	  if (PUFAS) 
       list(solution = sol,
       	    pufas    = if (n.free < 3) TRUE else 
                         isTRUE(PUFAS(objective=c(rep(0,K),1),                  # change this to FALSE if sol contains NA ???
                                const.mat=const.mat,const.rhs=const.rhs,
                                const.dir=const.dir,solution=lp_res$solution,
                                scale=scale,tol=tol_PUFAS))
           ) else 
       sol
  }
	
	for (i in 1:max(1,(K-2))) {                                                   # loop through orders 
    if (check.exists) 
      return(solve_fixed(fixed,return.min=TRUE)>=(lb-tol_lb))
    not_fixed      <- setdiff(1:K,fixed)

    # checking for uniqueness via PUFAS
    lptmp          <- solve_fixed(fixed,PUFAS=TRUE)
    if (any(is.na(lptmp$solution))) break
    res[1:K]       <- lptmp$solution                                            # update free (non-fixed) components
    if (lptmp$pufas) break                                                      # PUFAS condition fulfilled (or only 2 free variables left)? => DONE!

    ind <- intersect(which(res<=min(res[not_fixed])*(1+tol_min)),not_fixed)     # determine components which - up to some tolerance - produce the current order statistics
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
    fixed <- c(fixed,best_ind)                                                  # update set of fixed indices
  }
  res <- res/max(res)
  if (!isTRUE(min(res)>=(lb-tol_lb))) res[1:K] <- NA                            # if lb <= v[i], return res (vector of NAs else)
  cbind("max.order"=res)
}
