## internal functions, maybe document and export some of the functions in 
## later release

################################################################################
################################################################################
# Helper functions for TIME SERIES HANDLING                                    #
################################################################################
################################################################################

## Guess frequency of time series data
##
## \code{freq} guesses the frequency of time series data.
##
## \code{freq} guesses the frequency of time series data.
##
## @param datum XX.
## @param only.first XX
## @return character string (one of annually, quarterly, monthly)
## TO DO(?): export?
## @export freq
freq <- function(datum,only.first=FALSE) {
	res <- ifelse(grepl("Q",datum),"quarterly",
			ifelse(nchar(datum)<=4,ifelse(nchar(datum)==4,"annually",NA),"monthly"))
	if (only.first) {
		if (length(unique(res))>1) 
      warning("frequency ambiguous, only first entry is used")
		res[1]  
	} else res	
}

## sequence-function for data of various granularity (annual, quarterly, 
## monthly) format for monthly data: YYYY-MM (alternatively, e.g., YYYY/MM) 
## format for quarterly data: YYYYQX (alternatively, e.g., YYYY-QX)

## Constructing names for time series data of various frequencies
##
## \code{seqAQM} constructs sequences of names for annual, qua^rterly and 
## monthly time series data.
##
## \code{seqAQM} constructs sequences of names for annual, quarterly and 
## monthly time series data.
##
## @param from XX.
## @param to XX
## @param by XX
## @return vector of character strings XX
## TO DO(?): export?
## @export seqAQM
seqAQM <- function(from,to,by=1) {
	fr   <- freq(c(from,to),only.first=TRUE)
	QtoI <- function(x) 4*as.numeric(substr(x,1,4))+as.numeric(substr(x,nchar(x),
                                                                    nchar(x)))-1
	ItoQ <- function(x) paste0(x%/%4,"Q",(x%%4)+1)
	MtoI <- function(x) 12*as.numeric(substr(x,1,4))+as.numeric(substr(x,6,7))-1
	ItoM <- function(x) paste0(x%/%12,"-",sprintf("%02d",(x%%12)+1))
	switch(fr,
		"quarterly" = ItoQ(seq(from=QtoI(from),to=QtoI(to),by=by)),
		"monthly"   = ItoM(seq(from=MtoI(from),to=MtoI(to),by=by)),
		as.character(seq(from=as.numeric(from),to=as.numeric(to),by=by))
	)
}

## convert data with character names (annual, quarterly, monthly data)
## to class 'ts'

## Convert (actual) time series data to class 'ts'
##
## \code{AQM2ts} converts annual, quarterly and monthly time series data to 
## class 'ts'.
##
## \code{AQM2ts} converts annual, quarterly and monthly time series data to 
## class 'ts'.
##
## @param data XX.
## @return object of class \code{ts}
#' @importFrom stats ts
#' @importFrom utils head tail
## TO DO(?): export?
## @export AQM2ts
AQM2ts <- function(data) {
  data  <- data[order(names(data))]
  sanitize <- function(nam,fre) 
    switch(as.character(fre),
           "4"  = paste0(substr(nam,1,4),"Q",substr(nam,nchar(nam),nchar(nam))),
           "12" = paste0(substr(nam,1,4),"-",substr(nam,6,7)),
           substr(nam,1,4))
  fnam  <- head(names(data),1)
  fr    <- switch(freq(fnam),"quarterly"=4,"monthly"=12,1)
  syear <- as.numeric(substr(fnam,1,4))
  smq   <- switch(as.character(fr),
                    "4"  = as.numeric(substr(fnam,nchar(fnam),nchar(fnam))),
                    "12" = as.numeric(substr(fnam,6,7)),
                  NULL)
  fnam  <- tail(names(data),1)
  eyear <- as.numeric(substr(fnam,1,4))
  emq   <- switch(as.character(fr),
                    "4"  = as.numeric(substr(fnam,nchar(fnam),nchar(fnam))),
                    "12" = as.numeric(substr(fnam,6,7)),
                  NULL)
  AQM   <- seqAQM(head(names(data),1),tail(names(data),1))
  res   <- rep(NA,length(AQM))
  names(res) <- AQM
  res[sanitize(names(data),fr)] <- as.numeric(data)
  ts(as.numeric(res),start=c(syear,smq),end=c(eyear,emq),frequency=fr)
}

## @export AQM2Date
AQM2Date <- function(rng,from,to) {
  if (missing(from)) from <- rng[1]
  if (missing(to))   to   <- rng[2]
  from <- switch(freq(from),
            "annually"  = as.Date(paste0(from,"-01-01")),
            "monthly"   = as.Date(paste0(from,"-01")),
            "quarterly" = as.Date(paste0(substr(from,1,4),"-",
                                    as.numeric(substr(from,nchar(from),
                                                      nchar(from)))*3-2,"-01")))
  to   <- switch(freq(to),
            "annually"  = as.Date(paste0(to,"-12-31")),
            "monthly"   = as.Date(paste0(to,"-28")),
            "quarterly" = as.Date(paste0(substr(to,1,4),"-",
                                    as.numeric(substr(to,nchar(to),
                                                      nchar(to)))*3,"-28")))
  c(from=from,to=to)
}

#' @importFrom stats window
## @export AQMwindow
AQMwindow <- function(data,range) {
  fr   <- switch(freq(range,only.first=TRUE),"quarterly"=4,"monthly"=12,1)
  year <- as.numeric(substr(range,1,4))
  mq   <- switch(as.character(fr),
                   "4"  = as.numeric(substr(range,nchar(range),nchar(range))),
                   "12" = as.numeric(substr(range,6,7)),
                 NULL)
  window(data,start=c(year[1],mq[1]),end=c(year[2],mq[2]))
}

#' @importFrom stats window
## @export AQMtail
AQMtail <- function(data,after) {
  fr   <- switch(freq(after,only.first=TRUE),"quarterly"=4,"monthly"=12,1)
  year <- as.numeric(substr(after,1,4))
  mq   <- switch(as.character(fr),
                   "4"  = as.numeric(substr(after,nchar(after),nchar(after))),
                   "12" = as.numeric(substr(after,6,7)),
                 NULL)
  if (fr==1) year <- year + 1 else {
    if (fr==mq) {  
    mq <- 1; year <- year + 1
    } else mq <- mq + 1               
  }
  window(data,start=c(year,mq))
}

#' @importFrom stats frequency start end
ts2df <- function(tsl,maxfreq) {
  if (missing(maxfreq)) 
    maxfreq <- if (is.list(tsl)) max(sapply(tsl,frequency)) else frequency(tsl)
  if (is.list(tsl)) lapply(tsl,ts2df,maxfreq=maxfreq) else {
    freq  <- frequency(tsl)
    if (freq==1) {
        startts <- as.Date(paste0(start(tsl)[1],"-07-01"))
        endts   <- as.Date(paste0(end(tsl)[1],"-07-01"))
      byts <- "year"
    }
    if (freq==4) {
      startts <- as.Date(paste0(start(tsl)[1],"-",(start(tsl)[2]-1)*3+2,"-15"))
      endts   <- as.Date(paste0(end(tsl)[1],"-",(end(tsl)[2]-1)*3+2,"-15"))
      byts <- "quarter"
    }
    if (freq==12) {
      startts <- as.Date(paste0(start(tsl)[1],"-",start(tsl)[2],"-15"))
      endts   <- as.Date(paste0(end(tsl)[1],"-",end(tsl)[2],"-15"))
      byts <- "month"
    }
    res <- cbind(seq.Date(startts,endts,byts),as.data.frame(tsl))
    colnames(res) <- c("date",
                       if (is.null(colnames(tsl))) "value" else colnames(tsl))
    res                   
  }
}

################################################################################
################################################################################
# Helper functions for OPTIMIZATION PROCEDURE                                  #
################################################################################
################################################################################

## generates 'trafo'-matrix A
genA <- function(len) 
  matrix(unlist(sapply(len,function(x) c(rep.int(1,x),rep(0,sum(len))))),
         nrow=sum(len))[,seq_along(len)]

## generates a trafo for v
genTrafo <- function(n.v=NULL,names.v=NULL,len.v) {
  has.trafo <- !(missing(len.v)||is.null(len.v)||all(len.v==1))
  if (has.trafo) {                                                              # for SCM*T*, 'real' trafo
    trafo     <- function(v) rep.int(v,len.v) 
    A         <- genA(len.v)
    n.v       <- length(len.v)
  } else {
    if (is.null(n.v)) 
      n.v     <- length(names.v)
    trafo     <- function(v) v
    A         <- diag(n.v)
    len.v     <- rep.int(1,n.v)
  }    
  list(n.v=n.v,names.v=names.v,trafo=trafo,A=A,has.trafo=has.trafo,len.v=len.v)
}

## eliminate and re-impute (near) zero elements of a vector
ess  <- function(x,lb=0) x[x>lb]
blow <- function(x,full.names) {
  stopifnot(isTRUE(all(names(x) %in% full.names)))
  res           <- rep(0,length(full.names))
  names(res)    <- full.names
  res[names(x)] <- x
  res
}

## unit vectors
e <- function(i,n) { tmp <- numeric(n); tmp[i] <- 1/length(i); tmp }

## fast calculation of V'MV
fastVpMV <- function(M,V) .Call(C_fastVpMV,M,V)

## fast calculation of V'M'MV
fastVpMpMV <- function(M,V) .Call(C_fastVpMpMV,M,V)

## fast calculation of M'diag(V)M
fastMpdVM <- function(M,V) .Call(C_fastMpdVM,M,V)

## calculation of loss for dependent variables
lossDep <- function(Z,w) fastVpMpMV(Z,w)/nrow(Z)

## calculation of loss for predictor variables
lossPred <- function(X,w,v,trafo.v) 
  if (missing(trafo.v)) abs(fastVpMV(fastMpdVM(X,v),w)) else 
                        abs(fastVpMV(fastMpdVM(X,trafo.v$trafo(v)),w))

## is vector w0 sunny wrt. X?
#' @importFrom lpSolve lp
isSunny <- function(w0,X) 
  isTRUE(all.equal(suppressWarnings(
    lp(direction="min",
       objective.in=c(rep(0,length(w0)),1),
       const.mat=rbind(c(rep(1,length(w0)),0),
                       cbind(X,-(X %*% w0))),
       const.rhs=c(1,rep(0,nrow(X))),
       const.dir=rep("==",nrow(X)+1)            
      ))$objval,1))

## calculate 'global' optimum and check for feasibility
checkGlobalOpt <- function(X,Z,trafo.v,lb,single.v=FALSE,verbose=FALSE,
                           debug=FALSE) {
  tmp   <- wnnlsGetGlobalOpt(Z)
  w     <- tmp$w
  rmspe <- tmp$rmspe
  if (exists_v(w,X,Z,trafo.v,lb)) {
    list(w=w,
         v=if (is.na(single.v)) 
           cbind("backup"=single_v(w,X,Z,trafo.v,lb,check.exists=TRUE)) else
           if (isTRUE(single.v))
             cbind("max.order"=single_v(w,X,Z,trafo.v,lb,verbose=verbose,
                                        debug=debug)) else 
             all_v(w,X,Z,trafo.v,lb,verbose=verbose,debug=debug),
         loss.v=rmspe^2,rmspe=rmspe,conv=0,single.v=single.v)
  } else list(w=w,v=NA,loss.v=rmspe^2,rmspe=rmspe,conv=NA,single.v=NA)
}

## calculate 'global' optimum with wnnls
wnnlsGetGlobalOpt <- function(X) {
  w <- wnnls(A=X,B=matrix(0,nrow=nrow(X)),E=matrix(rep(1,ncol(X)),nrow=1),
             F=matrix(1),check.input=FALSE)$x
  list(w=w,rmspe=sqrt(lossDep(X,w)))
}    

## try to improve outer target function (if solution of inner opt. is ambig.)   # not used by default
improveZw <- function(Z,X,w,tol=.Machine$double.eps) {
  w <- wnnlsInt(W=.Call(C_prepareW1,Z,X,w),ME=nrow(X)+1,MA=nrow(Z),
                N=length(w))$x
  w[w<tol] <- 0
  w/sum(w)
}

## solve outer optimization when there are no sunny donors
solveNoSunny <- function(X,Z,trafo.v,lb.loss=.Machine$double.eps,
                         lb=sqrt(.Machine$double.eps),verbose=FALSE) {
  n <- ncol(X)
  w <- wnnls(A=Z,B=matrix(0,nrow=nrow(Z)),E=rbind(X,rep(1,n)),
             F=matrix(c(rep(0,nrow(X)),1),ncol=1))$x
  names(w) <- colnames(X)             
  w[w<lb]  <- 0
  w <- w/sum(w)
  v <- rep(1/trafo.v$n.v,trafo.v$n.v)
  names(v) <- trafo.v$names.v
  loss.w <- lossPred(X,w,v,trafo.v)
  if (loss.w>lb.loss) 
    if (verbose) catn("Warning: loss.w is ",loss.w,">0!")
  list(w=w,loss.w=loss.w,v=cbind("max.order"=v),loss.v=lossDep(Z,w),
       rmspe=sqrt(lossDep(Z,w)),conv=0,single.v=TRUE)
}

## 'solve' outer optimization when there is exactly one sunny donor
solveOneSunny <- function(X,Z,trafo.v) {
  w        <- 1
  names(w) <- colnames(X)             
  v        <- rep(1/trafo.v$n.v,trafo.v$n.v)
  names(v) <- trafo.v$names.v
  loss.w   <- lossPred(X,w,v,trafo.v)
  list(w=w,loss.w=loss.w,v=cbind("max.order"=v),loss.v=lossDep(Z,w),
       rmspe=sqrt(lossDep(Z,w)),conv=0,single.v=TRUE)
}

################################################################################
################################################################################
# Helper functions for NICE OUTPUT                                             #
################################################################################
################################################################################

## format output as comma separated list with maximum line width and indentation
csl <- function(title,x,max.width=options()$width,min.indent=0,newline=FALSE,
                sep=", ") {
  n1  <- max(min.indent,nchar(title)+2)
  na  <- max(0,n1 - nchar(title) - 2)
  n2  <- max.width - n1
  nc  <- nchar(x)
  i   <- n1
  out <- paste0(title,if (title!="") ": " else "  ",
                paste0(rep(" ",na),collapse=""))
  for (j in seq_along(x)) {
    if ((i+nc[j]>max.width-1)||((j>1)&&newline)) {
      out <- paste0(out,"\n",paste0(rep(" ",n1),collapse=""),x[j],
                    if (j<length(x)) sep else "\n")
      i   <- n1 + nc[j] + nchar(sep)
    } else {
      out <- paste0(out,x[j],if (j<length(x)) sep else "\n")
      i   <- i + nc[j] + nchar(sep)
    }   
  }
  out
}

## format tab-like
tab <- function(nam,val,sep=": ") {
  nc <- max(nchar(nam))
  paste0(format(nam,width=nc,justify="left"),sep,val)
}

## wrap text with title
wrap <- function(title,x,max.width=options()$width,min.indent=0) {
  x <- unlist(strsplit(x," "))
  csl(title,x,max.width,min.indent,sep=" ")
}

## wrap matrices with titles (as argument names)
wrapM <- function(...,max.width=options()$width,min.indent=0,matrix.width=12) {
  out <- ""
  arg <- list(`...`)
  argnames <- names(arg)
  if (is.null(argnames)) argnames <- rep("",length(arg))
  op  <- options()
  n1  <- max(min.indent,max(nchar(argnames))+2)
  n2  <- max.width - n1
  options(width=n2)
  for (i in seq_along(arg)) 
    if (is.null(rownames(arg[[i]]))) 
      rownames(arg[[i]]) <- rep("",nrow(arg[[i]]))
  lrn <- max(sapply(arg,function(x) max(nchar(rownames(x)))))
  for (i in seq_along(arg)) {
    rownames(arg[[i]]) <- format(rownames(arg[[i]]),width=lrn)
    na  <- max(0,n1 - nchar(argnames[i]) - 2)
    out <- paste0(out,argnames[i],if (argnames[i]!="") ": " else "  ",
                  paste0(rep(" ",na),collapse=""),collapse="")
    if (any(is.character(arg[[i]]))) {
      tmp <- as.data.frame(arg[[i]])
    } else {
      colnames(arg[[i]]) <- format(colnames(arg[[i]]),width=matrix.width,
                                   justify="right")
      tmp <- format(arg[[i]],justify="right",width=matrix.width,na.encode=TRUE)
      tmp[which(tmp==format("NA",justify="right",width=matrix.width))] <- NA
    }                  
    tmp <- capture.output(print(tmp,na.print="",quote=FALSE))
    tmp <- paste0(tmp,sep=paste0("\n",paste0(rep(" ",n1),collapse=""),sep=""))
    out <- paste0(out,paste0(tmp,collapse=""),"\n",collapse="")
  }
  options(op)
  out
}

## wrap text with title (... version)
catw <- function(title,...) {
  x <- paste0(unlist(list(`...`)),collapse="")
  cat0(wrap(title,x))
}

## functions for pretty printing
now  <- function() paste0(format(Sys.time(),"%H:%M:%S"),": ")
now1 <- function() paste0(format(Sys.time(),"%H:%M:%S"),":")
cat0 <- function(...) cat(...,sep="")

## cat without separator and with fill option
catf <- function(...) cat(...,fill=TRUE,sep="")

## cat without separator with fill and label (current time) option
catn <- function(...) 
  cat(unlist(strsplit(paste0(unlist(list(`...`)),collapse="")," ")),
      fill=TRUE,sep=" ",labels=now1())

## align
align <- function(...) {
  obj <- list(`...`)
  nc  <- sapply(obj,function(x) max(nchar(format(x))))
  tmp <- vector("list",length(obj))
  for (i in seq_along(tmp))
    tmp[[i]] <- format(obj[[i]],width=nc[i],justify="left")
  do.call("paste0",tmp)
}


################################################################################
################################################################################
# Helper functions for LINEAR PROGRAMS                                         #
################################################################################
################################################################################

#' @importFrom lpSolveAPI lp.control make.lp set.rhs set.objfn set.column
#' @importFrom lpSolveAPI set.constr.type set.bounds get.objective get.variables
lp3 <- function(direction,objective.in,const.mat,const.rhs,const.dir,
                #vmin=0,vmax=1,
                ...) {
  n <- ncol(const.mat)
#  const.mat <- const.mat[const.dir!="<=",]
#  const.rhs <- const.rhs[const.dir!="<="]
#  const.dir <- const.dir[const.dir!="<="]
  const.dir[const.dir=="=="] <- "="
  mylp <- make.lp(nrow(const.mat),n)
  set.rhs(mylp,const.rhs)
  set.constr.type(mylp,const.dir)
  for (i in 1:n) set.column(mylp,i,const.mat[,i])
#  set.bounds(mylp,lower=rep(vmin,n),upper=rep(vmax,n))
  set.objfn(mylp,objective.in)
  lp.control(mylp,presolve=c("lindep","rows","rowdominate"),
             sense=direction,timeout=10,verbose="neutral")
  status <- solve(mylp)
  list(solution=if (status==0) get.variables(mylp) else NA, status=status,
       objval=if (status==0) get.objective(mylp) else NA)
}

#' @importFrom Rglpk Rglpk_solve_LP
#' @importFrom lpSolve lp
lpsolve <- function(direction,objective.in,const.mat,const.rhs,const.dir,
                    not.fixed=NULL,scale=196L,tol.lp=1e-13,
                    verbose=TRUE,debug=FALSE,check.min=FALSE,...) {
  tol_mul <- 1 + sqrt(.Machine$double.eps)
  rglpk <- Rglpk_solve_LP(obj=objective.in,mat=const.mat,rhs=const.rhs,         # 1st try: solve lp via Rglpk
                          dir=const.dir,max=isTRUE(direction=="max"),
                          control=list(tm_limit=5000,presolve=TRUE))
  rglpk$objval <- rglpk$optimum
  dglpk <- lp_diag(objective.in,const.mat,const.rhs,const.dir,rglpk)
  mglpk <- if (check.min)                                                       # variable defined as minimum indeed minimal?
             isTRUE(rglpk$solution[length(rglpk$solution)] <=
                      min(rglpk$solution[not.fixed]) * tol_mul) else TRUE
  sglpk <- isTRUE(rglpk$status==0)
  score.glpk <- (!mglpk) + (!sglpk) + max(dglpk[2:3])
  if (score.glpk > tol.lp) {                                                    # minimum not ok or status not ok or constraints violated?
    rlpA <- lp3(objective.in,const.mat=const.mat,const.rhs=const.rhs,           # 2nd try: solve lp via lpSolveAPI
                const.dir=const.dir,direction=direction,scale=scale)
    dlpA <- lp_diag(objective.in,const.mat,const.rhs,const.dir,rlpA)
    mlpA <- if (check.min)                                                      # variable defined as minimum indeed minimal?
               isTRUE(rlpA$solution[length(rlpA$solution)] <=
                        min(rlpA$solution[not.fixed]) * tol_mul) else TRUE
    slpA <- isTRUE(rlpA$status==0)
    score.lpA <- (!mlpA) + (!slpA) + max(dlpA[2:3])
    if (score.lpA > tol.lp) {                                                   # minimum not ok or status not ok or constraints violated?
      rlp <- lp(objective.in,const.mat=const.mat,const.rhs=const.rhs,           # 3rd try: solve lp via lpSolve
                const.dir=const.dir,direction=direction,scale=scale)
      dlp <- lp_diag(objective.in,const.mat,const.rhs,const.dir,rlp)
      mlp <- if (check.min)                                                     # variable defined as minimum indeed minimal?
                 isTRUE(rlp$solution[length(rlp$solution)] <=
                          min(rlp$solution[not.fixed]) * tol_mul) else TRUE
      slp <- isTRUE(rlp$status==0)
      score.lp <- (!mlp) + (!slp) + max(dlp[2:3])
      if (score.lp > tol.lp) {                                                  # minimum not ok or status not ok or constraints violated?
        if (min(score.glpk,score.lpA,score.lp)>100*tol.lp) {
          if (debug) 
            catn("lp solvers Rglpk, lpSolveAPI and lpSolve failed severely. score: ", min(score.glpk,score.lpA,score.lp))
        } else if (debug) catn("lp solvers failed.")
        if (debug) print(rbind("glpk"=dglpk,"lpSolveAPI"=dlpA,"lpSolve"=dlp))
      } 
      switch(which.min(c(score.glpk,score.lpA,score.lp)),
              "1" = rglpk, "2" = rlpA, "3" = rlp)
    } else rlpA
  } else rglpk
}

## linear programs for ratios in target function
lprsolve <- function(direction = "min", a0, a, d0,d, const.mat, const.dir, 
                     const.rhs, tol.lp=1e-13, scale=196L,
                     verbose=FALSE,debug=FALSE,...) {
  const.mat    <- cbind(c(-const.rhs,d0),rbind(const.mat,d))
  objective.in <- c(a0,a)
  const.dir    <- c(const.dir,"==")
  const.rhs    <- c(rep(0,length(const.rhs)),1)
  tmp <- lpsolve(objective.in=objective.in,const.mat=const.mat,
                 const.dir=const.dir,const.rhs=const.rhs,direction=direction,
                 tol.lp=tol.lp,check.min=FALSE,scale=scale,verbose=verbose,
                 debug=debug,...)
  if (tmp$status==0) {
    if (tmp$solution[1]>0) tmp$solution[-1]/tmp$solution[1] else 
                           rep(0,length(a))
  } else rep(NA,length(a))
}

# new version of PUFAS which allows for inequality constraints as well as 
#   checking for unique optimal solution (default) and for unique feasible 
#   solution (by setting optimal=FALSE)
pufas <- function(objective.in,const.mat,const.rhs,const.dir,solution,tol=0,
                  do.check=FALSE,tol.check=tol,tol.cond=.Machine$double.eps,
                  optimal=TRUE,scale=196L,verbose=FALSE,debug=FALSE,...) {
  # x is an optimal basic (!) solution of max objective x 
  # s.t. const.mat x = const.rhs, x >= 0
  # Appa (2002): solve linear program max d x s.t. const.mat x = const.rhs, 
  #              objective x = objective solution, x >= 0, with d = 1 for 
  #              zero components of solution and 0 else
  
  if (do.check && !(rcond(const.mat[,solution>tol.check])>tol.cond)) return(NA) # check whether given solution is a basis solution, return NA else

  eqs       <- const.dir %in% c("=","==")
  ineq_less <- const.dir %in% c("<","<=")
  ineq_more <- const.dir %in% c(">",">=")
  mat_slack <- const.mat
  mat_slack[eqs,] <- 0
  mat_slack[ineq_more,] <- -mat_slack[ineq_more,]
  slack_zero      <- rep(NA,nrow(const.mat))
  slack_zero[eqs] <- FALSE
  slack_zero[ineq_less] <- ( const.mat[ineq_less,] %*% 
                                              solution >= const.rhs[ineq_less] )
  slack_zero[ineq_more] <- ( const.mat[ineq_more,] %*% 
                                              solution <= const.rhs[ineq_more] )
  
  obj.new <- 1*(solution<=tol)-colSums(mat_slack[slack_zero,,drop=FALSE])
  
  if (optimal) {
    mat.new <- rbind(const.mat,objective.in)
    new.rhs <- c(const.rhs,sum(objective.in*solution))
    if (length(const.dir)==1) const.dir <- rep(const.dir,nrow(const.mat))
    new.dir <- c(const.dir,"==")
  } else {
    mat.new <- const.mat
    new.rhs <- const.rhs
    new.dir <- const.dir
  }

  new_lp_res <- lpsolve(objective.in=obj.new,const.mat=mat.new,
                        const.rhs=new.rhs,const.dir=new.dir,direction="max",
                        verbose=verbose,debug=debug,scale=scale,...)

# if (new_lp_res$status==3) return(FALSE)                                       # if new_lp is unbounded, then there are many alternative solutions, FIX ME!!!
  if (new_lp_res$status!=0) return(NA)
  new_slack_zero      <- rep(NA,nrow(const.mat))
  new_slack_zero[eqs] <- FALSE
  new_slack_zero[ineq_less] <- ( const.mat[ineq_less,] %*% 
                                   new_lp_res$solution >= const.rhs[ineq_less] )
  new_slack_zero[ineq_more] <- ( const.mat[ineq_more,] %*% 
                                   new_lp_res$solution <= const.rhs[ineq_more] )
  all(new_lp_res$solution[solution<=tol]<=tol) && 
    all( !slack_zero | new_slack_zero)
}

################################################################################
################################################################################
# Helper functions for evaluating and improving solutions of linear programs   #
################################################################################
################################################################################
  
max_error <- function(a,b,subset=NULL,direction="==",correct="none") {
  if (!is.null(subset)) { a <- a[subset]; b <- b[subset] }
  if (length(a)>0) {
    abserr <- switch(direction,
                     "==" = abs(a-b),"<=" = pmax(a-b,0),">=" = pmax(b-a,0))
    base   <- switch(correct,
                     "none" = pmax(abs(a),abs(b)),"a"=abs(a),"b"=abs(b))
    relerr <- ifelse(base==0,Inf,abserr/base)
    max(pmin(abserr,relerr))
  } else 0
}

lp_diag <- function(objective.in,const.mat,const.rhs,const.dir,lpres) {
  if (is.null(lpres$solution)) {
    sapply(lpres, function(x) lp_diag(objective.in,const.mat,const.rhs,
                                      const.dir,x$solution)) 
  } else {
    if ((lpres$status==0)&&(length(lpres$solution)==ncol(const.mat))) {
      obj     <- sum(objective.in*lpres$solution)
      lhs     <- const.mat%*%lpres$solution
      eqmax   <- max_error(lhs,const.rhs,subset=(const.dir=="=="))
      ineqmax <- max(max_error(lhs,const.rhs,subset=(const.dir=="<="),
                               direction="<="),
                     max_error(lhs,const.rhs,subset=(const.dir==">="),
                               direction=">="))
      minimum <- min(lpres$solution[-length(lpres$solution)])
      target  <- lpres$solution[length(lpres$solution)]
    } else obj <- lhs <- eqmax <- ineqmax <- minimum <- target <- Inf  
    c("objective"=obj,"eqmaxinf"=eqmax,"ineqmaxinf"=ineqmax,
      "status"=lpres$status,"minimum"=minimum,"target"=target)
  }
}

sanitize_w <- function(w,tol=1e-14) {
  if (is.matrix(w)) apply(w,2,function(x) sanitize_w(x,tol=tol)) else {
    w[w<tol] <- 0
    w/sum(w)
  }
}

sanitize_v <- function(v,X,Z,trafo.v,lb=1e-8,rel.tol=1e-3,tol.lpq=1e-13,
                       tol.loss=1e-12,verbose=FALSE,debug=FALSE) {
  v <- v/max(v)
  check_lp <- function(v,w) {
    M0   <- M(w,X,trafo.v)
    K    <- ncol(M0)
    J    <- nrow(M0)
    const.mat <- M0
    const.rhs <- rep(0,nrow(M0))
    const.dir <- ifelse(w==0,">=","==")
    lhs     <- const.mat%*%v
    eqmax   <- max(c(0,abs(lhs-const.rhs)[const.dir=="=="]))
    ineqmax <- max(c(0,(lhs-const.rhs)[const.dir=="<="],
                       (const.rhs-lhs)[const.dir==">="]))
    max(eqmax,ineqmax)
  }
  
  if (any(is.na(v))) return(v)
  oldv.w <- wnnlsExt(v,X,NULL,trafo.v,return.w=TRUE)
  oldv.l <- lossDep(Z,oldv.w)
  newv <- v
  newv[v<lb*(1+rel.tol)] <- lb
  newv[v>1*(1-rel.tol)]  <- 1
  newv.w <- wnnlsExt(newv,X,NULL,trafo.v,return.w=TRUE)
  newv.l <- lossDep(Z,newv.w)

  if (debug) catn("sanitize_v: old v: ",oldv.l,", new v: ",newv.l,
                  ", change new vs. old: ",newv.l-oldv.l)

  if ((oldv.l>0)&&(min((newv.l-oldv.l)/oldv.l,newv.l-oldv.l)<tol.loss)&&
      (check_lp(newv,newv.w)-check_lp(v,oldv.w)<tol.lpq)) newv else v
}

sanitize_res <- function(res,X,Z,trafo.v,lb=1e-8,rel.tol=1e-3,tol.lpq=1e-13,
                         tol.loss=1e-12,verbose=FALSE,debug=FALSE) {
  res$v <- res$v/max(res$v)
  check_lp <- function(v,w) {
    M0   <- M(w,X,trafo.v)
    K    <- ncol(M0)
    J    <- nrow(M0)
    const.mat <- M0
    const.rhs <- rep(0,nrow(M0))
    const.dir <- ifelse(w==0,">=","==")
    lhs     <- const.mat%*%v
    eqmax   <- max(c(0,abs(lhs-const.rhs)[const.dir=="=="]))
    ineqmax <- max(c(0,(lhs-const.rhs)[const.dir=="<="],
                       (const.rhs-lhs)[const.dir==">="]))
    max(eqmax,ineqmax)
  }
  
  if (any(is.na(res$v))) return(res)
  wnames <- names(res$w)
  res$w  <- wnnlsExt(res$v,X,NULL,trafo.v,return.w=TRUE)
  res$loss.v <- lossDep(Z,res$w)
  res$rmspe  <- sqrt(res$loss.v)
  newv <- res$v
  newv[res$v<lb*(1+rel.tol)] <- lb
  newv[res$v>1*(1-rel.tol)]  <- 1
  newv.w <- wnnlsExt(newv,X,NULL,trafo.v,return.w=TRUE)
  newv.l <- lossDep(Z,newv.w)

  if (debug) catn("sanitize_res: old v: ",res$loss.v,", new v: ",newv.l,
                  ", change new vs. old: ",newv.l-res$loss.v)

  if ((res$loss.v>0)&&
      (min((newv.l-res$loss.v)/res$loss.v,newv.l-res$loss.v)<tol.loss)&&
      (check_lp(newv,newv.w)-check_lp(res$v,res$w)<tol.lpq)) {
    res$w <- newv.w
    res$v <- newv
    res$loss.v <- newv.l
    res$rmspe  <- sqrt(res$loss.v)
  }
  names(res$w) <- wnames
  res
}

vdiff <- function(v) {
  n <- ncol(v)
  if (n==1) return(0)
  idx   <- expand.grid(1:n,1:n)
  idx   <- idx[idx[,1]<idx[,2],,drop=FALSE]
  diffs <- apply(idx,1,function(x) 
                         max(pmax(abs(v[,x[1]]-v[,x[2]]),
                             abs(v[,x[1]]-v[,x[2]])/(v[,x[1]]+v[,x[2]]))))
  max(diffs)
}

wdiff <- function(v) {
  n <- ncol(v)
  if (n==1) return(0)
  idx   <- expand.grid(1:n,1:n)
  idx   <- idx[idx[,1]<idx[,2],,drop=FALSE]
  diffs <- apply(idx,1,function(x) max(abs(v[,x[1]]-v[,x[2]])))
  max(diffs)
}

bestifnear <- function(v,X,Z,trafo.v,tol=1e-3) {
  if(vdiff(v)<tol) {
    losses <- apply(v,2,function(x) {
      w <- wnnlsExt(x,X,NULL,trafo.v,return.w=TRUE)
      lossDep(Z,w)
    })
    v[,which.min(losses)]
  } else v[,1]
}
