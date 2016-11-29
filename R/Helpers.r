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
checkGlobalOpt <- function(X,Z,trafo.v,lb,single.v=FALSE) {
  tmp   <- wnnlsGetGlobalOpt(Z)
  w     <- tmp$w
  rmspe <- tmp$rmspe
  if (exists_v(w,X,trafo.v,lb)) {
    list(w=w,v=if (single.v) single_v(w,X,trafo.v,lb) else 
                             all_v(w,X,trafo.v,lb),
         loss.v=rmspe^2,rmspe=rmspe,conv=0,single.v=single.v)
  } else list(w=w,v=NA,loss.v=rmspe^2,rmspe=rmspe,conv=NA,single.v=NA)
}

## calculate 'global' optimum with wnnls
wnnlsGetGlobalOpt <- function(X) {                                              # TO DO(?): check whether w is unique or not (-> SK)
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
## @importFrom limSolve lsei
solveNoSunny <- function(X,Z,trafo.v,method=c("wnnls","lsei"),                  # was: method=c("lsei","wnnls"),
                         lb.loss=.Machine$double.eps,
                         lb=sqrt(.Machine$double.eps),verbose=FALSE) {
  method <- match.arg(method)                           
  n <- ncol(X)
  w <- switch(method,
              "lsei"  = lsei(A=Z,B=matrix(0,nrow=nrow(Z)),E=rbind(X,rep(1,n)),  # more stable than 'wnnls' (for *this* application). EDIT: Really???
                             F=matrix(c(rep(0,nrow(X)),1),ncol=1),
                             G=diag(n),H=matrix(0,nrow=n))$x,
              "wnnls" = wnnls(A=Z,B=matrix(0,nrow=nrow(Z)),E=rbind(X,rep(1,n)),
                              F=matrix(c(rep(0,nrow(X)),1),ncol=1))$x
       )
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

## linear programs for ratios in target function
#' @importFrom lpSolve lp
lpr <- function(direction = "min", a0, a, d0,d, const.mat, const.dir, 
                const.rhs, ...) {
  tmp <- lp(direction=direction,objective.in=c(a0,a),
            const.mat=cbind(c(-const.rhs,d0),rbind(const.mat,d)),
            const.dir=c(const.dir,"=="),
            const.rhs=c(rep(0,length(const.rhs)),1),...)
  if (tmp$status==0) tmp$solution[-1]/tmp$solution[1] else rep(NA,length(a))
}       


# new version of PUFAS which allows for inequality constraints as well as 
#   checking for unique optimal solution (default) and for unique feasible 
#   solution (by setting optimal=FALSE)
#' @importFrom lpSolve lp
PUFAS <- function(objective,const.mat,const.rhs,const.dir,solution,tol=0,
                  do.check=FALSE,tol.check=tol,tol.cond=.Machine$double.eps,
                  optimal=TRUE,...) {
	# x is an optimal basic (!) solution of max objective x 
	# s.t. const.mat x = const.rhs, x >= 0
	# Appa (2002): solve linear program max d x s.t. const.mat x = const.rhs, 
	#              objective x = objective solution, x >= 0, with d = 1 for 
	#              zero components of solution and 0 else
	
	if (do.check && !(rcond(const.mat[,solution>tol.check])>tol.cond)) return(NA) # check whether given solution is a basis solution, return NA else

	eqs <- const.dir %in% c("=","==")
	ineq_less <- const.dir %in% c("<","<=")
	ineq_more <- const.dir %in% c(">",">=")
	mat_slack <- const.mat
	mat_slack[eqs,] <- 0
	mat_slack[ineq_more,] <- -mat_slack[ineq_more,]
	slack_zero <- rep(NA,nrow(const.mat))
	slack_zero[eqs] <- FALSE
	slack_zero[ineq_less] <- ( const.mat[ineq_less,] %*% 
                                              solution >= const.rhs[ineq_less] )
	slack_zero[ineq_more] <- ( const.mat[ineq_more,] %*% 
                                              solution <= const.rhs[ineq_more] )
	
	obj.new <- 1*(solution<=tol)-colSums(mat_slack[slack_zero,,drop=FALSE])
	
	if (optimal)
	{
		mat.new <- rbind(const.mat,objective)
		new.rhs <- c(const.rhs,sum(objective*solution))
		if (length(const.dir)==1) const.dir <- rep(const.dir,nrow(const.mat))
		new.dir <- c(const.dir,"==")
	} else
	{
		mat.new <- const.mat
		new.rhs <- const.rhs
		new.dir <- const.dir
	}
	new_lp_res <- lp(direction="max",objective.in=obj.new,const.mat=mat.new,
						const.rhs=new.rhs,const.dir=new.dir,...)
	# print(new_lp_res$solution)			   
	if (new_lp_res$status==3) return(FALSE)                                       # if new_lp is unbounded, then there are many alternative solutions
	if (new_lp_res$status!=0) return(NA)
	# max(abs(solution-new_lp_res$solution))<=tol                                 # return TRUE if no new solution has been found
	new_slack_zero <- rep(NA,nrow(const.mat))
	new_slack_zero[eqs] <- FALSE
	new_slack_zero[ineq_less] <- ( const.mat[ineq_less,] %*% 
                                   new_lp_res$solution >= const.rhs[ineq_less] )
	new_slack_zero[ineq_more] <- ( const.mat[ineq_more,] %*% 
                                   new_lp_res$solution <= const.rhs[ineq_more] )
	all(new_lp_res$solution[solution<=tol]<=tol) && 
    all( !slack_zero | new_slack_zero)
}
  