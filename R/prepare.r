## internal function, export later?
## todo (?): @export prepare
## user interface made for compatibility with function dataprep of package Synth
## use listFromLong to convert 'long' format to list (as required by prepare)
#' @importFrom stats var
prepare <- function(dat, predictors, predictors.op = "mean", special.predictors,
                    dependent, treatment.identifier, controls.identifier, 
                    time.predictors.prior, time.optimize.ssr,
                    alpha=NULL, beta=NULL, gamma=NULL, scale.Z=FALSE)
{
  if (!is.list(dat)) stop("dat must be a list, consider using listFromLong")

  # main workhorse
	BuildMatrix <- function(dat,vars,times,donors,treated,aggfuns,alpha=NULL,
                          betagamma=NULL,scale=FALSE) {
		Select <- function(lst,nam,tim,cols,agg=NULL) {
			tmp <- lst[[nam]][intersect(tim,rownames(lst[[nam]])),cols,drop=FALSE]
            if (is.null(agg)||agg=="id") {                                      # w/o aggregation
                sapply(1:ncol(tmp),function(i) if (any(is.na(tmp[,i]))) 
                  stop(paste0(cols[i]," has some NAs in variable ",nam)))
                rownames(tmp) <- paste(nam,rownames(tmp),sep=".")
            } else {                                                            # with aggregation
                sapply(1:ncol(tmp),function(i) if (all(is.na(tmp[,i]))) 
                  stop(paste0(cols[i]," has only NAs in variable ",nam)))
                sapply(1:ncol(tmp),function(i) if (any(is.na(tmp[,i]))) 
                  warning(paste0(cols[i]," has some NAs in variable ",nam)))
                tmp <- apply(tmp,2,agg,na.rm=TRUE)
                if (!is.matrix(tmp)) tmp <- matrix(tmp,nrow=1)
                rownames(tmp) <- paste(nam,agg,tim[1],tim[length(tim)],sep=".")
            }    
			tmp
		}		
		if (!is.null(alpha)) stopifnot(length(alpha)==length(vars))
		if (!is.null(betagamma)) {
      stopifnot(length(betagamma)==length(vars))
      if (!is.list(betagamma)) betagamma <- as.list(betagamma)
    }  
		X0X1     <- NULL
    len      <- NULL
    names.v  <- NULL
    divisors <- NULL
    if ((length(times)==1)&&(length(vars)>1)) 
      for (i in 2:length(vars)) times[[i]] <- times[[1]]
		for (i in seq_along(vars)) {
      if (is.list(times[[i]])) {
        if (!is.null(betagamma)) 
          if (sum(sapply(times[[i]],length))!=length(betagamma[[i]])) 
            stop("beta/gamma probably has incorrect length")
        tmp <- NULL
        for (j in seq_along(times[[i]])) 
          tmp <- rbind(tmp,Select(dat,vars[i],times[[i]][[j]],c(donors,treated),
                                  if(missing(aggfuns)) NULL else aggfuns[i]))
      } else tmp <- Select(dat,vars[i],times[[i]],c(donors,treated),
                           if(missing(aggfuns)) NULL else aggfuns[i])
      if (scale) {
        divisor  <- sqrt(var(as.vector(tmp))*nrow(tmp))
        divisors <- c(divisors,divisor)
        tmp <- tmp/divisor
      }    
      if ((length(betagamma[[i]])>1)&&(length(betagamma[[i]])!=nrow(tmp)))
        stop("beta/gamma probably has incorrect length")
      if (!is.null(betagamma))   
        tmp <- tmp * if(length(betagamma[[i]])==1) 
                       betagamma[[i]]^seq(nrow(tmp)-1,0,by=-1) else
                       rev(betagamma[[i]])
      if (!is.null(alpha)) tmp <- tmp * sqrt(alpha[[i]])                       
      X0X1    <- rbind(X0X1,tmp)
      len     <- c(len,nrow(tmp))
      names.v <- c(names.v,
                   if (missing(aggfuns)||(aggfuns[i]=="id")) vars[i] else
                     if (nrow(tmp)==1) rownames(tmp) else vars[i])
    }
    colnames(X0X1) <- c(donors,treated)
    list(X0=X0X1[,-ncol(X0X1),drop=FALSE],X1=X0X1[,ncol(X0X1),drop=FALSE],
#         info.v=list(n.v=length(len),names.v=names.v,len.v=len),
         trafo.v=genTrafo(n.v=length(len),names.v=names.v,len.v=len),
         scaled=scale,divisors=divisors)
	}
    
  # check treated unit
  if (length(treatment.identifier) != 1 ) 
    stop("please specify exactly one treatment unit\n")

  # check control units
  if (length(controls.identifier) < 2 ) 
    stop("\n please specify at least two control units\n")
  if (any(duplicated(controls.identifier))) 
    stop("\n duplicate control units in controls.identifier\n")

  # more checks
  if (treatment.identifier %in% controls.identifier) 
    stop("\n treated unit among controls\n")

  allpredictors <- c(predictors,sapply(special.predictors,function(x) x[[1]]))
  allpredtimes  <- vector("list",length(predictors))
  for (i in seq_along(predictors)) 
    allpredtimes[[i]] <- time.predictors.prior
  for (i in seq_along(special.predictors)) 
    allpredtimes <- c(allpredtimes,list(special.predictors[[i]][[2]]))
  allpredaggfun <- rep(predictors.op,length(predictors))
  for (i in seq_along(special.predictors)) 
    allpredaggfun <- c(allpredaggfun,special.predictors[[i]][[3]])

  if (!is.list(time.optimize.ssr)) time.optimize.ssr <- list(time.optimize.ssr)
    
	tmpX <- BuildMatrix(dat,allpredictors,allpredtimes,controls.identifier,
                      treatment.identifier,allpredaggfun,scale=TRUE,
                      betagamma=gamma)
	tmpZ <- BuildMatrix(dat,dependent,time.optimize.ssr,controls.identifier,
                      treatment.identifier,scale=scale.Z,alpha=alpha,
                      betagamma=beta)
	tmpZu <- BuildMatrix(dat,dependent,time.optimize.ssr,controls.identifier,
                       treatment.identifier,scale=FALSE,alpha=alpha,
                       betagamma=beta)

  res <- list(X0=tmpX$X0,X1=tmpX$X1,Z0=tmpZ$X0,Z1=tmpZ$X1,trafo.v=tmpX$trafo.v,
              Z.scaled=scale.Z)
  if (scale.Z) res <- c(res,
                        list(Z0u=tmpZu$X0,Z1u=tmpZu$X1,
                        Z.len=tmpZu$trafo.v$len.v))
  res
}
