## R interface to the WNNLS Fortran algorithm 
## Fortran code from package limSolve, with kind permission of Karline Soetaert
## Code is an adapted version of function lsei of package limSolve
## PROGOPT not yet supported
## documenting and exporting this function may be useful...
## to do(?): @export wnnls
wnnls <- function(A=NULL, B=NULL, E=NULL, F=NULL, tol=sqrt(.Machine$double.eps),
                  L=0, verbose=TRUE, check.input=TRUE)  {

  # check input
  if (check.input) {
    if (is.vector(E) & length(F)==1)          E <- t(E)
      else if (! is.matrix(E) & ! is.null(E)) E <- as.matrix(E)
    if (is.vector(A) & length(B)==1)          A <- t(A)
      else if (! is.matrix(A) & ! is.null(A)) A <- as.matrix(A)
    if (! is.matrix(F) & ! is.null(F))        F <- as.matrix(F)
    if (! is.matrix(B) & ! is.null(B))        B <- as.matrix(B)
    if (is.null(A) && is.null (E)) {
      stop("there is no problem to solve: A and E are NULL")
    } else if (is.null(A)) {
      A <- E[1,,drop=FALSE]
      B <- F[1]
    }
  }
  
  # dimensions
  ME  <- nrow(E)
  MA  <- nrow(A)
  MDW <- ME + MA
  N   <- if (!is.null(ncol(A))) ncol(A) else ncol(E)

  if (check.input) {
    # No equalities? -> stop (and hint to nnls), check dimensions otherwise
    if (is.null (ME)) {
      stop("no equalities, please use function nnls in package nnls instead")
    } else  {
      if (ncol(E)   != N)  stop("wrong input: A and E are not compatible")
      if (length(F) != ME) stop("wrong input: E and F are not compatible")
    }
    if (length(B) !=  MA) stop("wrong input: A and B are not compatible")
  }
  
  # Initialize workspace and big matrix input for Fortran code
  IWORK      <- integer(ME+MA+N)
  WORK       <- numeric(ME+MA+5*N)
  IWORK[1:2] <- c(length(WORK),length(IWORK))
  W          <- rbind(cbind(E,F),cbind(A,B))

  # obtain solution
  sol <-.Fortran(C_wnnls,W=as.double(W),MDW=as.integer(MDW),ME=as.integer(ME),
                 MA=as.integer(MA),N=as.integer(N),L=as.integer(L),
                 PRGOPT=as.double(1.0),X=double(N),RNORM=double(1),
                 MODE=integer(1),IWORK=IWORK,WORK=WORK)

  Error <- any(is.infinite(sol$X)) || (sol$MODE>0) 

  X        <- sol$X
  names(X) <- if (!is.null(colnames(A))) colnames(A) else colnames(E)

  # prepare residuals/solution
  if (any(is.infinite(X))) {
    residual <- solution <- Inf
  } else  {
    # residual (l1) norm
    residual <- if (ME > 0) sum(abs(E %*% X - F)) else 0
    if (residual>tol) Error <- TRUE
    # solution value (l2 norm)
    solution <- if (MA > 0) sum ((A %*% X - B)^2) else 0
  }

  list(x=X,value=solution,residual=residual,error=Error,status=sol$MODE)
}

## internally used version of wnnls (no input checks!)
wnnlsInt <- function(W, ME, MA, N, tol=sqrt(.Machine$double.eps))  {
  MDW        <- ME + MA
  IWORK      <- integer(MDW+N)
  WORK       <- double(MDW+5*N)
  IWORK[1:2] <- c(length(WORK),length(IWORK))
  sol <- .Fortran(C_wnnls,W=W,MDW=as.integer(MDW),ME=as.integer(ME),
                  MA=as.integer(MA),N=as.integer(N),L=0L,PRGOPT=as.double(1.0),
                  X=double(N),RNORM=double(1),MODE=integer(1),IWORK=IWORK,
                  WORK=WORK)
  Error <- any(is.infinite(sol$X)) || (sol$MODE>0) 
  list(x=sol$X,error=Error,status=sol$MODE)
}
