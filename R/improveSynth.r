#' Check (and Improve) Results of Package Synth
#'
#' \code{improveSynth} checks the results of \code{synth}
#' for feasibility and optimality and tries to find a better solution.
#'
#' Performing SCM means solving a nested optimization problem. Depending on
#' the validity of the results of the inner optimization, SCM may produce
#' \itemize{
#' \item invalid or infeasible results, if the vector \code{w} of donor 
#' weights reported as the result of the inner optimization 
#' is in fact not optimal, ie. produces too large \code{loss.w},
#' \item suboptimal results, if the vector \code{v} of predictor weights 
#' reported as the result of the outer optimization is in fact not
#' optimal (which may be caused by shortcomings of the inner optimization).
#' }
#'
#' \code{improveSynth} first checks \code{synth.out} for feasibility and
#' then tries to find a feasible and optimal solution by applying the
#' optimization methods of package \code{MSCMT} to \code{dataprep.out} 
#' (with default settings, more flexibility will probably be added in a 
#' future release).
#'
#' @param synth.out A result of \code{synth} from package 
#' \code{'Synth'}.
#' @param dataprep.out The input of function \code{synth} which
#' led to \code{synth.out}.
#' @param lb A numerical scalar (default: \code{1e-8}), corresponding to the
#' lower bound for the outer optimization.
#' @param tol A numerical scalar (default: \code{1e-5}). If the relative *and*
#' absolute improvement of \code{loss.v} and \code{loss.w}, respectively,
#' exceed \code{tol}, this is reported 
#' (if \code{verbose} is \code{TRUE}). Better optima are always looked for 
#' (independent of the value of \code{tol}).
#' @param verbose A logical scalar. Should the ouput be verbose (defaults to 
#' \code{TRUE}).
#' @param seed A numerical vector or \code{NULL}. See the corresponding
#' documentation for \code{\link[MSCMT]{mscmt}}. Defaults to 1 in order to
#' provide reproducibility of the results. 
#' @param ... Further arguments to \code{\link[MSCMT]{mscmt}}. Supported 
#' arguments are \code{check.global}, \code{inner.optim}, \code{inner.opar}, 
#' \code{outer.optim}, and \code{outer.opar}.
#' @return An updated version of \code{synth.out}, where \code{solution.v},
#' \code{solution.w}, \code{loss.v}, and \code{loss.w} are replaced by the
#' optimum obtained by package \code{'MSCMT'} and all other components 
#' of \code{synth.out} are removed.
#' @export improveSynth
#' @examples 
#' \dontrun{
#' ## example has been removed because package 'Synth' has been archived
#' ## See vignette 'Checking and Improving Results of package Synth'
#' ## for an example working with a cached copy
#' }
#' @importFrom stats sd
improveSynth <- function(synth.out,dataprep.out,lb=1e-8,tol=1e-5,
                         verbose=TRUE,seed=1,...) {
  storage.mode(dataprep.out$X0) <- storage.mode(dataprep.out$X1) <- 
    storage.mode(dataprep.out$Z0) <- storage.mode(dataprep.out$Z1) <- "double"
  v   <- as.numeric(synth.out$solution.v)
  Xu  <- cbind(dataprep.out$X0, dataprep.out$X1)                                # generate (scaled!) X from dataprep object
  Xs  <- Xu/apply(Xu, 1, sd)
  X   <- Xs[, -ncol(Xs)] - drop(Xs[, ncol(Xs)])
  Z   <- dataprep.out$Z0 - drop(dataprep.out$Z1)
  ME  <- 1L
  MA  <- length(v)
  N   <- ncol(X)
  MDW <- ME+MA
  globals$Ipar  <- as.integer(c(ME=ME,MA=MA,MDW=MDW,N=N)) 
  globals$IWORK <- integer(MDW+N)
  globals$WORK  <- double(MDW+5*N)
  globals$IWORK[1:2] <- c(length(globals$WORK),length(globals$IWORK))
  globals$RNORM <- double(1)
  globals$MODE  <- integer(1)
  globals$X     <- double(N)
  tv  <- as.numeric(v)
  sol <-.Fortran(C_wnnls,W=.Call(C_prepareW4,X,tv),
                 MDW=globals$Ipar[3],ME=globals$Ipar[1],MA=globals$Ipar[2],
                 N=globals$Ipar[4],L=0L,PRGOPT=1.0,X=globals$X,
                 RNORM=globals$RNORM,MODE=globals$MODE,IWORK=globals$IWORK,
                 WORK=globals$WORK)
  if ((any(is.infinite(sol$X))) || (sol$MODE>0)) 
    warning("error in inner optimization (wnnls)")                
  w <- sol$X
  if (verbose) {
    cat0("Results reported by package Synth\n",
         "=================================\n\n")
    catw("Optimal V    ",paste0(synth.out$solution.v,sep=" "))
    catw("Optimal W*(V)",paste0(synth.out$solution.w,sep=" "))
    catf("with corresponding predictor loss ('loss W') of ",
         as.numeric(synth.out$loss.w)," ",
         "and corresponding dependent loss ('loss V') of ",
         as.numeric(synth.out$loss.v),".\n\n")
  }         
  
  w.orig <- as.numeric(synth.out$solution.w)
  if ((sum(w.orig)<1)&&verbose) 
    catf("Components of W*(V) do not sum to 1, dependent loss ('loss V') of ",
         "rescaled W*(V) is ",lossDep(Z,w.orig/sum(w.orig)),".\n\n")

  new.loss.w <- lossPred(X,w,tv)
  new.loss.v <- lossDep(Z,w)
  infeasible.w <- ((synth.out$loss.w-new.loss.w)/new.loss.w > tol) #&&
                    #(synth.out$loss.w-new.loss.w > tol)
  if (verbose) {
    cat0("Feasibility of W*(V)\n",
         "====================\n\n")
    if (infeasible.w) {
      catw("WARNING","W*(V) is NOT optimal and thus infeasible!")
      catw("'True' W*(V)",paste0(w,sep=" "))
      catf("with corresponding predictor loss ('loss W') of ",
           new.loss.w," ",
           "and corresponding dependent loss ('loss V') of ",
           new.loss.v,".\n\n")
    } else {       
      catw("GOOD","W*(V) is (essentially) optimal (new loss W: ",new.loss.w,
           ").\n\n")
    }       
  }     

  tmp    <- multiOpt(X,0,Z,0,outer.par=list(lb=lb),verbose=FALSE,
                     seed=seed,...)
  v      <- favourite_v(tmp$v,tmp$w,X,Z,tmp$trafo.v,lb=lb)
  v      <- v/sum(v)
  loss.w <- lossPred(X,tmp$w,v,tmp$trafo.v)
  loss.v <- tmp$rmspe^2
  suboptimal.v <- ((as.numeric(synth.out$loss.v)-loss.v)/loss.v > tol) #&&
                     #(as.numeric(synth.out$loss.v)-loss.v > tol)
                     
  if (verbose) {
    cat0("Optimality of V\n",
         "===============\n\n")
    
    if (suboptimal.v||infeasible.w) {
      catw("WARNING","'Optimal' V (as reported by package Synth) is not ",
           "optimal",if (infeasible.w) " (W*(V) was infeasible)",", ",
           "(one of potentially many) 'true' optimal V* (with sum(V*)=1):")
      catw("Optimal V*    ",paste0(v,sep=" "))
      catw("Optimal W*(V*)",paste0(tmp$w,sep=" "))
      catf("with corresponding predictor loss ('loss W') of ",
           loss.w," ","and corresponding dependent loss ('loss V') of ",
           loss.v,".\n")
    } else {       
      catw("GOOD","V is (essentially) optimal (new loss V: ",loss.v,").\n")     
    }  
  }    
  
  synth2.out                 <- synth.out
  if (is.numeric(synth2.out$solution.v)) {
    names(v) <- names(synth2.out$solution.v)
    synth2.out$solution.v  <- v 
  } else synth2.out$solution.v[1,]  <- v
  synth2.out$solution.w[,1]  <- tmp$w
  synth2.out$loss.w[1,1]     <- loss.w
  synth2.out$loss.v[1,1]     <- loss.v
  synth2.out$new.loss.v      <- new.loss.v
  synth2.out$new.loss.w      <- new.loss.w
  synth2.out$custom.v        <- NULL
  synth2.out$rgV.optim       <- NULL
       
  invisible(synth2.out)
}

checkInner <- function(X0,X1,v,w,Z0=NULL,Z1=NULL,tol=sqrt(.Machine$double.eps),
                       verbose=TRUE,...) {
  storage.mode(X0) <- storage.mode(X1) <- 
    storage.mode(Z0) <- storage.mode(Z1) <- "double"
  Xu  <- cbind(X0,X1)
  Xs  <- Xu/apply(Xu, 1, sd)
  X   <- Xs[, -ncol(Xs)] - drop(Xs[, ncol(Xs)])
  Z   <- if (!is.null(Z0)&&!is.null(Z1)) Z0 - drop(Z1) else NULL
  ME  <- 1L
  MA  <- length(v)
  N   <- ncol(X)
  MDW <- ME+MA
  globals$Ipar  <- as.integer(c(ME=ME,MA=MA,MDW=MDW,N=N)) 
  globals$IWORK <- integer(MDW+N)
  globals$WORK  <- double(MDW+5*N)
  globals$IWORK[1:2] <- c(length(globals$WORK),length(globals$IWORK))
  globals$RNORM <- double(1)
  globals$MODE  <- integer(1)
  globals$X     <- double(N)
  tv  <- as.numeric(v)
  sol <-.Fortran(C_wnnls,W=.Call(C_prepareW4,X,tv),
                 MDW=globals$Ipar[3],ME=globals$Ipar[1],MA=globals$Ipar[2],
                 N=globals$Ipar[4],L=0L,PRGOPT=1.0,X=globals$X,
                 RNORM=globals$RNORM,MODE=globals$MODE,IWORK=globals$IWORK,
                 WORK=globals$WORK)
  if ((any(is.infinite(sol$X))) || (sol$MODE>0)) 
    warning("error in inner optimization (wnnls)")                
  wnew <- sol$X
  new.loss.w <- lossPred(X,wnew,tv)
  old.loss.w <- lossPred(X,w,tv)
  new.loss.v <- if (is.null(Z)) NULL else lossDep(Z,wnew)
  old.loss.v <- if (is.null(Z)) NULL else lossDep(Z,w)
  abs.diff   <- old.loss.w-new.loss.w
  rel.diff   <- abs.diff/new.loss.w
  changed    <- (abs(rel.diff) > tol) #|| (abs(abs.diff) > tol)
  if (isTRUE(verbose&&changed))  
    cat0("Old: ",old.loss.w,", New: ",new.loss.w,", Relative Difference: ",
         (old.loss.w-new.loss.w)/new.loss.w,"\n")
  c(old.loss.w=old.loss.w, new.loss.w=new.loss.w, old.loss.v=old.loss.v,
    new.loss.v=new.loss.v)
}

rescaleLoss <- function(synth.out,dataprep.out) {
  storage.mode(dataprep.out$X0) <- storage.mode(dataprep.out$X1) <- 
    storage.mode(dataprep.out$Z0) <- storage.mode(dataprep.out$Z1) <- "double"
  Z   <- dataprep.out$Z0 - drop(dataprep.out$Z1)
  w.orig <- as.numeric(synth.out$solution.w)
  lossDep(Z,w.orig/sum(w.orig))
}

is.inner.wrong <- function(synth.out,dataprep.out,lb=1e-8,tol=1e-8) {
  storage.mode(dataprep.out$X0) <- storage.mode(dataprep.out$X1) <- 
    storage.mode(dataprep.out$Z0) <- storage.mode(dataprep.out$Z1) <- "double"
  v   <- as.numeric(synth.out$solution.v)
  Xu  <- cbind(dataprep.out$X0, dataprep.out$X1)                                # generate (scaled!) X from dataprep object
  Xs  <- Xu/apply(Xu, 1, sd)
  X   <- Xs[, -ncol(Xs)] - drop(Xs[, ncol(Xs)])
  Z   <- dataprep.out$Z0 - drop(dataprep.out$Z1)
  ME  <- 1L
  MA  <- length(v)
  N   <- ncol(X)
  MDW <- ME+MA
  globals$Ipar  <- as.integer(c(ME=ME,MA=MA,MDW=MDW,N=N)) 
  globals$IWORK <- integer(MDW+N)
  globals$WORK  <- double(MDW+5*N)
  globals$IWORK[1:2] <- c(length(globals$WORK),length(globals$IWORK))
  globals$RNORM <- double(1)
  globals$MODE  <- integer(1)
  globals$X     <- double(N)
  tv  <- as.numeric(v)
  sol <-.Fortran(C_wnnls,W=.Call(C_prepareW4,X,tv),
                 MDW=globals$Ipar[3],ME=globals$Ipar[1],MA=globals$Ipar[2],
                 N=globals$Ipar[4],L=0L,PRGOPT=1.0,X=globals$X,
                 RNORM=globals$RNORM,MODE=globals$MODE,IWORK=globals$IWORK,
                 WORK=globals$WORK)
  if ((any(is.infinite(sol$X))) || (sol$MODE>0)) 
    warning("error in inner optimization (wnnls)")                
  w <- sol$X
  w.orig <- as.numeric(synth.out$solution.w)

  real.loss.w <- lossPred(X,w.orig/sum(w.orig),tv)
  new.loss.w <- lossPred(X,w,tv)
  ((real.loss.w-new.loss.w)/new.loss.w > tol)
}

updateInner <- function(X0,X1,v,w,Z0=NULL,Z1=NULL,tol=sqrt(.Machine$double.eps),
                        verbose=TRUE,...) {
  storage.mode(X0) <- storage.mode(X1) <- 
    storage.mode(Z0) <- storage.mode(Z1) <- "double"
  Xu  <- cbind(X0,X1)
  Xs  <- Xu/apply(Xu, 1, sd)
  X   <- Xs[, -ncol(Xs)] - drop(Xs[, ncol(Xs)])
  Z   <- if (!is.null(Z0)&&!is.null(Z1)) Z0 - drop(Z1) else NULL
  ME  <- 1L
  MA  <- length(v)
  N   <- ncol(X)
  MDW <- ME+MA
  globals$Ipar  <- as.integer(c(ME=ME,MA=MA,MDW=MDW,N=N)) 
  globals$IWORK <- integer(MDW+N)
  globals$WORK  <- double(MDW+5*N)
  globals$IWORK[1:2] <- c(length(globals$WORK),length(globals$IWORK))
  globals$RNORM <- double(1)
  globals$MODE  <- integer(1)
  globals$X     <- double(N)
  tv  <- as.numeric(v)
  sol <-.Fortran(C_wnnls,W=.Call(C_prepareW4,X,tv),
                 MDW=globals$Ipar[3],ME=globals$Ipar[1],MA=globals$Ipar[2],
                 N=globals$Ipar[4],L=0L,PRGOPT=1.0,X=globals$X,
                 RNORM=globals$RNORM,MODE=globals$MODE,IWORK=globals$IWORK,
                 WORK=globals$WORK)
  if ((any(is.infinite(sol$X))) || (sol$MODE>0)) 
    warning("error in inner optimization (wnnls)")                
  wnew <- sol$X
  names(wnew) <- names(w)
  wnew
}


