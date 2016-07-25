#' Check (and Improve) Results of Package Synth
#'
#' \code{improveSynth} checks the results of \code{\link[Synth]{synth}}
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
#' @param synth.out A result of \code{\link[Synth]{synth}} from package 
#' \code{'Synth'}.
#' @param dataprep.out The input of function \code{\link[Synth]{synth}} which
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
#' @return An updated version of \code{synth.out}, where \code{solution.v},
#' \code{solution.w}, \code{loss.v}, and \code{loss.w} are replaced by the
#' optimum obtained by package \code{'MSCMT'} and all other components 
#' of \code{synth.out} are removed.
#' @export improveSynth
#' @examples 
#' \dontrun{
#' ## check whether package 'Synth' is available
#' if (require("Synth")) {
#'
#' ## process first example of function "synth" in package 'Synth' 
#' ## (comments are removed):
#'
#'   data(synth.data)
#'   dataprep.out<-
#'     dataprep(
#'      foo = synth.data,
#'      predictors = c("X1", "X2", "X3"),
#'      predictors.op = "mean",
#'      dependent = "Y",
#'      unit.variable = "unit.num",
#'      time.variable = "year",
#'      special.predictors = list(
#'         list("Y", 1991, "mean"),
#'         list("Y", 1985, "mean"),
#'         list("Y", 1980, "mean")
#'                               ),
#'      treatment.identifier = 7,
#'      controls.identifier = c(29, 2, 13, 17, 32, 38),
#'      time.predictors.prior = c(1984:1989),
#'      time.optimize.ssr = c(1984:1990),
#'      unit.names.variable = "name",
#'      time.plot = 1984:1996
#'      )
#' 
#'   synth.out <- synth(dataprep.out)
#'
#' ## check and (try to) improve these results:
#'   synth2.out <- improveSynth(synth.out,dataprep.out)
#' }
#' }
#' @importFrom stats sd
improveSynth <- function(synth.out,dataprep.out,lb=1e-8,tol=1e-5,verbose=TRUE) {
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

  new.loss.w <- lossPred(X,w,tv)
  new.loss.v <- lossDep(Z,w)
  infeasible.w <- ((synth.out$loss.w-new.loss.w)/new.loss.w > tol) &&
                    (synth.out$loss.w-new.loss.w > tol)
  if (verbose) {
    cat0("Feasibility of W*(V)\n",
         "====================\n\n")
    if (infeasible.w) {
      catw("WARNING","W*(V) is NOT optimal and thus infeasible!")
      catw("'True' W*(V)",paste0(synth.out$solution.w,sep=" "))
      catf("with corresponding predictor loss ('loss W') of ",
           new.loss.w," ",
           "and corresponding dependent loss ('loss V') of ",
           new.loss.v,".\n\n")
    } else {       
      catw("GOOD","W*(V) is (essentially) optimal (new loss W: ",new.loss.w,
           ").\n\n")
    }       
  }     

  tmp    <- modOptim(X,0,Z,0,outer.par=list(lb=lb),verbose=FALSE,
                     outer.optim="genoud")     
  v      <- favourite_v(tmp$w,X,tmp$trafo.v,lb=lb)
  v      <- v/sum(v)
  loss.w <- lossPred(X,tmp$w,v,tmp$trafo.v)
  loss.v <- tmp$rmspe^2
  suboptimal.v <- ((as.numeric(synth.out$loss.v)-loss.v)/loss.v > tol) &&
                     (as.numeric(synth.out$loss.v)-loss.v > tol)
                     
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
  synth2.out$solution.v[1,]  <- v
  synth2.out$solution.w[,1]  <- tmp$w
  synth2.out$loss.w[1,1]     <- loss.w
  synth2.out$loss.v[1,1]     <- loss.v
  synth2.out$custom.v        <- NULL
  synth2.out$rgV.optim       <- NULL
       
  invisible(synth2.out)
}
