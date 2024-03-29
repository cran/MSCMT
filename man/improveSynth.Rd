% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/improveSynth.r
\name{improveSynth}
\alias{improveSynth}
\title{Check (and Improve) Results of Package Synth}
\usage{
improveSynth(
  synth.out,
  dataprep.out,
  lb = 1e-08,
  tol = 1e-05,
  verbose = TRUE,
  seed = 1,
  ...
)
}
\arguments{
\item{synth.out}{A result of \code{\link[Synth]{synth}} from package 
\code{'Synth'}.}

\item{dataprep.out}{The input of function \code{\link[Synth]{synth}} which
led to \code{synth.out}.}

\item{lb}{A numerical scalar (default: \code{1e-8}), corresponding to the
lower bound for the outer optimization.}

\item{tol}{A numerical scalar (default: \code{1e-5}). If the relative *and*
absolute improvement of \code{loss.v} and \code{loss.w}, respectively,
exceed \code{tol}, this is reported 
(if \code{verbose} is \code{TRUE}). Better optima are always looked for 
(independent of the value of \code{tol}).}

\item{verbose}{A logical scalar. Should the ouput be verbose (defaults to 
\code{TRUE}).}

\item{seed}{A numerical vector or \code{NULL}. See the corresponding
documentation for \code{\link[MSCMT]{mscmt}}. Defaults to 1 in order to
provide reproducibility of the results.}

\item{...}{Further arguments to \code{\link[MSCMT]{mscmt}}. Supported 
arguments are \code{check.global}, \code{inner.optim}, \code{inner.opar}, 
\code{outer.optim}, and \code{outer.opar}.}
}
\value{
An updated version of \code{synth.out}, where \code{solution.v},
\code{solution.w}, \code{loss.v}, and \code{loss.w} are replaced by the
optimum obtained by package \code{'MSCMT'} and all other components 
of \code{synth.out} are removed.
}
\description{
\code{improveSynth} checks the results of \code{\link[Synth]{synth}}
for feasibility and optimality and tries to find a better solution.
}
\details{
Performing SCM means solving a nested optimization problem. Depending on
the validity of the results of the inner optimization, SCM may produce
\itemize{
\item invalid or infeasible results, if the vector \code{w} of donor 
weights reported as the result of the inner optimization 
is in fact not optimal, ie. produces too large \code{loss.w},
\item suboptimal results, if the vector \code{v} of predictor weights 
reported as the result of the outer optimization is in fact not
optimal (which may be caused by shortcomings of the inner optimization).
}

\code{improveSynth} first checks \code{synth.out} for feasibility and
then tries to find a feasible and optimal solution by applying the
optimization methods of package \code{MSCMT} to \code{dataprep.out} 
(with default settings, more flexibility will probably be added in a 
future release).
}
\examples{
\dontrun{
## check whether package 'Synth' is available
if (require("Synth")) {

## process first example of function "synth" in package 'Synth' 
## (comments are removed):

  data(synth.data)
  dataprep.out<-
    dataprep(
     foo = synth.data,
     predictors = c("X1", "X2", "X3"),
     predictors.op = "mean",
     dependent = "Y",
     unit.variable = "unit.num",
     time.variable = "year",
     special.predictors = list(
        list("Y", 1991, "mean"),
        list("Y", 1985, "mean"),
        list("Y", 1980, "mean")
                              ),
     treatment.identifier = 7,
     controls.identifier = c(29, 2, 13, 17, 32, 38),
     time.predictors.prior = c(1984:1989),
     time.optimize.ssr = c(1984:1990),
     unit.names.variable = "name",
     time.plot = 1984:1996
     )

  synth.out <- synth(dataprep.out)

## check and (try to) improve these results:
  synth2.out <- improveSynth(synth.out,dataprep.out)
}
}
}
