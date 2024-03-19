#' Multivariate Synthetic Control Method Using Time Series
#'
#' MSCMT implements the Multivariate Synthetic Control Method Using
#' Time Series.
#'
#' MSCMT implements three generalizations of the synthetic control method (which 
#' has already an implementation in package 'Synth'): \enumerate{
#' \item it allows for using multiple outcome variables, 
#' \item time series can be supplied as economic predictors,
#' \item a well-defined cross-validation approach can be used.
#' } 
#' Much effort has been taken to make the implementation as stable as possible 
#' (including edge cases) without losing computational efficiency.
#'
#' @rdname MSCMTpackage
#' @aliases MSCMT-package
#' @keywords internal 
#' @name MSCMT
#' @useDynLib MSCMT, .registration = TRUE, .fixes = "C_"
#' @importFrom Rdpack reprompt
#' @importFrom Rdpack insert_ref
#' @examples
#' \dontrun{
#' ## for examples, see the package vignettes:
#' browseVignettes(package="MSCMT")
#' }
#' @references
#' \insertRef{Abadie2003}{MSCMT}
#'
#' \insertRef{Abadie2010}{MSCMT}
#'
#' \insertRef{FastReliable}{MSCMT}
#'
#' \insertRef{CV}{MSCMT}
#'
#' \insertRef{KP16}{MSCMT}
#'
"_PACKAGE"
