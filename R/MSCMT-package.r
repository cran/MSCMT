#' Multivariate Synthetic Control Method Using Time Series
#'
#' MSCMT implements the Multivariate Synthetic Control Method Using
#' Time Series.
#'
#' MSCMT implements two generalizations of the synthetic control method (which 
#' has already an implementation in package 'Synth'): first, it allows for 
#' using multiple outcome variables, second, time series can be supplied as 
#' economic predictors. 
#' Much effort has been taken to make the implementation as stable as possible 
#' (including edge cases) without losing computational efficiency.
#'
#' @rdname MSCMTpackage
#' @docType package
#' @name MSCMT
#' @useDynLib MSCMT, .registration = TRUE, .fixes = "C_"
#' @examples
#' \dontrun{
#' ## for examples, see the package vignettes:
#' browseVignettes(package="MSCMT")
#' }
NULL
