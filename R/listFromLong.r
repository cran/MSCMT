#' Convert Long Format to List Format
#'
#' \code{listFromLong} converts long to list format.
#'
#' \code{listFromLong} is a convenience function to convert long format
#' (in a \code{\link[base]{data.frame}}, as used by package 'Synth') to list 
#' format, where data is stored as a list of matrices.
#'
#' Most parameter names are named after their equivalents in the 
#' \code{\link[Synth]{dataprep}} function of package 'Synth'.
#'
#' @param foo A \code{data.frame} containing the data in "long" format.
#' @param unit.variable Either a numeric scalar with the column number (in 
#' \code{foo}) containing the units or a character scalar with the corresponding
#' column name in \code{foo}.
#' @param time.variable Either a numeric scalar with the column number (in 
#' \code{foo}) containing the times or a character scalar with the corresponding
#' column name in \code{foo}.
#' @param unit.names.variable Optional. If not \code{NULL}, either a numeric 
#' scalar with the column number (in \code{foo}) containing the unit names or 
#' a character scalar with  the corresponding column name in \code{foo}. Must
#' match with the units defined by \code{unit.variable} (if not \code{NULL}).
#' @param exclude.columns Optional (defaults to \code{NULL}). Numeric vector 
#' with column numbers of \code{foo} to be excluded from the conversion.
#' @return A list of matrices with rows corresponding to the times and columns
#' corresponding to the unit (or unit names, respectively) for all columns of
#' \code{foo} which are neither excluded nor have a special role as time, unit,
#' or unit names variable.
#' @importFrom stats na.omit
#' @export listFromLong
#' @examples
#' if (require("Synth")) {
#'   data(basque)
#'   Basque <- listFromLong(basque, unit.variable="regionno", 
#'                          time.variable="year", 
#'                          unit.names.variable="regionname")
#'   names(Basque)
#'   head(Basque$gdpcap)
#' }
listFromLong <- function(foo, unit.variable, time.variable, 
                         unit.names.variable=NULL,exclude.columns=NULL) {
  if(!is.data.frame(foo)) stop("foo must be a data.frame")

  # main helper function
  DFtoList <- function(input,rowcol,colcol,colnamecol=NULL,exclude=NULL) {
    stopifnot(length(dim(input))==2)
    datcols    <- setdiff(seq_len(ncol(input)), 
                          c(rowcol,colcol,colnamecol,exclude))
    res        <- vector("list",length(datcols))
    names(res) <- if (!is.null(colnames(input))) colnames(input)[datcols] else 
                                                 as.character(datcols)
    if (!is.null(colnamecol)) {
      c2n        <- na.omit(unique(input[,colnamecol]))
      names(c2n) <- na.omit(unique(input[,colcol]))
    }  
                                                 
    for (i in seq_along(res)) {
      idx  <- !is.na(input[,datcols[i]])
      rown <- unique(input[idx,rowcol])
      coln <- unique(input[idx,colcol])
      res[[i]] <- matrix(NA,nrow=length(rown),ncol=length(coln))
      rownames(res[[i]]) <- rown
      colnames(res[[i]]) <- coln
      for (j in which(idx)) 
        res[[i]][as.character(input[j,rowcol]),as.character(input[j,colcol])] <- 
          input[j,datcols[i]]
      if (!is.null(colnamecol)) colnames(res[[i]]) <- c2n[as.character(coln)]
      res[[i]] <- res[[i]][order(rownames(res[[i]])),,drop=FALSE]
    }
    res
  }
    
  # helper function for data check and column lookup
  nam2idx <- function(...,id,type="numeric") {
    if (missing(id)) id <- as.character(match.call())[2]
    obj <- list(`...`)[[1]]
    if (length(obj)!=1) stop(paste0("\n Please specify only one ",id,"\n"))
    if (!(mode(foo[,obj]) %in% type)) stop("\n ",id," not found as ",type,
                                           " variable in foo.\n")
    if (mode(obj) == "character") which(names(foo)==obj) else obj
  }

  # check and look for unit.variable and time.variable
  unit.variable <- nam2idx(unit.variable)
  time.variable <- nam2idx(time.variable,type=c("numeric","character"))

  # check and look for unit.name.variable (if present)
  if (!is.null(unit.names.variable)) {
    idx <- !(is.na(foo[,unit.variable])|is.na(foo[,unit.names.variable]))
    unit.names.variable <- nam2idx(unit.names.variable,type="character")
    if (length(unique(foo[idx,unit.names.variable])) !=
          length(unique(foo[idx,unit.variable])))
        stop("lengths of unit.names and unit.names.variable do not match")
    if (length(unique(paste(foo[idx,unit.variable],foo[idx,unit.names.variable],
                            sep="----------"))) != 
          length(unique(foo[idx,unit.variable])))
        stop("unit.names and unit.names.variable do not match")
  }
  
  # do the actual conversion
  DFtoList(foo,rowcol=time.variable,colcol=unit.variable,
           colnamecol=unit.names.variable,exclude=exclude.columns)
}
