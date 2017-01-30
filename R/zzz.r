#' @importFrom utils packageVersion
.onAttach <- function(libname, pkgname) {
  version <- utils::packageVersion(pkgname)
  if(isTRUE(interactive())) packageStartupMessage(paste0(
    "# MSCMT version ",version,"\n#\n",
    "# When publishing results obtained using this package,\n",
    "# the original authors are to be cited as:\n",
    "# Martin Becker and Stefan Kloessner (2017). \n",
    "# MSCMT: Multivariate Synthetic Control Method Using Time Series. \n",
    "# R package version ",version,"."
  ))
}

.onUnload <- function(libpath) {
  library.dynam.unload(utils::packageName(), libpath)
}
