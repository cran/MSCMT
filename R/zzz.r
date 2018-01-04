#' @importFrom utils packageVersion
.onAttach <- function(libname, pkgname) {
  version <- utils::packageVersion(pkgname)
  if(isTRUE(interactive())) packageStartupMessage(paste0(
    "# MSCMT version ",version,"\n#\n",
    "# When publishing results obtained using this package,\n",
    "# the original authors are to be cited as:\n",
    "# Martin Becker and Stefan Kl\u00F6\u00DFner (2018). \n",
    "# Fast and Reliable Computation of Generalized Synthetic Controls. \n",
    "# Econometrics and Statistics, 5, 1-19."
  ))
}

.onUnload <- function(libpath) {
  library.dynam.unload(utils::packageName(), libpath)
}
