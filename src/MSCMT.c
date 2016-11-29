#include <R_ext/Rdynload.h>
#include <R_ext/Applic.h>
#include <R_ext/Callbacks.h>

#include "MSCMT.h"

// Prepare registering native C routines
static R_CallMethodDef callMethods[] = {
  {"fastVpMV", (DL_FUNC) &fastVpMV, 2},
  {"fastVpMpMV", (DL_FUNC) &fastVpMpMV, 2},
  {"fastMpdVM", (DL_FUNC) &fastMpdVM, 2},
  {"prepareW1", (DL_FUNC) &prepareW1, 3},
  {"prepareW4", (DL_FUNC) &prepareW4, 2},
  {"DE", (DL_FUNC) &DE, 14},
  {NULL, NULL, 0}
};

// Prepare registering native Fortran routines
static R_NativePrimitiveArgType DWNNLS_t[] = {
  REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, 
  REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, REALSXP
};

static R_NativePrimitiveArgType lsei_t[] = {
  INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, 
  REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, INTSXP,
  REALSXP, LGLSXP, LGLSXP
};

static R_FortranMethodDef fortranMethods[] = {
  {"wnnls", (DL_FUNC) &dwnnls_, 12, DWNNLS_t},
  {"lsei",  (DL_FUNC) &lsei_, 21, lsei_t},
  {NULL, NULL, 0}
};

// Register native C and Fortran routines
void R_init_MSCMT(DllInfo *info) {
  R_registerRoutines(info, NULL, callMethods, fortranMethods, NULL);
}
