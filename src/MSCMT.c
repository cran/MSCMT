#include <R_ext/Rdynload.h>
#include <R_ext/Applic.h>
#include <R_ext/Callbacks.h>

#include "MSCMT.h"

// Register native routines
static R_CallMethodDef callMethods[] = {
  {"fastVpMV", (DL_FUNC) &fastVpMV, 2},
  {"fastVpMpMV", (DL_FUNC) &fastVpMpMV, 2},
  {"fastMpdVM", (DL_FUNC) &fastMpdVM, 2},
  {"prepareW1", (DL_FUNC) &prepareW1, 3},
  {"prepareW4", (DL_FUNC) &prepareW4, 2},
  {NULL, NULL, 0}
};

static R_NativePrimitiveArgType DWNNLS_t[] = {
  REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, 
  REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, REALSXP
};

void F77_SUB(dwnnls)(double* W, int* MDW, int* ME, int* MA, int* N, int* L,
  double* PRGOPT, double* X, double* RNORM, int* MODE, int* IWORK, double* WORK);

static R_NativePrimitiveArgType lsei_t[] = {
  INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, 
  REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, INTSXP,
  REALSXP, LGLSXP, LGLSXP
};

void F77_SUB(lsei)(int* NUnknowns, int* NEquations, int* NConstraints,
  int* NApproximate, double* A, double* B, double* E, double* F, double* G,
  double* H, double* X, int* mIP, int* mdW, int* mWS, int* IP, double* W,
  double* WS, int* lpr, double* ProgOpt, int* verbose, int* IsError);
  
static R_FortranMethodDef fortranMethods[] = {
  {"wnnls", (DL_FUNC) &dwnnls_, 12, DWNNLS_t},
  {"lsei",  (DL_FUNC) &lsei_, 21, lsei_t},
  {NULL, NULL, 0}
};

void R_init_MSCMT(DllInfo *info) {
  R_registerRoutines(info, NULL, callMethods, fortranMethods, NULL);
}
