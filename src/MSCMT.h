// From inverse.f
void F77_SUB(dwnnls)(double* W, int* MDW, int* ME, int* MA, int* N, int* L,
                     double* PRGOPT, double* X, double* RNORM, int* MODE, 
                     int* IWORK, double* WORK);
void F77_SUB(lsei)(int* NUnknowns, int* NEquations, int* NConstraints,
                   int* NApproximate, double* A, double* B, double* E, 
                   double* F, double* G, double* H, double* X, int* mIP, 
                   int* mdW, int* mWS, int* IP, double* W, double* WS, 
                   int* lpr, double* ProgOpt, int* verbose, int* IsError);

// From Helpers.c
SEXP fastVpMV(SEXP M, SEXP V);
SEXP fastVpMpMV(SEXP M, SEXP V);
SEXP fastMpdVM(SEXP M, SEXP V);
SEXP prepareW1(SEXP Z, SEXP X, SEXP W);
SEXP prepareW4(SEXP X, SEXP V);

// FROM DE.c
SEXP DE(SEXP X, SEXP Z, SEXP LenV, SEXP NP, SEXP NG, SEXP F, SEXP CR, SEXP Min, 
        SEXP Max, SEXP Minimpr, SEXP Waitgen, SEXP OptSep, SEXP CheckAmbiguity, 
        SEXP Width);
