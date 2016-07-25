#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Applic.h>
#include <Rmath.h>

// calculates V' M V for a vector V and a corresponding matrix M
SEXP fastVpMV(SEXP M, SEXP V) {
	const double * const m = REAL(M);
	const double * const v = REAL(V);
  const int            n = length(V);

	SEXP Res;
	PROTECT (Res = allocVector(REALSXP, 1));
	double * const res = REAL(Res);
	double * const mv  = (double *) R_alloc(n,sizeof(double));
    
  int i,j;
  res[0] = 0.0;
  for (i=0;i<n;i++) mv[i] = 0.0;
  for (i=0;i<n;i++) for (j=0;j<n;j++) mv[i] += m[i+j*n]*v[j];
  for (i=0;i<n;i++) res[0] += v[i]*mv[i];

  UNPROTECT(1);
  return(Res);
}

// calculates V' M' M V for a vector V and a corresponding matrix M
SEXP fastVpMpMV(SEXP M, SEXP V) {
	const double * const m = REAL(M);
	const double * const v = REAL(V);
  const int            n = ncols(M);
  const int            l = nrows(M);

	SEXP Res;
	PROTECT (Res = allocVector(REALSXP, 1));
	double * const res = REAL(Res);
	double * const mv  = (double *) R_alloc(l,sizeof(double));
    
  int i,j;
  res[0] = 0.0;
  for (i=0;i<l;i++) mv[i] = 0.0;
  for (i=0;i<l;i++) for (j=0;j<n;j++) mv[i] += m[i+j*l]*v[j];
  for (i=0;i<l;i++) res[0] += mv[i]*mv[i];

  UNPROTECT(1);
  return(Res);
}

// calculates M' diag(V) M for a matrix M and a corresponding vector V
SEXP fastMpdVM(SEXP M, SEXP V) {
	const double * const m = REAL(M);
	const double * const v = REAL(V);
  const int            n = ncols(M);
  const int            l = nrows(M);

	SEXP Res;
	PROTECT (Res = allocMatrix(REALSXP, n, n));
	double * const res = REAL(Res);
    
  int i,j,k;
  res[0] = 0.0;
  for (i=0;i<n;i++) {
    for (j=0;j<=i;j++) {
      res[i+j*n] = 0.0;
      for (k=0;k<l;k++) res[i+j*n] += v[k]*m[k+i*l]*m[k+j*l];
      if (i!=j) res[j+i*n] = res[i+j*n];
    }
  }

  UNPROTECT(1);
  return(Res);
}

// prepares Matrix W for wnnls alrogithm in improve_Zw
SEXP prepareW1(SEXP Z, SEXP X, SEXP W) {
	const double * const z = REAL(Z);
	const double * const x = REAL(X);
	const double * const w = REAL(W);
  const int            l = nrows(Z);
  const int            m = nrows(X);
  const int            n = ncols(Z);

  SEXP Res;
	PROTECT (Res = allocVector(REALSXP, (l+m+1)*(n+1)));
  double * const res = REAL(Res);
  int i,j;
    
  for (i=0;i<m;i++) {
    for (j=0;j<n;j++) res[i+j*(l+m+1)] = x[i+j*m];       // copy X
    res[i+n*(l+m+1)] = 0.0;                              // append X'w
    for (j=0;j<n;j++) res[i+n*(l+m+1)] += x[i+j*m]*w[j];
  }
  for (j=0;j<n+1;j++) res[m+j*(l+m+1)] = 1.0;            // append row of 1's (for iota'w=1)
  for (i=0;i<l;i++) {
    for (j=0;j<n;j++) res[m+1+i+j*(l+m+1)] = z[i+j*l];   // copy Z
    res[m+1+i+n*(l+m+1)] = 0.0;                          // append column of 0's
  }

  UNPROTECT(1);
  return(Res);
}

// prepares Matrix W for wnnls algorithm in opt.wnnls
SEXP prepareW4(SEXP X, SEXP V) {
	const double * const x = REAL(X);
	const double * const v = REAL(V);
  const int            m = nrows(X);
  const int            n = ncols(X);

  SEXP Res;
	PROTECT (Res = allocVector(REALSXP, (m+1)*(n+1)));
  double * const res = REAL(Res);
  int i,j;
    
  for (j=0;j<n+1;j++) res[j*(m+1)] = 1.0;                     // row of 1's (for iota'w=1)
  for (i=0;i<m;i++) {
    for (j=0;j<n;j++) res[i+1+j*(m+1)] = sqrt(v[i])*x[i+j*m]; // build sqrt(v)*X
    res[i+1+n*(m+1)] = 0.0;                                   // append column of 0's
  }

  UNPROTECT(1);
  return(Res);
}
