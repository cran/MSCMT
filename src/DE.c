// Combination of outer optimization via Differential Evolution
// and inner optimization via dwnnls (-> big speed boost!)
//
// For the Differential Evolution algorithm, see, e.g.:
// Gilli, Manfred, Maringer, Dietmar and Schumann, Enrico:
// Numerical Methods and Optimization in Finance, Elsevier, 2011,
// algorithm 49, p. 347.
// 
// A (pure) R implementation of algorithm 49 is contained in
// function DEopt of package NMOF. Some of the variable names in
// the following (integrated) C implementation have been adopted 
// from DEopt to make both implementations (somehow) comparable.

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Applic.h>
#include <Rmath.h>

#include "MSCMT.h"

SEXP DE(SEXP X, SEXP Z, SEXP LenV, SEXP NP, SEXP NG, SEXP F, SEXP CR, SEXP Min, 
        SEXP Max, SEXP Minimpr, SEXP Waitgen, SEXP OptSep, SEXP CheckAmbiguity, SEXP Width) {

  int np   = 0;
  int best = 0;
  int i,j,k,l,m,vi,li;

  // make SEXP function arguments accessible
  const double * const x       = REAL(X);
  const double * const z       = REAL(Z);
  const int *    const lenV    = INTEGER(LenV);
  const int            ndonors = ncols(Z);
  const int            ndepend = nrows(Z);
  const int            nv      = nrows(X);
  const int            nP      = asInteger(NP); 
  const int            nG      = asInteger(NG); 
  const double         f       = asReal(F);
  const double         cr      = asReal(CR);
  const double         vmin    = asReal(Min);
  const double         vmax    = asReal(Max);
  const double         minimpr = asReal(Minimpr);
  const int            waitgen = asInteger(Waitgen);
  const int            optSep  = asLogical(OptSep);
  const int            checkA  = asLogical(CheckAmbiguity);
  const int            pwidth  = asInteger(Width);

  int                  d       = length(LenV);                                   // initialize dimension of outer optimization problem
  
  if (optSep) {                                                                  // do we fix one of the v[i] to vmax?
    li = d;                                                                      //   yes, so outer loop has length d and ... 
    d  = d - 1;                                                                  //     ... 'true' problem dimension is reduced by 1
  } else li = 1;                                                                 //   no, so outer loop has length 1

  // initialize SEXP objects for the (final) result
  SEXP Res, Names, Par, Value, Count, Nruns;
  Res   = PROTECT(allocVector(VECSXP,4));  np++;
  Names = PROTECT(allocVector(STRSXP,4));  np++;
  Nruns = PROTECT(allocVector(INTSXP,1));  np++;
  Par   = PROTECT(allocVector(REALSXP,d+optSep)); np++;
  Value = PROTECT(allocVector(REALSXP,1)); np++;
  Count = PROTECT(allocVector(INTSXP,1));  np++;
  SET_VECTOR_ELT(Res,0,Par);
  SET_VECTOR_ELT(Res,1,Value);
  SET_VECTOR_ELT(Res,2,Count);
  SET_VECTOR_ELT(Res,3,Nruns);
  SET_STRING_ELT(Names, 0, mkChar("par"));
  SET_STRING_ELT(Names, 1, mkChar("value"));
  SET_STRING_ELT(Names, 2, mkChar("counts"));
  SET_STRING_ELT(Names, 3, mkChar("nruns"));
  setAttrib(Res, R_NamesSymbol, Names);
  int *    nruns = INTEGER(Nruns);
  double * value = REAL(Value);
  int      nimpr = 0;
  int      impr  = 0;
  nruns[0] = 0;
  value[0] = R_PosInf;
  
  // prepare some variables and allocate memory for WNNLS
  double Rnorm   = 0.0;
  int    Mode    = 0;
  double PRGOPT  = 1.0;
  int    L       = 0;
  int    ME      = 1;
  int    MA      = nv;
  int    MDW     = ME + MA;
  int    N       = ndonors;
  double * Wmat  = (double *) R_alloc((ndonors+1)*(nv+1),sizeof(double));
  double * Work  = (double *) R_alloc(MDW+5*N, sizeof(double));
  double * W     = (double *) R_alloc(N, sizeof(double));
  int *    Iwork = (int *)    R_alloc(MDW+N ,sizeof(int));  
  
  // prepare some variables and allocate memory for WNNLS (for 'ambiguity check')
  int    ME2      = nv+1;
  int    MA2      = ndepend;
  int    MDW2     = ME2 + MA2;
  double * Wmat2  = (double *) R_alloc((ndonors+1)*(ndepend+nv+1),sizeof(double));
  double * Work2  = (double *) R_alloc(MDW2+5*N, sizeof(double));
  double * W2     = (double *) R_alloc(N, sizeof(double));
  int *    Iwork2 = (int *)    R_alloc(MDW2+N ,sizeof(int));  
  
  // allocate memory for DE algorithm
  double * v   = (double *) R_alloc(nv,sizeof(double));
  double * mv  = (double *) R_alloc(ndepend,sizeof(double));
  double * bV  = (double *) R_alloc(nG,sizeof(double));
  double * vF  = (double *) R_alloc(nP,sizeof(double));
  double * vFv = (double *) R_alloc(nP,sizeof(double));
  double * mP  = (double *) R_alloc(nP*d,sizeof(double));
  double * mPv = (double *) R_alloc(nP*d,sizeof(double));
  int *    vI  = (int *)    R_alloc(nP*d,sizeof(int));
  int *    sI  = (int *)    R_alloc(nP,sizeof(int));
  
  // prepare progress bar
  int    width;
  if (pwidth>=15) width = pwidth; else width = 0;                                // If we do not have enough space, cancel progress bar
  int    done    = 0;
  int    ip      = 0;
  int    tbd     = nG*li;
  int    perc    = -1;
  char * buffer  = (char *) R_alloc(width+13,sizeof(char));
  if (width) {
    buffer[0] = '\r';
    buffer[1] = '[';
    if (100*done/tbd>perc) {
      perc = 100*done/tbd;
      for (ip=2;ip<2+(width-7)*done/tbd;ip++) buffer[ip] = '=';
      for (;ip<width-5;ip++) buffer[ip] = '-';
      snprintf(buffer+ip,16,"] %3i%s",perc,"%%");
      REprintf(buffer);
    }  
  }
  
  GetRNGstate();	

  for (vi=0;vi<li;vi++) {
    // set up initial population
    j = 0;
    /*
    if ((d<15)&&((1<<d)<nP)) {                                                   // if nP is big enough, begin with a 'grid' ...
      for (k=0;k<(1<<d);k++) {                                                   // ... where each predictor weight is either vmin or vmax 
        m = k; 
        for (l=0;l<d;l++) {
          if (m&1) mP[j] = vmin; else mP[j] = vmax;
          m = m >> 1;
          j++;
        }  
      }
    }  
    */
    for (;j<nP*d;j++) mP[j] = vmin + unif_rand() * (vmax-vmin);                  // generate (remainder of) initial population at random
//    for (;j<nP*d;j++) mP[j] = vmin + pow(1.0-pow(unif_rand(),5.0),20.0) * (vmax-vmin);                  // generate (remainder of) initial population at random

    // evaluate initial population
    for (j=0;j<nP;j++) { 
      nruns[0] = nruns[0] + 1;                                                   // increase evaluation counter for benchmarks
      // create predictor weights
      m=0; 
      for (k=0;k<d+optSep;k++) 
        for (l=0;l<lenV[k];l++) 
          if (optSep&&(k==vi)) v[m++] = pow(10.0,vmax); else                     // difficult to read, but necessary to handle log-space
                               v[m++] = pow(10.0,mP[j*d+k-(optSep&&(k>vi))]);    // and potentially fixed (to 10^vmax) predictor weights.
      // prepare Wmat and Iwork for dwnnls and call dwnnls
      for (k=0;k<ndonors+1;k++) Wmat[k*(nv+1)] = 1.0;                            // row of 1's (for iota'w=1)
      for (k=0;k<nv;k++) {
        Wmat[k+1+ndonors*(nv+1)] = 0.0;                                          // append column of 0's
        for (m=0;m<ndonors;m++) Wmat[k+1+m*(nv+1)] = sqrt(v[k])*x[k+m*nv];       // build diag(sqrt(v)) %*% X
      }  
      Iwork[0] = MDW + 5*N; Iwork[1] = MDW + N;
      F77_CALL(dwnnls)(Wmat, &MDW, &ME, &MA, &N, &L, &PRGOPT, W, &Rnorm,         // solve inner optimization W = argmin_w w'X'diag(V)Xw  
                       &Mode, Iwork, Work);
      // check for ambiguity of (and possible 'improvements' for) W                 
      if (checkA) {
        for (k=0;k<nv;k++) {                                                     // prepare Wmat2
          for (m=0;m<ndonors;m++) Wmat2[k+m*(ndepend+nv+1)] = x[k+m*nv];         // copy X
          Wmat2[k+ndonors*(ndepend+nv+1)] = 0.0;                                 // append X'w
          for (m=0;m<ndonors;m++) 
            Wmat2[k+ndonors*(ndepend+nv+1)] += x[k+m*nv]*W[m];
        }
        for (m=0;m<ndonors+1;m++) Wmat2[nv+m*(ndepend+nv+1)] = 1.0;              // append row of 1's (for iota'w=1)
        for (k=0;k<ndepend;k++) {
          for (m=0;m<ndonors;m++)                                                // copy Z
            Wmat2[nv+1+k+m*(ndepend+nv+1)] = z[k+m*ndepend];   
          Wmat2[nv+1+k+ndonors*(ndepend+nv+1)] = 0.0;                            // append column of 0's
        }
        Iwork2[0] = MDW2 + 5*N; Iwork2[1] = MDW2 + N;
        F77_CALL(dwnnls)(Wmat2, &MDW2, &ME2, &MA2, &N, &L, &PRGOPT, W2,          // try to 'improve' solution by solving W2 = argmin_w* w*'Z'Zw* s.t. Xw*=Xw
                         &Rnorm, &Mode, Iwork2, Work2);
        impr = 0; 
        for (k=0;k<ndonors;k++)                                                  // are there any non-negligible differences between W and W2?
          if (fabs(W[k]-W2[k])>0.0000001) { impr = 1; break; }
        if (impr) {                                                              // yes, so set W = W2
          for (k=0;k<ndonors;k++) W[k] = W2[k]; 
          nimpr++;
        }  
      }                 
      // calculate dependent loss (W'Z'ZW)/ndepend
      vF[j] = 0.0;
      for (k=0;k<ndepend;k++) {
        mv[k] = 0.0;
        for (m=0;m<ndonors;m++) mv[k] += z[k+m*ndepend]*W[m];
      }  
      for (k=0;k<ndepend;k++) vF[j] += mv[k]*mv[k];
      vF[j] = vF[j] / ndepend;
    }

    // loop over nG (at max) generations
    for (i=0;i<nG;i++) {
      R_CheckUserInterrupt();                                                    // abort initiated by user?
      bV[i] = R_PosInf;                                                          // initialize best value of ith generation
      
      // update (and, if needed, repair) population
      l = nP;
      for (j=0;j<nP;j++) sI[j] = j;                                              // initialize sequence 1:nP ...
      for (j=0;j<nP;j++) {                                                       // ... and randomly permutate this sequence ...
        k     = (int)(l * unif_rand());                                          // ... by drawing nP times without replacement ...
        vI[j] = sI[k];                                                           // ... with result vI
        sI[k] = sI[--l];
      }

      for (j=0;j<nP;j++)                                                         // this is the core of differential evolution!
        for (k=0;k<d;k++)
          if (unif_rand() > cr)                                                  // no crossover, so no auxiliary solution necessary
            mPv[j*d+k] = mP[j*d+k]; 
          else {                                                                 // crossover, auxiliary solution calculated ...
            mPv[j*d+k] = mP[vI[(j+nP-1)%nP]*d+k]                                 // ... with modulo calculus for shifting the ... 
                + f * (mP[vI[(j+nP-2)%nP]*d+k] - mP[vI[(j+nP-3)%nP]*d+k]);       // ... random sequence vI
            if (mPv[j*d+k]<vmin) mPv[j*d+k]=vmin;                                // repair new population member (if necessary)
              else if (mPv[j*d+k]>vmax) mPv[j*d+k]=vmax;
          }	

      // evaluate updated population
      for (j=0;j<nP;j++) { 
        nruns[0] = nruns[0] + 1;                                                 // increase evaluation counter for benchmarks
        // create predictor weights
        m=0; 
        for (k=0;k<d+optSep;k++) 
          for (l=0;l<lenV[k];l++) 
            if (optSep&&(k==vi)) v[m++] = pow(10.0,vmax); else                   // difficult to read, but necessary to handle log-space
                                 v[m++] = pow(10.0,mPv[j*d+k-(optSep&&(k>vi))]); // and potentially fixed (to 10^vmax) predictor weights.
        // prepare Wmat and Iwork for dwnnls and call dwnnls
        for (k=0;k<ndonors+1;k++) Wmat[k*(nv+1)] = 1.0;                          // row of 1's (for iota'w=1)
        for (k=0;k<nv;k++) {
          Wmat[k+1+ndonors*(nv+1)] = 0.0;                                        // append column of 0's
          for (m=0;m<ndonors;m++) Wmat[k+1+m*(nv+1)] = sqrt(v[k])*x[k+m*nv];     // build diag(sqrt(v)) %*% X
        }  
        Iwork[0] = MDW + 5*N; Iwork[1] = MDW + N;
        F77_CALL(dwnnls)(Wmat, &MDW, &ME, &MA, &N, &L, &PRGOPT, W, &Rnorm,       // solve inner optimization W = argmin_w w'X'diag(V)Xw 
                         &Mode, Iwork, Work);
        // check for ambiguity of (and possible 'improvements' for) W                 
        if (checkA) {
          for (k=0;k<nv;k++) {                                                   // prepare Wmat2
            for (m=0;m<ndonors;m++) Wmat2[k+m*(ndepend+nv+1)] = x[k+m*nv];       // copy X
            Wmat2[k+ndonors*(ndepend+nv+1)] = 0.0;                               // append X'w
            for (m=0;m<ndonors;m++) 
              Wmat2[k+ndonors*(ndepend+nv+1)] += x[k+m*nv]*W[m];
          }
          for (m=0;m<ndonors+1;m++) Wmat2[nv+m*(ndepend+nv+1)] = 1.0;            // append row of 1's (for iota'w=1)
          for (k=0;k<ndepend;k++) {
            for (m=0;m<ndonors;m++)                                              // copy Z
              Wmat2[nv+1+k+m*(ndepend+nv+1)] = z[k+m*ndepend];   
            Wmat2[nv+1+k+ndonors*(ndepend+nv+1)] = 0.0;                          // append column of 0's
          }
          Iwork2[0] = MDW2 + 5*N; Iwork2[1] = MDW2 + N;
          F77_CALL(dwnnls)(Wmat2, &MDW2, &ME2, &MA2, &N, &L, &PRGOPT, W2,        // try to 'improve' solution by solving W2 = argmin_w* w*'Z'Zw* s.t. Xw*=Xw
                           &Rnorm, &Mode, Iwork2, Work2);
          impr = 0; 
          for (k=0;k<ndonors;k++)                                                // are there any non-negligible differences between W and W2?
            if (fabs(W[k]-W2[k])>0.0000001) { impr = 1; break; }
          if (impr) {                                                            // yes, so set W = W2
            for (k=0;k<ndonors;k++) W[k] = W2[k]; 
            nimpr++;
          }  
        }                 
        // calculate dependent loss (W'Z'ZW)/ndepend
        vFv[j] = 0.0;
        for (k=0;k<ndepend;k++) {
          mv[k] = 0.0;
          for (m=0;m<ndonors;m++) mv[k] += z[k+m*ndepend]*W[m];
        }  
        for (k=0;k<ndepend;k++) vFv[j] += mv[k]*mv[k];
        vFv[j] = vFv[j] / ndepend;
      }
      // update population members with improvements and look for best member
      for (j=0;j<nP;j++) {
        if(vFv[j]<vF[j]) {
          vF[j] = vFv[j];
          for (k=0;k<d;k++) mP[j*d+k] = mPv[j*d+k];
        }
        if (vF[j]<bV[i]) {
          best  = j;
          bV[i] = vF[j];
        }
      }
      // update progress bar
      if (width) {
        done++;
        if (100*done/tbd>perc) {
          perc = 100*done/tbd;
          for (ip=2;ip<2+(width-7)*done/tbd;ip++) buffer[ip] = '=';
          for (;ip<width-5;ip++) buffer[ip] = '-';
          snprintf(buffer+ip,16,"] %3i%s",perc,"%%");
          REprintf(buffer);
        }  
      }
      if (i>=waitgen) if (bV[i-waitgen]<(1.0+minimpr)*bV[i]) {                   // can we stop early because we did not improve?
        if (width) done += nG-i-1;                                               // yes, so let progress bar make a jump ...
        break;                                                                   // ... and leave loop over nP generations
      }
    } 
    if (vF[best]<value[0]) {                                                     // needed for fixed v[vi]'s. Is this (by now) the 'best' vi?
      REAL(Value)[0]    = vF[best];                                              // yes, so update (final) results of outer optimization ... 
      for (k=0;k<d+optSep;k++) 
        if (optSep&&(k==vi)) REAL(Par)[k] = vmax; else                           // ... respecting (potentially) fixed v[vi]'s ... 
                             REAL(Par)[k] = mP[best*d+k-(optSep&&(k>vi))];
      INTEGER(Count)[0] = i;                                                     // ... and save the number of generations for benchmarks
    }
  }
  // finish progress bar and report number of 'ambiguities' (if applicable)
  if (width) {
    perc = 100;
    for (ip=2;ip<width-5;ip++) buffer[ip] = '=';
    snprintf(buffer+ip,16,"] %3i%s",perc,"%%");
    REprintf(buffer);
    REprintf("\n");
    if (checkA) REprintf("Number of observed improvements (ambiguities): %d\n",nimpr);
  }
  
  // tidy up and return!	
  PutRNGstate();
  UNPROTECT(np);
  return(Res);
}
