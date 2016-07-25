## R interface to the LSEI Fortran algorithm 
## Fortran code from package limSolve, with kind permission of Karline Soetaert
## only a rudimentary interface, see package limSolve for a complete interface!
lsei <- function(A=NULL, B=NULL, E=NULL, F=NULL, G=NULL, H=NULL,
                 tol=sqrt(.Machine$double.eps))  {
  # fix (potentially incorrect) storage modes
  storage.mode(A) <- storage.mode(B) <- storage.mode(E) <- storage.mode(F) <- 
    storage.mode(G) <- storage.mode(H) <- "double"

  # calculate problem dimensions
  Neq  <- nrow(E)
  Napp <- nrow(A)
  Nx   <- ncol(A)
  Nin  <- nrow(G)
  ineq <- Nin+Nx
  mIP  <- ineq+2*Nx+2
  mdW  <- Neq + Napp + ineq
  mWS  <- 2*(Neq+Nx)+max(Napp+ineq,Nx)+(ineq+2)*(Nx+7)

  # call the optimizer
  sol <-.Fortran(C_lsei,NUnknowns=Nx,NEquations=Neq,NConstraints=Nin,
                 NApproximate=Napp,A=A,B=B,E=E,F=F,G=G,H=H,
                 X=numeric(Nx),mIP=as.integer(mIP),mdW=as.integer(mdW),
                 mWS=as.integer(mWS),IP=integer(mIP),W=numeric(mdW*(Nx+1)),
                 WS=numeric(mWS),lpr=1L,ProgOpt=as.double(1.0),
                 verbose=FALSE,IsError=FALSE)

  sol$X[which(abs(sol$X)<tol)] <- 0                                             # zero very tiny values
  list(x=sol$X)
}
