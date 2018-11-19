llik_reml_fn2d <- function(x,sig2u,nsp=1,Wsp=NULL,y,X,Z,G_eff,np_eff,
                         bfixed=NULL,rho_fixed)
{
  rho <- x[1]
  ############ ASSUMPTION: nsp>0, otherwise FAILS ###############
  In <- Matrix::Diagonal(nsp)
  if (!is.null(Wsp) & rho != 0 ) {
    A <- Matrix::Diagonal(ncol(Wsp)) - rho*Wsp
  } else A <- In
  V <- Matrix::tcrossprod(Z %*% Matrix::Matrix(diag(sqrt(G_eff)))) +
           sig2u*Matrix::Diagonal(nsp)
  Vinv <- lltAinv(as(V,"matrix"))
  if(class(Vinv) == "try-error") {
    Vinv <- Matrix::Matrix(MASS::ginv(as(V,"matrix")))
  }
  Vinv.X <- Vinv %*% X
  XtVinvX.inv <- lltAinv(as(Matrix::t(X) %*% Vinv.X,"matrix"))
  P <- try(Vinv - Vinv.X %*% (XtVinvX.inv %*% Matrix::t(Vinv.X)))
  if(class(P) == "try-error") {
    P <- Vinv - Vinv.X %*%
      MASS::ginv(as(Matrix::t(X) %*% Vinv.X %*% Matrix::t(Vinv.X),
                    "matrix"))
  }
  # Compute log-lik REML
  ldetV <- detlltA(as(V,"matrix"))$ldet
  ldet_XtVinvX <- detlltA(as(Matrix::t(X) %*% Vinv.X,"matrix"))$ldet
  #rm(chol.V,chol.XtVinvX)
  A_y <- Matrix::Matrix(A %*% y)
  Ay_P_Ay <- Matrix::t(A_y) %*% (P %*% A_y)
  ldetA <- Matrix::determinant(A)$modulus
  log_lik_reml <- -0.5*(ldetV + ldet_XtVinvX + Ay_P_Ay) + ldetA
  return(as.numeric(-log_lik_reml))
}


#######################################################################################
llik_reml_fn3d <- function(x,sig2u,nsp=1,ntime=1,Wsp=NULL,y,X,Z,G_eff,np_eff,
                         bfixed=NULL,rho_fixed,phi_fixed)
{
    rho <- x[1]
    if (length(x) > 1) phi <- x[2] else phi <- 0
    ############ ASSUMPTION: nsp>0, otherwise FAILS ###############
    In <- Matrix::Diagonal(nsp)
    It <- Matrix::Diagonal(ntime)
    if (!is.null(Wsp) & rho != 0 ) {
      A <- Matrix::Diagonal(ncol(Wsp)) - rho*Wsp
    } else A <- In
    if (phi != 0) {
      call_Omega <- build_Omega_ar1(phi,ntime)
      Omega <- Matrix::Matrix(call_Omega$Omega)
      rm(call_Omega)
      # Create a list with the replicates of Omega matrix
      # to build a blockdiag matrix
      Omega.rep <- plyr::alply(replicate(nsp,as.matrix(Omega)),3)
      V <- Matrix::tcrossprod(Z %*% Matrix::Matrix(diag(sqrt(G_eff)))) +
                sig2u*Matrix::bdiag(Omega.rep)
      rm(Omega.rep)
    }  else V <- Matrix::tcrossprod(Z %*% Matrix::Matrix(diag(sqrt(G_eff)))) +
                     sig2u*Matrix::Diagonal(nsp*ntime)
    Vinv <- lltAinv(as(V,"matrix"))
    if(class(Vinv) == "try-error") {
      Vinv <- Matrix::Matrix(MASS::ginv(as(V,"matrix")))
    }
    Vinv.X <- Vinv %*% X
    XtVinvX.inv <- lltAinv(as(Matrix::t(X) %*% Vinv.X,"matrix"))
    P <- try(Vinv - Vinv.X %*% (XtVinvX.inv %*% Matrix::t(Vinv.X)))
    if(class(P) == "try-error") {
        P <- Vinv - Vinv.X %*%
                  MASS::ginv(as(Matrix::t(X) %*% Vinv.X %*% Matrix::t(Vinv.X),
                                "matrix"))
    }
    # Compute log-lik y log-lik REML
    ldet.V <- detlltA(as(V,"matrix"))$ldet
    ldet.XtVinvX <- detlltA(as(Matrix::t(X) %*% Vinv.X,"matrix"))$ldet
    #rm(chol.V,chol.XtVinvX)
    A_It_y <- Matrix::Matrix( matrix(RH(A,
                    RH(It, array(y,dim=c(ncol(It),ncol(A))))),ncol=1) )
    # Es lo mismo que esta operación (lo he comprobado):
    #      A_It_y <- matrix(kronecker(A,I_t)%*%y,ncol=1)
    AIty.P.AIty <- Matrix::t(A_It_y) %*% (P %*% A_It_y)
    ldet.A <- Matrix::determinant(A)$modulus
    ldet.AIt <- ntime*ldet.A
    log_lik_reml <- -0.5*(ldet.V + ldet.XtVinvX + AIty.P.AIty) + ldet.AIt
    #log_lik <- -0.5*(ldet.V + AIty.P.AIty) + ldet.AIt
    return(as.numeric(-log_lik_reml))
}


