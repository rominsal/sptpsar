##############################################################################
llik_reml_var2d <- function(x,sig2u,nsp,Wsp,y,X,Z,G_eff,np_eff,
                            bfixed=NULL,rho_fixed)
{
  rho <- x[1]
  In <- Matrix::Diagonal(nsp)
  if (!is.null(Wsp) & rho != 0) A <- In - rho*Wsp else A <- In
  V <- Matrix::tcrossprod(Z %*% Matrix::Matrix(diag(sqrt(G_eff)))) +
             sig2u*Matrix::Diagonal(nsp)
  Vinv <- Matrix::Matrix(lltAinv(as(V,"matrix")))
  if(class(Vinv) == "try-error") {
    Vinv <- Matrix::Matrix(MASS::ginv(as(V,"matrix")))
  }
  Vinv.X <- Vinv %*% X
  XtVinvX.inv <- Matrix::Matrix(lltAinv(as(Matrix::t(X) %*% Vinv.X,"matrix")))
  P <- try(Vinv - Vinv.X %*% (XtVinvX.inv %*% Matrix::t(Vinv.X)))
  if(class(P) == "try-error") {
    P <- Matrix::Matrix(Vinv - Vinv.X %*%
                          MASS::ginv(as(Matrix::t(X) %*% Vinv.X %*% Matrix::t(Vinv.X),
                                        "matrix")))
  }
  ldetV <- detlltA(as(V,"matrix"))$ldet
  ldet_XtVinvX <- detlltA(as(Matrix::t(X) %*% Vinv.X,"matrix"))$ldet
  #rm(chol.V,chol.XtVinvX)
  A_y <- Matrix::Matrix(A %*% y)
  Ay_P_Ay <- Matrix::t(A_y) %*% (P %*% A_y)
  ldetA <- Matrix::determinant(A,logarithm=TRUE)$modulus
  log_lik_reml <- as.numeric( -0.5*(ldetV + ldet_XtVinvX + Ay_P_Ay) + ldetA )
  log_lik <- as.numeric( -0.5*(ldetV + Ay_P_Ay) + ldetA )
  if (!rho_fixed ){
    # Compute the analytic hessian with respect to rho parameter REPASAR
    der2_reml_rho <- as.numeric( - Matrix::t(y) %*%
                                   (Matrix::t(Wsp) %*% (P %*% (Wsp  %*% y))) -
                                   sum(Matrix::diag(Matrix::solve(A,Wsp)^2)) )
    hessian_reml_rho <- c(der2_reml_rho)
    var_reml_rho <- 1/-hessian_reml_rho
    se_rho <- sqrt(var_reml_rho)
  } else se_rho <- 0

  res <- list(llik = log_lik,llik_reml = log_lik_reml,
              se_rho = se_rho)
  return(res)
}



##############################################################################
llik_reml_var3d <- function(x,sig2u,nsp,ntime,Wsp,y,X,Z,G_eff,np_eff,
                         bfixed=NULL,rho_fixed,phi_fixed)
{
    rho <- x[1]
    phi <- x[2]
    It <- Matrix::Diagonal(ntime)
    In <- Matrix::Diagonal(nsp)
    if (!is.null(Wsp) & rho != 0) A <- In - rho*Wsp else A <- In
    call_Omega <- build_Omega_ar1(phi,ntime)
    Omega <- Matrix::Matrix(call_Omega$Omega)
    der_sig2u_Omega_phi <- Matrix::Matrix(call_Omega$der_sig2u_Omega_phi)
    der2_sig2u_Omega_phi <- Matrix::Matrix(call_Omega$der2_sig2u_Omega_phi)
    rm(call_Omega)
    # Create a list with the replicates of Omega matrix
    # to build a blockdiag matrix
    Omega.rep <- plyr::alply(replicate(nsp,as(Omega,"matrix")),3)
    if (phi != 0) {
        #V <- tcrossprod(Z %*% diag(sqrt(G_eff))) +
        #        sig2u*Matrix::bdiag(Omega.rep)
        V <- Matrix::tcrossprod(Z %*% Matrix::Matrix(diag(sqrt(G_eff)))) +
                sig2u*Matrix::bdiag(Omega.rep)

    } else {
            #V <- tcrossprod(Z %*% diag(sqrt(G_eff))) +
            #    sig2u*Matrix::Diagonal(nsp*ntime)
            V <- Matrix::tcrossprod(Z %*% Matrix::Matrix(diag(sqrt(G_eff)))) +
                              sig2u*Matrix::Diagonal(nsp*ntime)

    }
    rm(Omega.rep)
    #chol.V <- chol(V)
    #chol.V <- lltA(V)
    #Vinv <- try(chol2inv(chol.V))
    Vinv <- Matrix::Matrix(lltAinv(as(V,"matrix")))
    if(class(Vinv) == "try-error") {
      Vinv <- Matrix::Matrix(MASS::ginv(as(V,"matrix")))
    }
    Vinv.X <- Vinv %*% X
    #chol.XtVinvX <- chol(t(X)%*%Vinv.X)
    #chol.XtVinvX <- lltA(t(X)%*%Vinv.X)
    #XtVinvX.inv <- chol2inv(chol.XtVinvX)
    XtVinvX.inv <- Matrix::Matrix(lltAinv(as(Matrix::t(X) %*% Vinv.X,"matrix")))
    P <- try(Vinv - Vinv.X %*% (XtVinvX.inv %*% Matrix::t(Vinv.X)))
    if(class(P) == "try-error") {
        P <- Matrix::Matrix(Vinv - Vinv.X %*%
              MASS::ginv(as(Matrix::t(X) %*% Vinv.X %*% Matrix::t(Vinv.X),
                             "matrix")))
    }
    #rm(V,Vinv,Vinv.X,XtVinvX.inv)
    # Compute log-lik y log-lik REML
    #ldet.V <- sum(2*log(diag(chol.V)))
    #ldet.XtVinvX <- sum(2*log(diag(chol.XtVinvX )))
    ldet.V <- detlltA(as(V,"matrix"))$ldet
    ldet.XtVinvX <- detlltA(as(Matrix::t(X) %*% Vinv.X,"matrix"))$ldet
    #rm(chol.V,chol.XtVinvX)
    A_It_y <- Matrix::Matrix(matrix( RH(A,RH(It,
                                array(y,dim=c(ncol(It),ncol(A))))),
                                    ncol=1) )
    # Es lo mismo que esta operación (lo he comprobado):
    #      A_It_y <- matrix(kronecker(A,I_t)%*%y,ncol=1)
    AIty.P.AIty <- Matrix::t(A_It_y) %*% (P %*% A_It_y)
    ldet.A <- Matrix::determinant(A,logarithm=TRUE)$modulus
    ldet.AIt <- ntime*ldet.A
    log_lik_reml <- as.numeric( -0.5*(ldet.V + ldet.XtVinvX + AIty.P.AIty) +
                            ldet.AIt )
    log_lik <- as.numeric( -0.5*(ldet.V + AIty.P.AIty) + ldet.AIt )
    # Compute the analytic hessian with respect to rho and phi parameters
    if (!phi_fixed) {
      sig2eps <- sig2u*(1 - phi^2)
      der_V_phi <- sig2eps*kronecker(In,der_sig2u_Omega_phi)
      der2.V_phi <- sig2eps*kronecker(In,der2_sig2u_Omega_phi)
      der_Vinv_phi <- -Vinv %*% der_V_phi %*% Vinv
      der_Xt.Vinv.X.inv <- -XtVinvX.inv %*%
        (Matrix::t(X) %*% der_Vinv_phi %*% X) %*% XtVinvX.inv
      der_Vinv.derVphi_Vinv_phi <- der_Vinv_phi %*% der_V_phi %*% Vinv +
        Vinv %*% der2.V_phi %*% Vinv + Vinv %*% der_V_phi %*% der_Vinv_phi

      der_P_phi <- der_Vinv_phi - ((der_Vinv_phi %*% X) %*%
                                     XtVinvX.inv %*% Matrix::t(Vinv.X) +
                                     Vinv.X %*% der_Xt.Vinv.X.inv %*% Matrix::t(Vinv.X) +
                                     Vinv.X %*% XtVinvX.inv %*%
                                     (Matrix::t(X) %*% der_Vinv_phi))
      der2.reml_phi <- as.numeric(
      -0.5*sum(Matrix::diag(der_P_phi %*% der_V_phi + P %*% der2.V_phi)) +
        0.5*Matrix::t((A %x% It) %*% (y - X %*% bfixed)) %*%
          (der_Vinv.derVphi_Vinv_phi %*%
             ((A %x% It) %*% (y - X %*% bfixed))) )
    } else { sigeps <- sig2u; der2.reml_phi <- 0}
    if (!rho_fixed) {
      der2.reml_rho <- as.numeric( -Matrix::t(y) %*%
                                (Matrix::t(Wsp) %x% It) %*% (P %*%
                                      ((Wsp %x% It) %*% y)) -
                          ntime*sum(Matrix::diag(Matrix::solve(A,Wsp)^2)) )
    } else der2.reml_rho <- 0
    if (!rho_fixed & !phi_fixed){
        der2.reml_phi_rho <- as.numeric(
               (Matrix::t(y) %*% (Matrix::t(Wsp) %x% It)) %*%
          ((Vinv %*% der_V_phi %*% Vinv) %*%
               ((A %x% It) %*% (y - X %*% bfixed))) )

    } else der2.reml_phi_rho <- 0
    if (!rho_fixed & !phi_fixed){
        hessian.reml_rho_phi <- matrix(c(der2.reml_rho,der2.reml_phi_rho,
                           der2.reml_phi_rho,der2.reml_phi),nrow=2,ncol=2)
        var_reml_rho_phi <- solve(-hessian.reml_rho_phi)
        se_rho <- sqrt(var_reml_rho_phi[1,1])
        se_phi <- sqrt(var_reml_rho_phi[2,2])
        cov_rho_phi <- var_reml_rho_phi[1,2] }
    if (!rho_fixed & phi_fixed){
        hessian.reml_rho_phi <- matrix(der2.reml_rho,nrow=1,ncol=1)
        var_reml_rho_phi <- solve(-hessian.reml_rho_phi)
        se_rho <- sqrt(var_reml_rho_phi[1,1])
        se_phi <- 0
        cov_rho_phi <- 0 }
    if (rho_fixed & !phi_fixed){
        hessian.reml_rho_phi <- matrix(der2.reml_phi,nrow=1,ncol=1)
        var_reml_rho_phi <- solve(-hessian.reml_rho_phi)
        se_rho <- 0
        se_phi <- sqrt(var_reml_rho_phi[1,1])
        cov_rho_phi <- 0 }
    if (rho_fixed & phi_fixed){se_rho <- 0; se_phi <- 0; cov_rho_phi <- 0}
    res <- list(llik = log_lik,llik_reml = log_lik_reml,
                se_phi = se_phi,se_rho = se_rho, cov_rho_phi = cov_rho_phi)
    return(res)
}


