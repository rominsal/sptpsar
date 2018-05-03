# Functions to compute scores with respect to rho and phi.

score_rho <- function(param,sig2u,nsp,ntime,Wsp,y,X,Z,G.eff,np.eff,
                      b_fixed=NULL,trace=FALSE)
{
  rho <- param[1]/(1+abs(param[1])) #Ensures the interval (-1,1)
  #Wsp <- as(Wsp,"CsparseMatrix") # OJO: NO HE ELEGIDO SIMÉTRICA PORQUE
  # Wsp NO ES SIMÉTRICA SI ESTÁ ESTANDARIZADA POR COLUMNAS.
  A <- Matrix::Diagonal(nsp)-rho*Wsp
  It <- Matrix::Diagonal(ntime)
  A_It_y <- matrix(RH(A,RH(It, matrix(y,nrow=ntime))),ncol=1)
  # Es lo mismo que (comprobado): A_It_y <- matrix(kronecker(A,It)%*%y,ncol=1)
  W_It_y <-matrix(RH(Wsp,RH(It, matrix(y,nrow=ntime))),ncol=1)
  # Es lo mismo que (comprobado): W_It_y <- matrix(kronecker(Wsp,It)%*%y,ncol=1)
  # OJO: V NO TIENE PORQUÉ SER SPARSE
  #V <- as(tcrossprod(Z%*%diag(sqrt(G.eff))) +
  #            sig2u*(Matrix::Diagonal(nsp*ntime)),"dsyMatrix")
  V <- tcrossprod(Z%*%diag(sqrt(G.eff))) + sig2u*diag(nsp*ntime)
  #V <- AAt(Z%*%diag(sqrt(G.eff))) + sig2u*diag(nsp*ntime)
  #chol.V <- chol(V)
  #chol.V <- lltA(V)
  #Vinv <- try(chol2inv(chol.V))
  Vinv <- lltAinv(V)
  if(class(Vinv) == "try-error")  Vinv <- MASS::ginv(as(V,"matrix"))
  Vinv.X <- Vinv %*% X
  #chol.XtVinvX <- chol(t(X)%*%Vinv.X)
  #chol.XtVinvX <- lltA(t(X)%*%Vinv.X)
  #XtVinvX.inv <- chol2inv(chol.XtVinvX)
  XtVinvX.inv <- lltAinv(t(X)%*%Vinv.X)
  P <- try(Vinv - Vinv.X %*% (XtVinvX.inv %*% t(Vinv.X)))
  if(class(P) == "try-error") {
      P <- Vinv - Vinv.X %*% (MASS::ginv(as(t(X)%*%Vinv.X,"matrix"))
                              %*% t(Vinv.X))
  }
  # if(is.null(weights)) w <- as.vector(matrix(1,nrow=length(y))) else w <- weights
  # mat <- construct_matrices(X,Z,A_It_y,w,GLAM)
  # # ES LA MATRIZ C DE COEFICIENTES DEL SISTEMA (12) EN PAPER SAP
  # # NO SE PUEDE RESOLVER SI SE INCLUYEN LOS ZEROS EN MATRICES X,Z,G,
  # # POR ESO SE USA G.EFF
  # C <- construct_block(mat$XtX, t(mat$ZtX*G.eff), mat$ZtX,
  #                      t(mat$ZtZ*G.eff))
  # D <- diag(c(rep(0,np.eff[1]), rep(1,sum(np.eff[-1]))))
  # Hinv <- try(solve((1/sig2u)*C + D))
  #  if(class(Hinv) == "try-error") Hinv <- MASS::ginv((1/sig2u)*C + D)
  # P.Hv <- (1/sig2u)*(diag(length(y))-
  #                      (1/sig2u)*(cbind(X,Z%*%diag(G.eff)))%*%
  #                      (Hinv%*%t(cbind(X,Z))))
  score_rho <-  t(P%*%A_It_y) %*% W_It_y - ntime*sum(diag(solve(A,Wsp)))
  if (trace) cat("rho ",rho," score_rho ",as.numeric(score_rho),"\n")
  return(as.numeric(score_rho))
}
###############################################################################
score_phi <- function(param,sig2u,nsp,ntime,Wsp=NULL,y,X,Z,G.eff,
                      np.eff,b_fixed,trace=FALSE)
{
  phi <- param[1]/(1+abs(param[1])) #Ensures the interval (-1,1)
  call.Omega <- build_Omega_ar1(phi,ntime)
  Omega <- call.Omega$Omega
  der_sig2u_Omega_phi <- call.Omega$der_sig2u_Omega_phi
  #Omegainv <- call.Omega$Omegainv
  rm(call.Omega)
  Omega.rep <- plyr::alply(replicate(nsp,Omega),3)
  V <- tcrossprod(Z %*% diag(sqrt(G.eff))) + sig2u*Matrix::bdiag(Omega.rep)
  rm(Omega.rep)
  if(!is.matrix(V)) V <- as.matrix(V)
  Vinv <- lltAinv(V)
  if(class(Vinv) == "try-error")  Vinv <- MASS::ginv(as(V,"matrix"))
  Vinv.X <- Vinv %*% X
  XtVinvX.inv <- lltAinv(t(X)%*%Vinv.X)
  P <- try(Vinv - Vinv.X %*% (XtVinvX.inv %*% t(Vinv.X)))
  if(class(P) == "try-error") {
      P <- Vinv - Vinv.X %*% (MASS::ginv(as(t(X)%*%Vinv.X,"matrix"))
                              %*% t(Vinv.X))
  }
  A_It_y <- y
  Vinv_resids <- Vinv%*%(A_It_y-X%*%b_fixed)
  rm(V,Vinv,Vinv.X,XtVinvX.inv)
  #der_V_phi <- sig2u*kronecker(Matrix::Diagonal(nsp),der_sig2u_Omega_phi)
  #score_phi <- -0.5*sum(diag(as.matrix(P %*% der_V_phi))) +
  #    0.5*t(as.matrix(Vinv_resids)) %*% der_V_phi %*% Vinv_resids

  Insp <- Matrix::Diagonal(nsp)
  prod_mat <- apply(aperm(RH(Insp,RH(der_sig2u_Omega_phi,
                array(P,dim=c(ncol(der_sig2u_Omega_phi),ncol(Insp),ncol(P))))),
                perm=c(2,1,3)),2,rbind)
  score_phi <- -0.5*sig2u*sum(diag(prod_mat)) +
      0.5*sig2u*t(as.matrix(Vinv_resids)) %*%
      matrix(RH(Insp,RH(der_sig2u_Omega_phi,
                matrix(Vinv_resids,nrow=nrow(der_sig2u_Omega_phi)))),ncol=1)
  rm(prod_mat,Insp)

  if (trace){
     cat("score_phi: ",as.numeric(score_phi),"\n")
     cat("phi: ",as.numeric(phi),"\n")
   }

  return(as.numeric(score_phi))
}
###############################################################################
score_rho_phi <- function(param,sig2u,nsp,ntime,Wsp,y,X,Z,G.eff,
                          np.eff,b_fixed,trace=FALSE)
{
    rho <- param[1]/(1+abs(param[1])) #Ensures the interval (-1,1)
    phi <- param[2]/(1+abs(param[2]))

  call.Omega <- build_Omega_ar1(phi,ntime)
  Omega <- call.Omega$Omega
  der_sig2u_Omega_phi <- call.Omega$der_sig2u_Omega_phi
  Omegainv <- call.Omega$Omegainv
  rm(call.Omega)
  # score phi
  Omega.rep <- plyr::alply(replicate(nsp,Omega),3)
  V <- as.matrix(tcrossprod(Z %*% diag(sqrt(G.eff))) +
            sig2u*Matrix::bdiag(Omega.rep))
  #V <- as.matrix(AAt(Z %*% diag(sqrt(G.eff))) +
  #                 sig2u*Matrix::bdiag(Omega.rep))
  rm(Omega.rep)
  #chol.V <- chol(V)
  #chol.V <- lltA(V)
  #Vinv <- try(chol2inv(chol.V))
  Vinv <- lltAinv(V)
  if(class(Vinv) == "try-error")  Vinv <- MASS::ginv(as(V,"matrix"))
  Vinv.X <- Vinv %*% X
  #chol.XtVinvX <- chol(t(X)%*%Vinv.X)
  #chol.XtVinvX <- lltA(t(X)%*%Vinv.X)
  #XtVinvX.inv <- chol2inv(chol.XtVinvX)
  XtVinvX.inv <- lltAinv(t(X)%*%Vinv.X)
  P <- try(Vinv - Vinv.X %*% (XtVinvX.inv %*% t(Vinv.X)))
  if(class(P) == "try-error") {
      P <- Vinv - Vinv.X %*% (MASS::ginv(as(t(X)%*%Vinv.X,"matrix"))
                              %*% t(Vinv.X))
  }
  A <- as(Matrix::Diagonal(nsp)-rho*Wsp,"CsparseMatrix")
  It <- Matrix::Diagonal(ntime)
  A_It_y <- matrix(RH(A,RH(It, matrix(y,nrow=ntime))),ncol=1)
  W_It_y <- matrix(RH(Wsp,RH(It, matrix(y,nrow=ntime))),ncol=1)

  Vinv_resids <- Vinv%*%(A_It_y-X%*%b_fixed)
  Insp <- Matrix::Diagonal(nsp)
  prod_mat <- apply(aperm(RH(Insp,RH(der_sig2u_Omega_phi,
                array(P,dim=c(ncol(der_sig2u_Omega_phi),ncol(Insp),ncol(P))))),
                perm=c(2,1,3)),2,rbind)
  score_phi <- -0.5*sig2u*sum(diag(prod_mat)) +
      0.5*sig2u*t(as.matrix(Vinv_resids)) %*%
      matrix(RH(Insp,RH(der_sig2u_Omega_phi,
               matrix(Vinv_resids,nrow=nrow(der_sig2u_Omega_phi)))),ncol=1)
  rm(prod_mat,Insp)
  # der_V_phi <- sig2u*kronecker(Matrix::Diagonal(nsp),der_sig2u_Omega_phi)
  # score_phi <- -0.5*sum(diag(as.matrix(P%*%der_V_phi))) +
  #                 0.5*t(Vinv_resids)%*%der_V_phi%*%Vinv_resids
  score_rho <-  t(P%*%A_It_y) %*% W_It_y - ntime*sum(diag(solve(A,Wsp)))
  if (trace){
    cat("score_phi: ",as.numeric(score_phi),"\n")
    cat("score_rho: ",as.numeric(score_rho),"\n")
    cat("phi: ",phi,"\n")
    cat("rho: ",rho,"\n")
  }
  return(c(as.numeric(score_rho),as.numeric(score_phi)))
}
