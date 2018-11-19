<<<<<<< HEAD
fit_pspt2d_sar <- function(y,vary_init,sp1,sp2,Xfull,Zfull,Wsp = NULL,
                           nvarpar,nvarnopar,nvarspt,
                           #weights = NULL, GLAM = FALSE,
                           cspt,dsptlist,bdegspt,pordspt,nknotsspt,
                           cnopar,dnoparlist,bdegnopar,pordnopar,nknotsnopar,
                           names_varnopar,names_varpar,
                           psanova, f1_main, f2_main,
                           f12_int, sar = FALSE, rho_init = 0,
                           rho_fixed = FALSE,
                           bold = NULL, maxit = 30, thr = 1e-2, trace=FALSE,
                           var_num = FALSE)
{
  nsp_full <- length(sp1)
  if (!is.null(Wsp)) {
    Wsp <- Matrix::Matrix(Wsp)
    nsp_short <- nrow(Wsp)
    if ((nsp_full %% nsp_short) != 0) stop("Dimensions of W and spatial coordinates not compatibles")
    Wsp_full <-   Wsp %x% Matrix::Diagonal(nsp_full/nsp_short)
  } else Wsp_full <- NULL
  In <- Matrix::Diagonal(nsp_full)
  X <- Matrix::Matrix(Xfull)
  Z <- Matrix::Matrix(Zfull)
  if (!is.null(nknotsnopar)) nvarnopar <- length(nknotsnopar)
  if (!sar) rho_fixed <- TRUE

  # Build vector and matrices for variance components in mixed model
  var_comp <- par_var_comp2d(la = var(as.numeric(y)), np_fixed = ncol(X),
                            pordspt = pordspt, cspt = cspt, dsptlist = dsptlist,
                            nvarnopar = nvarnopar,cnopar = cnopar,
                            pordnopar = pordnopar, dnoparlist = dnoparlist,
                            psanova = psanova, f1_main = f1_main,
                            f2_main = f2_main,f12_int = f12_int)

   # Number of parameters
   np <- var_comp$np; np_eff <- var_comp$np_eff
   # Vector of parameters
   la <- var_comp$la
   # Other vectors for SAP
   g1u <- var_comp$g1u; g2u <- var_comp$g2u
   g1b <- var_comp$g1b; g2b <- var_comp$g2b
   G1inv.n <- var_comp$G1inv.n; G2inv.n <- var_comp$G2inv.n
   if (psanova) {
       g1v <- var_comp$g1v; g2v <- var_comp$g2v
       G3inv.n <- var_comp$G3inv.n; G4inv.n <- var_comp$G4inv.n
   }
  # Do not remove var_comp, it is used in the next loop...


  # 0 0
  # 0 I
  D <- Matrix::Matrix(diag(c(rep(0,np_eff[1]), rep(1,sum(np_eff[-1])))))
  ed <- rep(0,length(la)-1)
  # add rho and phi to la vector of parameters
  la <- c(la,rho_init)
  # Initialise the parameters
  if (is.null(bold)) bold = rep(0,sum(np_eff))
  eta <- X %*% bold[1:np_eff[1]] + Z %*% bold[-(1:np_eff[1])] #+ offset

  start <- proc.time()[3]
  for (iq in 1:maxit) { # Nested loops for SAP and rho (for SAR case)
      for (it in 1:maxit) {
	      rho <- la[length(la)]
	      if (is.null(Wsp_full)) { A <- In } else { A <- In - rho*Wsp_full }
	      A_y <- matrix( A %*% y)
	      # Build covariance G matrix for random effects: block diagonal matrix
	      lG <- build_G2d(la = la, lg = var_comp,
	             nvarnopar = nvarnopar, dnoparlist = dnoparlist,
	             psanova = psanova, f1_main = f1_main, f2_main = f2_main,
	             f12_int = f12_int)
	      G <- lG$G; Ginv <- lG$Ginv;
	      G_eff <- lG$G_eff; Ginv_eff <- lG$Ginv_eff
	      rm(lG)
	      sig2u <- la[1]
	      #if (is.null(weights)) {
	      #    w <- as.vector(matrix(1,nrow=length(y))) } else { w <- weights }
	      mat <- construct_matrices(X,Z,A_y)
	      C <- Matrix::Matrix( construct_block(mat$XtX, Matrix::t(mat$ZtX*G_eff),
	                           mat$ZtX, Matrix::t(mat$ZtZ*G_eff)) )

	      # ES LA MATRIZ C DE COEFICIENTES DEL SISTEMA (12) EN PAPER SAP
	      # NO SE PUEDE RESOLVER SI SE INCLUYEN LOS ZEROS EN MATRICES X,Z,G
	      Hinv <- try(Matrix::solve((1/la[1])*C + D,tol = 1e-40))
	      if (class(Hinv) == "try-error")
	        Hinv <- MASS::ginv((1/la[1])*as.matrix(C) + as.matrix(D))
	      b <- as.vector( (1/la[1])*Hinv %*% mat$u )
	      bfixed <- b[1:np_eff[1]]
	      names(bfixed) <- gsub("X_","",colnames(X))
	      names(bfixed) <- paste("fixed_",names(bfixed),sep="")
	      brandom <- G_eff*b[-(1:np_eff[1])]
	      names(brandom) <- gsub("Z_","",colnames(Z))
	      names(brandom) <- paste("random_",names(brandom),sep="")

	      # Compute effective dimensions and variances
	      # Only the diagonal of ZtPZ
	      dZtPZ <- 1/la[1]*apply((Matrix::t(Hinv[-(1:np_eff[1]),])*mat$ZtXtZ),2,sum)
	      dZtPZ_wide <- rep(0,length(G))
	      brandom_wide <- rep(0,length(G))
	      index.zeros.G <- G==0
	      index <- 1
	      for (ide in 1:length(G)) {
	        if (!index.zeros.G[ide]) {
	          brandom_wide[ide] <- brandom[index]
	          dZtPZ_wide[ide] <- dZtPZ[index]
	          index <- index + 1
	        }
	      }
	      # MIRAR AÑADIR tau_init, tau_fixed, taunopar_fixed y taunoparinit
	      ltau_edf <- update_tau2d(la = la, lg = var_comp, G = G,
	                  dZtPZ_wide = dZtPZ_wide, brandom_wide = brandom_wide,
	                  psanova = psanova, f1_main = f1_main, f2_main = f2_main,
	                  f12_int = f12_int,
	                  nvarnopar = nvarnopar, dnoparlist = dnoparlist)

	      tau1 <- ltau_edf$tau1; tau2 <- ltau_edf$tau2
	      ed1 <- ltau_edf$ed1; ed2 <- ltau_edf$ed2
	      tauspt <- c(tau1,tau2); edfspt <- c(ed1,ed2)
	      names(tauspt) <- names(edfspt) <- c("sp1","sp2")
	      taunopar <- ltau_edf$taunopar; edfnopar <- ltau_edf$edfnopar
	      if (psanova){
	        tau3 <- ltau_edf$tau4; tau4 <- ltau_edf$tau4
	        ed3 <- ltau_edf$ed3; ed4 <- ltau_edf$ed4
	        tauspt <- c(tauspt,tau3,tau4)
	        edfspt <- c(edfspt,ed3,ed4)
	        names(tauspt) <- names(edfspt) <- c("f1_main","f2_main",
	         "f12.1","f12.2")
	      }
	      # Regression (Fahrmeir et al.) pp. 180
	      ssr <- as.numeric(mat$yty - t(c(bfixed, brandom)) %*% (2*mat$u - C %*% b))
	      # New variance components
	      lanew <- c(sig2u, tauspt)
	      edftot <- sum(edfspt) + ncol(X)
	      if (nvarnopar>0){
	          lanew <- c(lanew,taunopar)
	          edftot <- edftot + sum(edfnopar)
	      }
	      if ((!rho_fixed) | (rho!=0))  edftot <- edftot + 1
	      sig2u <- as.numeric((ssr/(length(y) - edftot)))
	      # Update first component of la
	      lanew[1] <- sig2u
	      # Update rho
	      lanew <- c(lanew,rho)
	      dla <- mean(abs(la - lanew))
	      la <- lanew
	      if (trace) {
	        cat('\n Iteration SAP: ',it)
	        cat('\n sig2u ',la[1])
	        cat('\n edfspt:',round(edfspt,2))
	        if (!is.null(edfnopar)) cat('\n edfnopar: ',round(edfnopar,2))
	        }
	      #  convergence check
	      if (dla < thr) break
	    } # end for (it in 1:maxit)
      if (!rho_fixed) {
         par_init <- c(rho)
         lb <- c(-Inf)
         ub <- c(1)
      } else { rho <- rho_init }
      if (!rho_fixed){
          rho_phi_optim <- nloptr::nloptr(x0=par_init,
                                   eval_f=llik_reml_fn2d,
                                   #eval_grad=score_phi,
                                   lb=lb, ub=ub,
                                   opts=list("algorithm"="NLOPT_LN_BOBYQA",
                                              "xtol_rel"=1.0e-1,
                                              "ftol_abs"=1.0e-2,
                                              "maxtime"=100,
                                              "print_level"=0),
                                    sig2u=sig2u,nsp=nsp_full,
                                    Wsp=Wsp_full,y=y,X=X,Z=Z,G_eff=G_eff,
                                    np_eff=np_eff,
                                    bfixed=bfixed,
                                    rho_fixed=rho_fixed)

        rho <- rho_phi_optim$solution[1]
     }
	   # Check Convergence in rho and phi parameters
	   rho_old <- la[length(la)]
	   drho <- abs(rho_old - rho)
	   la[length(la)] <- rho
	   if (trace & sar){
	       cat("\n Iteration SAR: ",iq)
	       cat("\n rho ",rho)
	   }
	   if (drho < thr) break
	} # End loop for SAR
  #if (trace) {
     end <- proc.time()[3]
		 cat("\n Time to fit the model: ", (end-start), "seconds")
	#}
  eta <- X %*% bfixed + Z %*% brandom #+ offset
  if (is.null(names_varnopar)) { edfnopar <- NULL; taunopar<-NULL }

#  FINAL ESTIMATES OF PARAMETERS
  sig2u <- la[1]
  if(!psanova) tauspt <- la[2:3] else tauspt <- la[2:5]
  rho <- la[length(la)]
  # Valor de log.lik y log.lik.reml en el óptimo
  foptim <- llik_reml_var2d(c(rho),
                          sig2u=sig2u,nsp=nsp_full,
                          Wsp=Wsp_full,y=y,X=X,Z=Z,G_eff=G_eff,
                          np_eff=np_eff,bfixed=bfixed,rho_fixed=rho_fixed)
  llik_reml <- foptim$llik_reml
  llik <- foptim$llik
  se_rho_an <- foptim$se_rho
  rm(foptim)
  if (var_num){
    hessian_optim <- numDeriv::hessian(llik_reml_fn2d,c(rho),
                                       sig2u=sig2u,nsp=nsp_full,
                                       Wsp=Wsp_full,y=y,X=X,Z=Z,G_eff=G_eff,
                                       np_eff=np_eff,
                                       bfixed=bfixed,
                                       rho_fixed=rho_fixed)
    var_rho_num <- as.numeric( 1 / hessian_optim )
    se_rho_num <- sqrt(var_rho_num)
  } else { var_rho_num <- se_rho_num <- NULL }
  ########################## CÁLCULO MATRIZ COVARIANZAS EFECTOS FIJOS Y ALEATORIOS
  # 		  # pp.375 Fahrmeir et al.
  #   # Se reescala A multiplicándola por sig2u toda la matriz
  #   # Matriz Covarianzas Bayesiana.
  if (trace) start <- proc.time()[3]
  A_cov1 <- rbind(cbind(Matrix::crossprod(X),
                        Matrix::t(X) %*% Z),
                  cbind(Matrix::t(Z) %*% X,
                        Matrix::crossprod(Z) +
                        sig2u*Matrix::Matrix(diag(Ginv_eff)) ))
  cov1_eff <- try( sig2u*Matrix::solve(A_cov1) )
  if (class(cov1_eff) == "try-error") {
        cov1_eff <- try( sig2u*solve(as.matrix(A_cov1),tol=1e-40) )
        if (class(cov1_eff) == "try-error") {
          cov1_eff <- sig2u*MASS::ginv(as.matrix(A_cov1))
          cov1_eff.ginv <- TRUE
        } else { cov1_eff.ginv <- FALSE }
  }
  cov1_eff <- Matrix::Matrix(cov1_eff)
  rownames(cov1_eff) <- colnames(cov1_eff) <- c(names(bfixed),names(brandom))
  se_bfixed <- sqrt(diag(as.matrix(cov1_eff[names(bfixed),names(bfixed)])))
  names(se_bfixed) <- names(bfixed)
  se_brandom <- sqrt(diag(as.matrix(cov1_eff[names(brandom),names(brandom)])))
  names(se_brandom) <- names(brandom)
  #if (trace) {
    end <- proc.time()[3]
    cat("\n Time to compute covariances: ", (end-start), "seconds")
  #}
  # 		  # Matriz Covarianzas Frequentist tipo Sandwich (Algo menor las varianzas)
  # 		  C.Rinv.C <- (1/sig2u)*rbind(cbind(Matrix::crossprod(X),
  # 		                                    Matrix::t(X) %*% Z),
  # 		                              cbind(Matrix::t(Z) %*% X,
  # 		                                    Matrix::crossprod(Z)))
  # 		  cov2.eff <- cov1_eff %*% C.Rinv.C %*% cov1_eff
  # CALCULAR VALORES ESTIMADOS Y RESIDUOS
  fit_Ay <- X %*% bfixed + Z %*% brandom # + offset
  fit <- try( (Matrix::solve(A,fit_Ay)) )
  if (class(fit) == "try-error") fit <- MASS::ginv(as.matrix(A)) %*% eta
  ### Commented code: much slower same result
  # var_fit_Ay <- cbind(X,Z) %*% cov1_eff %*% Matrix::t(cbind(X,Z))
  # var_fit <- try( (Matrix::solve(A,var_fit_Ay) %*%
  #  		                   (Matrix::t(Matrix::solve(A)))) )
  # if (class(var_fit) == "try-error") {
  #         MASS::ginv(as.matrix(A)) %*% var_fit_Ay %*%
  #                   (t(MASS::ginv(as.matrix(A))))
  # }
  # se_fit_Ay <- sqrt(diag(as.matrix(var_fit_Ay)))
  # se_fit <- sqrt(diag(as.matrix(var_fit)))
  se_fit_Ay <- Matrix::rowSums((cbind(X,Z) %*% cov1_eff) * cbind(X,Z))^0.5
  se_fit <- Matrix::rowSums((Matrix::solve(A,cbind(X,Z)) %*% cov1_eff) *
                              (Matrix::solve(A,cbind(X,Z))))^0.5

  # resids N(0,sigma_u^2*Omega) and resids.norm N(0,sigma_eps*I)
  residuals <- as.vector( (A %*% y) - fit_Ay)
  # Compute AIC y BIC based on loglik functions (Fahrmeir, pp. 664 and 677)
  aic <- -2*llik + 2*edftot
  bic <- -2*llik + log(length(y))*edftot

  res <- list(edfspt = edfspt, edfnopar = edfnopar, edftot = edftot,
              tauspt = tauspt, taunopar = taunopar,
              psanova = psanova, sar = sar,
              fitted.values = as.vector(fit),
              fit_Ay = as.vector(fit_Ay),
              se_fitted.values = as.vector(se_fit),
              se_fit_Ay = as.vector(se_fit_Ay),
              residuals = as.vector(residuals),
              sig2 = sig2u,
              rho = rho, se_rho = se_rho_an,
              se_rho_num = se_rho_num,
              bfixed = bfixed, brandom = brandom,
              se_bfixed = se_bfixed, se_brandom = se_brandom,
              llik = llik, llik_reml = llik_reml,
              aic = aic, bic = bic,
              vcov_b = cov1_eff,
              sp1 = sp1,sp2 = sp2, time = NULL)
} # end of function
=======
fit_pspt2d_sar <- function(y,vary_init,sp1,sp2,Xfull,Zfull,Wsp = NULL,
                           nvarpar,nvarnopar,nvarspt,
                           #weights = NULL, GLAM = FALSE,
                           cspt,dsptlist,bdegspt,pordspt,nknotsspt,
                           cnopar,dnoparlist,bdegnopar,pordnopar,nknotsnopar,
                           names_varnopar,names_varpar,
                           psanova, f1_main, f2_main,
                           f12_int, sar = FALSE, rho_init = 0,
                           rho_fixed = FALSE,
                           bold = NULL, maxit = 30, thr = 1e-2, trace=FALSE,
                           var_num = FALSE)
{
  nsp_full <- length(sp1)
  if (!is.null(Wsp)) {
    Wsp <- Matrix::Matrix(Wsp)
    nsp_short <- nrow(Wsp)
    if ((nsp_full %% nsp_short) != 0) stop("Dimensions of W and spatial coordinates not compatibles")
    Wsp_full <-   Wsp %x% Matrix::Diagonal(nsp_full/nsp_short)
  } else Wsp_full <- NULL
  In <- Matrix::Diagonal(nsp_full)
  X <- Matrix::Matrix(Xfull)
  Z <- Matrix::Matrix(Zfull)
  if (!is.null(nknotsnopar)) nvarnopar <- length(nknotsnopar)
  if (!sar) rho_fixed <- TRUE

  # Build vector and matrices for variance components in mixed model
  var_comp <- par_var_comp2d(la = var(as.numeric(y)), np_fixed = ncol(X),
                            pordspt = pordspt, cspt = cspt, dsptlist = dsptlist,
                            nvarnopar = nvarnopar,cnopar = cnopar,
                            pordnopar = pordnopar, dnoparlist = dnoparlist,
                            psanova = psanova, f1_main = f1_main,
                            f2_main = f2_main,f12_int = f12_int)

   # Number of parameters
   np <- var_comp$np; np_eff <- var_comp$np_eff
   # Vector of parameters
   la <- var_comp$la
   # Other vectors for SAP
   g1u <- var_comp$g1u; g2u <- var_comp$g2u
   g1b <- var_comp$g1b; g2b <- var_comp$g2b
   G1inv.n <- var_comp$G1inv.n; G2inv.n <- var_comp$G2inv.n
   if (psanova) {
       g1v <- var_comp$g1v; g2v <- var_comp$g2v
       G3inv.n <- var_comp$G3inv.n; G4inv.n <- var_comp$G4inv.n
   }
  # Do not remove var_comp, it is used in the next loop...


  # 0 0
  # 0 I
  D <- Matrix::Matrix(diag(c(rep(0,np_eff[1]), rep(1,sum(np_eff[-1])))))
  ed <- rep(0,length(la)-1)
  # add rho and phi to la vector of parameters
  la <- c(la,rho_init)
  # Initialise the parameters
  if (is.null(bold)) bold = rep(0,sum(np_eff))
  eta <- X %*% bold[1:np_eff[1]] + Z %*% bold[-(1:np_eff[1])] #+ offset

  start <- proc.time()[3]
  for (iq in 1:maxit) { # Nested loops for SAP and rho (for SAR case)
      for (it in 1:maxit) {
	      rho <- la[length(la)]
	      if (is.null(Wsp_full)) { A <- In } else { A <- In - rho*Wsp_full }
	      A_y <- matrix( A %*% y)
	      # Build covariance G matrix for random effects: block diagonal matrix
	      lG <- build_G2d(la = la, lg = var_comp,
	             nvarnopar = nvarnopar, dnoparlist = dnoparlist,
	             psanova = psanova, f1_main = f1_main, f2_main = f2_main,
	             f12_int = f12_int)
	      G <- lG$G; Ginv <- lG$Ginv;
	      G_eff <- lG$G_eff; Ginv_eff <- lG$Ginv_eff
	      rm(lG)
	      sig2u <- la[1]
	      #if (is.null(weights)) {
	      #    w <- as.vector(matrix(1,nrow=length(y))) } else { w <- weights }
	      mat <- construct_matrices(X,Z,A_y)
	      C <- Matrix::Matrix( construct_block(mat$XtX, Matrix::t(mat$ZtX*G_eff),
	                           mat$ZtX, Matrix::t(mat$ZtZ*G_eff)) )

	      # ES LA MATRIZ C DE COEFICIENTES DEL SISTEMA (12) EN PAPER SAP
	      # NO SE PUEDE RESOLVER SI SE INCLUYEN LOS ZEROS EN MATRICES X,Z,G
	      Hinv <- try(Matrix::solve((1/la[1])*C + D,tol = 1e-40))
	      if (class(Hinv) == "try-error")
	        Hinv <- MASS::ginv((1/la[1])*as.matrix(C) + as.matrix(D))
	      b <- as.vector( (1/la[1])*Hinv %*% mat$u )
	      bfixed <- b[1:np_eff[1]]
	      names(bfixed) <- gsub("X_","",colnames(X))
	      names(bfixed) <- paste("fixed_",names(bfixed),sep="")
	      brandom <- G_eff*b[-(1:np_eff[1])]
	      names(brandom) <- gsub("Z_","",colnames(Z))
	      names(brandom) <- paste("random_",names(brandom),sep="")

	      # Compute effective dimensions and variances
	      # Only the diagonal of ZtPZ
	      dZtPZ <- 1/la[1]*apply((Matrix::t(Hinv[-(1:np_eff[1]),])*mat$ZtXtZ),2,sum)
	      dZtPZ_wide <- rep(0,length(G))
	      brandom_wide <- rep(0,length(G))
	      index.zeros.G <- G==0
	      index <- 1
	      for (ide in 1:length(G)) {
	        if (!index.zeros.G[ide]) {
	          brandom_wide[ide] <- brandom[index]
	          dZtPZ_wide[ide] <- dZtPZ[index]
	          index <- index + 1
	        }
	      }
	      # MIRAR AÑADIR tau_init, tau_fixed, taunopar_fixed y taunoparinit
	      ltau_edf <- update_tau2d(la = la, lg = var_comp, G = G,
	                  dZtPZ_wide = dZtPZ_wide, brandom_wide = brandom_wide,
	                  psanova = psanova, f1_main = f1_main, f2_main = f2_main,
	                  f12_int = f12_int,
	                  nvarnopar = nvarnopar, dnoparlist = dnoparlist)

	      tau1 <- ltau_edf$tau1; tau2 <- ltau_edf$tau2
	      ed1 <- ltau_edf$ed1; ed2 <- ltau_edf$ed2
	      tauspt <- c(tau1,tau2); edfspt <- c(ed1,ed2)
	      names(tauspt) <- names(edfspt) <- c("sp1","sp2")
	      taunopar <- ltau_edf$taunopar; edfnopar <- ltau_edf$edfnopar
	      if (psanova){
	        tau3 <- ltau_edf$tau4; tau4 <- ltau_edf$tau4
	        ed3 <- ltau_edf$ed3; ed4 <- ltau_edf$ed4
	        tauspt <- c(tauspt,tau3,tau4)
	        edfspt <- c(edfspt,ed3,ed4)
	        names(tauspt) <- names(edfspt) <- c("f1_main","f2_main",
	         "f12.1","f12.2")
	      }
	      # Regression (Fahrmeir et al.) pp. 180
	      ssr <- as.numeric(mat$yty - t(c(bfixed, brandom)) %*% (2*mat$u - C %*% b))
	      # New variance components
	      lanew <- c(sig2u, tauspt)
	      edftot <- sum(edfspt) + ncol(X)
	      if (nvarnopar>0){
	          lanew <- c(lanew,taunopar)
	          edftot <- edftot + sum(edfnopar)
	      }
	      if ((!rho_fixed) | (rho!=0))  edftot <- edftot + 1
	      sig2u <- as.numeric((ssr/(length(y) - edftot)))
	      # Update first component of la
	      lanew[1] <- sig2u
	      # Update rho
	      lanew <- c(lanew,rho)
	      dla <- mean(abs(la - lanew))
	      la <- lanew
	      if (trace) {
	        cat('\n Iteration SAP: ',it)
	        cat('\n sig2u ',la[1])
	        cat('\n edfspt:',round(edfspt,2))
	        if (!is.null(edfnopar)) cat('\n edfnopar: ',round(edfnopar,2))
	        }
	      #  convergence check
	      if (dla < thr) break
	    } # end for (it in 1:maxit)
      if (!rho_fixed) {
         par_init <- c(rho)
         lb <- c(-Inf)
         ub <- c(1)
      } else { rho <- rho_init }
      if (!rho_fixed){
          rho_phi_optim <- nloptr::nloptr(x0=par_init,
                                   eval_f=llik_reml_fn2d,
                                   #eval_grad=score_phi,
                                   lb=lb, ub=ub,
                                   opts=list("algorithm"="NLOPT_LN_BOBYQA",
                                              "xtol_rel"=1.0e-1,
                                              "ftol_abs"=1.0e-2,
                                              "maxtime"=100,
                                              "print_level"=0),
                                    sig2u=sig2u,nsp=nsp_full,
                                    Wsp=Wsp_full,y=y,X=X,Z=Z,G_eff=G_eff,
                                    np_eff=np_eff,
                                    bfixed=bfixed,
                                    rho_fixed=rho_fixed)

        rho <- rho_phi_optim$solution[1]
     }
	   # Check Convergence in rho and phi parameters
	   rho_old <- la[length(la)]
	   drho <- abs(rho_old - rho)
	   la[length(la)] <- rho
	   if (trace & sar){
	       cat("\n Iteration SAR: ",iq)
	       cat("\n rho ",rho)
	   }
	   if (drho < thr) break
	} # End loop for SAR
  #if (trace) {
     end <- proc.time()[3]
		 cat("\n Time to fit the model: ", (end-start), "seconds")
	#}
  eta <- X %*% bfixed + Z %*% brandom #+ offset
  if (is.null(names_varnopar)) { edfnopar <- NULL; taunopar<-NULL }

#  FINAL ESTIMATES OF PARAMETERS
  sig2u <- la[1]
  if(!psanova) tauspt <- la[2:3] else tauspt <- la[2:5]
  rho <- la[length(la)]
  # Valor de log.lik y log.lik.reml en el óptimo
  foptim <- llik_reml_var2d(c(rho),
                          sig2u=sig2u,nsp=nsp_full,
                          Wsp=Wsp_full,y=y,X=X,Z=Z,G_eff=G_eff,
                          np_eff=np_eff,bfixed=bfixed,rho_fixed=rho_fixed)
  llik_reml <- foptim$llik_reml
  llik <- foptim$llik
  se_rho_an <- foptim$se_rho
  rm(foptim)
  if (var_num){
    hessian_optim <- numDeriv::hessian(llik_reml_fn2d,c(rho),
                                       sig2u=sig2u,nsp=nsp_full,
                                       Wsp=Wsp_full,y=y,X=X,Z=Z,G_eff=G_eff,
                                       np_eff=np_eff,
                                       bfixed=bfixed,
                                       rho_fixed=rho_fixed)
    var_rho_num <- as.numeric( 1 / hessian_optim )
    se_rho_num <- sqrt(var_rho_num)
  } else { var_rho_num <- se_rho_num <- NULL }
  ########################## CÁLCULO MATRIZ COVARIANZAS EFECTOS FIJOS Y ALEATORIOS
  # 		  # pp.375 Fahrmeir et al.
  #   # Se reescala A multiplicándola por sig2u toda la matriz
  #   # Matriz Covarianzas Bayesiana.
  if (trace) start <- proc.time()[3]
  A_cov1 <- rbind(cbind(Matrix::crossprod(X),
                        Matrix::t(X) %*% Z),
                  cbind(Matrix::t(Z) %*% X,
                        Matrix::crossprod(Z) +
                        sig2u*Matrix::Matrix(diag(Ginv_eff)) ))
  cov1_eff <- try( sig2u*Matrix::solve(A_cov1) )
  if (class(cov1_eff) == "try-error") {
        cov1_eff <- try( sig2u*solve(as.matrix(A_cov1),tol=1e-40) )
        if (class(cov1_eff) == "try-error") {
          cov1_eff <- sig2u*MASS::ginv(as.matrix(A_cov1))
          cov1_eff.ginv <- TRUE
        } else { cov1_eff.ginv <- FALSE }
  }
  cov1_eff <- Matrix::Matrix(cov1_eff)
  rownames(cov1_eff) <- colnames(cov1_eff) <- c(names(bfixed),names(brandom))
  se_bfixed <- sqrt(diag(as.matrix(cov1_eff[names(bfixed),names(bfixed)])))
  names(se_bfixed) <- names(bfixed)
  se_brandom <- sqrt(diag(as.matrix(cov1_eff[names(brandom),names(brandom)])))
  names(se_brandom) <- names(brandom)
  #if (trace) {
    end <- proc.time()[3]
    cat("\n Time to compute covariances: ", (end-start), "seconds")
  #}
  # 		  # Matriz Covarianzas Frequentist tipo Sandwich (Algo menor las varianzas)
  # 		  C.Rinv.C <- (1/sig2u)*rbind(cbind(Matrix::crossprod(X),
  # 		                                    Matrix::t(X) %*% Z),
  # 		                              cbind(Matrix::t(Z) %*% X,
  # 		                                    Matrix::crossprod(Z)))
  # 		  cov2.eff <- cov1_eff %*% C.Rinv.C %*% cov1_eff
  # CALCULAR VALORES ESTIMADOS Y RESIDUOS
  fit_Ay <- X %*% bfixed + Z %*% brandom # + offset
  fit <- try( (Matrix::solve(A,fit_Ay)) )
  if (class(fit) == "try-error") fit <- MASS::ginv(as.matrix(A)) %*% eta
  ### Commented code: much slower same result
  # var_fit_Ay <- cbind(X,Z) %*% cov1_eff %*% Matrix::t(cbind(X,Z))
  # var_fit <- try( (Matrix::solve(A,var_fit_Ay) %*%
  #  		                   (Matrix::t(Matrix::solve(A)))) )
  # if (class(var_fit) == "try-error") {
  #         MASS::ginv(as.matrix(A)) %*% var_fit_Ay %*%
  #                   (t(MASS::ginv(as.matrix(A))))
  # }
  # se_fit_Ay <- sqrt(diag(as.matrix(var_fit_Ay)))
  # se_fit <- sqrt(diag(as.matrix(var_fit)))
  se_fit_Ay <- Matrix::rowSums((cbind(X,Z) %*% cov1_eff) * cbind(X,Z))^0.5
  se_fit <- Matrix::rowSums((Matrix::solve(A,cbind(X,Z)) %*% cov1_eff) *
                              (Matrix::solve(A,cbind(X,Z))))^0.5

  # resids N(0,sigma_u^2*Omega) and resids.norm N(0,sigma_eps*I)
  residuals <- as.vector( (A %*% y) - fit_Ay)
  # Compute AIC y BIC based on loglik functions (Fahrmeir, pp. 664 and 677)
  aic <- -2*llik + 2*edftot
  bic <- -2*llik + log(length(y))*edftot

  res <- list(edfspt = edfspt, edfnopar = edfnopar, edftot = edftot,
              tauspt = tauspt, taunopar = taunopar,
              psanova = psanova, sar = sar,
              fitted.values = as.vector(fit),
              fit_Ay = as.vector(fit_Ay),
              se_fitted.values = as.vector(se_fit),
              se_fit_Ay = as.vector(se_fit_Ay),
              residuals = as.vector(residuals),
              sig2 = sig2u,
              rho = rho, se_rho = se_rho_an,
              se_rho_num = se_rho_num,
              bfixed = bfixed, brandom = brandom,
              se_bfixed = se_bfixed, se_brandom = se_brandom,
              llik = llik, llik_reml = llik_reml,
              aic = aic, bic = bic,
              vcov_b = cov1_eff,
              sp1 = sp1,sp2 = sp2, time = NULL)
} # end of function
>>>>>>> 3dd5daf8ee882992716e3954f4de5c576169b2fa
