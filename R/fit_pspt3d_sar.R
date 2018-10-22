fit_pspt3d_sar <- function(y,vary_init,sp1,sp2,time,Xfull,Zfull,Wsp = NULL,
                           nvarpar,nvarnopar,nvarspt,
                           #weights = NULL, GLAM = FALSE,
                           cspt,dsptlist,bdegspt,pordspt,nknotsspt,
                           cnopar,dnoparlist,bdegnopar,pordnopar,nknotsnopar,
                           names_varnopar,names_varpar,
                           psanova = TRUE,
                           f1_main = TRUE, f2_main = TRUE, ft_main = TRUE,
                           f12_int = FALSE, f1t_int = FALSE,
                           f2t_int = FALSE, f12t_int = FALSE,
                           sar = FALSE, rho_init = 0, rho_fixed = FALSE,
                           ar1 = FALSE, phi_init = 0, phi_fixed = FALSE,
                           bold = NULL, maxit = 20, thr = 1e-2,
                           trace = FALSE, var_num = FALSE)
{
  nsp <- length(sp1);  ntime <- length(time)
  if (!is.null(Wsp)) Wsp <- Matrix::Matrix(Wsp)
  In <- Matrix::Diagonal(nsp);  It <- Matrix::Diagonal(ntime)
  X <- Matrix::Matrix(Xfull)
  Z <- Matrix::Matrix(Zfull)
  if (!is.null(nknotsnopar)) nvarnopar <- length(nknotsnopar)
  if (!sar){
    rho_fixed <- TRUE
    rho_init <- 0
  } else if (is.null(rho_fixed)) rho_fixed <- FALSE
  if (!ar1){
    phi_fixed <- TRUE
    phi_init <- 0
  } else if (is.null(phi_fixed)) phi_fixed <- FALSE

  # Build vector and matrices for variance components in mixed model
  var_comp <- par_var_comp3d(la = var(as.numeric(y)), np_fixed = ncol(X),
                            pordspt = pordspt, cspt = cspt, dsptlist = dsptlist,
                            nvarnopar = nvarnopar,cnopar = cnopar,
                            pordnopar = pordnopar, dnoparlist = dnoparlist,
                            psanova = psanova, f1_main = f1_main,
                            f2_main = f2_main, ft_main = ft_main,
                            f12_int = f12_int, f1t_int = f1t_int,
                            f2t_int = f2t_int, f12t_int = f12t_int)

   # Number of parameters
   np <- var_comp$np; np_eff <- var_comp$np_eff
   # Vector of parameters
   la <- var_comp$la
   # Other vectors for SAP
   g1u <- var_comp$g1u; g2u <- var_comp$g2u; g3u <- var_comp$g3u
   g11b <- var_comp$g11b; g21b <- var_comp$g21b
   g12b <- var_comp$g12b; g31b <- var_comp$g31b
   g22b <- var_comp$g22b; g32b <- var_comp$g32b
   g1t <- var_comp$g1t; g2t <- var_comp$g2t; g3t <- var_comp$g3t
   G1inv.n <- var_comp$G1inv.n; G2inv.n <- var_comp$G2inv.n
   G3inv.n <- var_comp$G3inv.n
   if (psanova) {
       g12u <- var_comp$g12u; g21u <- var_comp$g21u
       g12b <- var_comp$g12b; g21b <- var_comp$g21b
       g13u <- var_comp$g13u; g31u <- var_comp$g31u
       g13b <- var_comp$g13b; g31b <- var_comp$g31b
       g23u <- var_comp$g23u; g32u <- var_comp$g32u
       g23b <- var_comp$g23b; g32b <- var_comp$g32b
       g123u <- var_comp$g123u; g213u <- var_comp$g213u
       g321u <- var_comp$g321u
       g123b <- var_comp$g123b; g213b <- var_comp$g213b
       g132b <- var_comp$g132b; g312b <- var_comp$g312b
       g231b <- var_comp$g231b; g321b <- var_comp$g321b
       G4inv.n <- var_comp$G4inv.n; G5inv.n <- var_comp$G5inv.n
       G6inv.n <- var_comp$G6inv.n; G7inv.n <- var_comp$G7inv.n
       G8inv.n <- var_comp$G8inv.n; G9inv.n <- var_comp$G9inv.n
       G10inv.n <- var_comp$G10inv.n; G11inv.n <- var_comp$G11inv.n
       G12inv.n <- var_comp$G12inv.n
   }
  # Do not remove var_comp, it is used in the next loop...


  # 0 0
  # 0 I
  D <- Matrix::Matrix(diag(c(rep(0,np_eff[1]), rep(1,sum(np_eff[-1])))))
  ed <- rep(0,length(la)-1)
  # add rho and phi to la vector of parameters
  la <- c(la,rho_init,phi_init)
  # Initialise the parameters
  if (is.null(bold)) bold = rep(0,sum(np_eff))
  eta <- X %*% bold[1:np_eff[1]] + Z %*% bold[-(1:np_eff[1])] #+ offset

  start <- proc.time()[3]
  for (iq in 1:maxit) { # Nested loops for SAP and phi (AR1) and rho (SAR)
      for (it in 1:maxit) {
	      rho <- la[length(la)-1]
	      phi <- la[length(la)]
	      if (is.null(Wsp)) {
	          A <- In
	          } else {
	          A <- In - rho*Wsp
	      }
	      # fast index ntime, slow index nsp
	      A_It_y <- matrix(RH(A,RH(It, array(y,dim=c(ncol(It),ncol(A))))),ncol=1)
	      # Es lo mismo que esta operación (lo he comprobado):
	      #      A_It_y <- matrix(kronecker(A,I_t)%*%y,ncol=1)
	      # Build covariance G matrix for random effects: block diagonal matrix
	      lG <- build_G3d(la = la, lg = var_comp,
	             nvarnopar = nvarnopar, dnoparlist = dnoparlist,
	             psanova = psanova, f1_main = f1_main, f2_main = f2_main,
	             ft_main = ft_main, f12_int = f12_int, f1t_int = f1t_int,
	             f2t_int = f2t_int, f12t_int = f12t_int)
	      G <- lG$G; Ginv <- lG$Ginv;
	      G_eff <- lG$G_eff; Ginv_eff <- lG$Ginv_eff
	      rm(lG)
	      sig2u <- la[1]
	      if (phi_fixed & phi_init==0) { # (Omega=I)
	        Xstar <- X
	        Zstar <- Z
	        A_It_ystar <- A_It_y } else {
	        # CONSTRUIMOS MATRIZ OMEGA (SI phi=0 ENTONCES OMEGA=I)
	        call.Omega <- build_Omega_ar1(phi,ntime)
	        Omega <- call.Omega$Omega
	        Omegainv <- call.Omega$Omegainv
	        chol.Omegainv <- lltA(Omegainv)
	        Xstar <- apply(aperm(RH(In,RH(chol.Omegainv,
	                 array(X,dim=c(ncol(chol.Omegainv),ncol(In),ncol(X))))),
	                 perm=c(2,1,3)),2,rbind)
	        Zstar <- apply(aperm(RH(In,RH(chol.Omegainv,
	                 array(Z,dim=c(ncol(chol.Omegainv),ncol(In),ncol(Z))))),
	                 perm=c(2,1,3)),2,rbind)
	        A_It_ystar <- matrix(RH(In,RH(chol.Omegainv,
	               array(A_It_y,dim=c(ncol(chol.Omegainv),ncol(In))))),ncol=1)
	        rm(call.Omega,Omegainv,chol.Omegainv)
	      }
	      if (is.null(weights)) {
	          w <- as.vector(matrix(1,nrow=length(y))) } else { w <- weights }
	      mat <- construct_matrices(Xstar,Zstar,A_It_ystar)
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
	      ltau_edf <- update_tau3d(la = la, lg = var_comp, G = G,
	                  dZtPZ_wide = dZtPZ_wide, brandom_wide = brandom_wide,
	                  psanova = psanova, f1_main = f1_main, f2_main = f2_main,
	                  ft_main = ft_main, f12_int = f12_int, f1t_int = f1t_int,
	                  f2t_int = f2t_int, f12t_int = f12t_int,
	                  nvarnopar = nvarnopar, dnoparlist = dnoparlist)

	      tau1 <- ltau_edf$tau1; tau2 <- ltau_edf$tau2; tau3 <- ltau_edf$tau3
	      ed1 <- ltau_edf$ed1; ed2 <- ltau_edf$ed2; ed3 <- ltau_edf$ed3
	      tauspt <- c(tau1,tau2,tau3); edfspt <- c(ed1,ed2,ed3)
	      names(tauspt) <- names(edfspt) <- c("sp1","sp2","time")
	      taunopar <- ltau_edf$taunopar; edfnopar <- ltau_edf$edfnopar
	      if (psanova){
	        tau4 <- ltau_edf$tau4; tau5 <- ltau_edf$tau5; tau6 <- ltau_edf$tau6;
	        tau7 <- ltau_edf$tau7; tau8 <- ltau_edf$tau8; tau9 <- ltau_edf$tau9;
	        tau10 <- ltau_edf$tau10; tau11 <- ltau_edf$tau11; tau12 <- ltau_edf$tau12
	        ed4 <- ltau_edf$ed4; ed5 <- ltau_edf$ed5; ed6 <- ltau_edf$ed6;
	        ed7 <- ltau_edf$ed7; ed8 <- ltau_edf$ed8; ed9 <- ltau_edf$ed9;
	        ed10 <- ltau_edf$ed10; ed11 <- ltau_edf$ed11; ed12 <- ltau_edf$ed12
	        tauspt <- c(tauspt,tau4,tau5,tau6,tau7,tau8,tau9,tau10,tau11,tau12)
	        edfspt <- c(edfspt,ed4,ed5,ed6,ed7,ed8,ed9,ed10,ed11,ed12)
	        names(tauspt) <- names(edfspt) <- c("f1_main","f2_main","ft_main",
	         "f12.1","f12.2","f1t.1","f1t.2","f2t.1","f2t.2",
	         "f12t.1","f12t.2","f12t.3")
	      }
	      # OJO: USA LAS VARIABLES STAR PARA CALCULAR SSR Y EDF'S
	      # Regression (Fahrmeir et al.) pp. 180
	      ssr <- as.numeric(mat$yty - t(c(bfixed, brandom)) %*% (2*mat$u - C %*% b))
	      # New variance components
	      lanew <- c(sig2u, tauspt)
	      edftot <- sum(edfspt) + ncol(Xstar)
	      if (nvarnopar>0){
	          lanew <- c(lanew,taunopar)
	          edftot <- edftot + sum(edfnopar)
	      }
	      if ((!rho_fixed) || !(rho==0))  edftot <- edftot + 1
	      if ((!phi_fixed) || !(phi==0))  edftot <- edftot + 1
	      sig2u <- as.numeric((ssr/(length(y) - edftot)))
	      # Update first component of la
	      lanew[1] <- sig2u
	      # Update last two components of la
	      lanew <- c(lanew,rho,phi)
	      dla <- mean(abs(la - lanew))
	      la <- lanew
	      if (trace) {
	        cat('\n Iteration SAP: ',it)
	        cat('\n sig2u ',la[1])
	        cat('\n edfspt:',round(edfspt,2))
	        if (!is.null(edfnopar)) cat('\n edfnopar: ',round(edfnopar,2),'\n')
	      }
	      #  convergence check
	      if (dla < thr) break
	    } # end for (it in 1:maxit)
        if(!rho_fixed & !phi_fixed) {
          par_init <- c(rho,phi)
          lb <- c(-Inf,-1)
          ub <- c(1,1)
        } else if (!rho_fixed & phi_fixed) {
          par_init <- c(rho,phi_init)
          lb <- c(-Inf,phi_init)
          ub <- c(1,phi_init)
        } else if (rho_fixed & !phi_fixed) {
          par_init <- c(rho_init,phi)
          lb <- c(rho_init,-1)
          ub <- c(rho_init,1)
        } else {
          rho <- rho_init
          phi <- phi_init
        }
      if (!rho_fixed | !phi_fixed){
          if (trace) print_level = 1 else print_level = 0
          rho_phi_optim <- nloptr::nloptr(x0=par_init,
                                   eval_f=llik_reml_fn3d,
                                   #eval_grad=score_phi,
                                   lb=lb, ub=ub,
                                   opts=list("algorithm"="NLOPT_LN_BOBYQA",
                                              "xtol_rel"=thr,
                                              "ftol_abs"=thr,
                                              "maxtime"=100,
                                              "print_level"=print_level),
                                    sig2u=sig2u,nsp=nsp,ntime=ntime,
                                    Wsp=Wsp,y=y,X=X,Z=Z,G_eff=G_eff,
                                    np_eff=np_eff,
                                    bfixed=bfixed,
                                    rho_fixed=rho_fixed,
                                    phi_fixed=phi_fixed)

        rho <- rho_phi_optim$solution[1]
        phi <- rho_phi_optim$solution[2]
     }
	    # Check Convergence in rho and phi parameters
	    rho_old <- la[length(la)-1]
	    phi_old <- la[length(la)]
	    drhophi <- mean(abs(c(rho_old,phi_old)-c(rho,phi)))
	    la[length(la)-1] <- rho
	    la[length(la)] <- phi
	    if (trace & any(c(sar,ar1))) {
	        cat("\n Iteration SAR-AR1: ",iq)
	        if (sar) cat("\n rho ",rho)
	        if (ar1) cat("\n phi",phi)
	     }
	     if (drhophi < thr) break
	  } # End loop for SAR-AR1
  end <- proc.time()[3]
  cat("\n Time to fit the model: ", (end-start), "seconds")
  eta <- Xstar%*%bfixed + Zstar%*%brandom #+ offset
  if (is.null(names_varnopar)) { edfnopar <- NULL; taunopar<-NULL }
#  FINAL ESTIMATES OF PARAMETERS
  sig2u <- la[1]
  # CREO QUE HABRÍA QUE ACTUALIZAR TAMBIÉN LOS EDF... REPASAR...
  if(!psanova) tauspt <- la[2:4] else tauspt <- la[2:13]
  rho <- la[length(la)-1]
  phi <- la[length(la)]

  # Valor de log.lik y log.lik.reml en el óptimo
  foptim <- llik_reml_var3d(c(rho,phi),
                          sig2u=sig2u,nsp=nsp,ntime=ntime,
                          Wsp=Wsp,y=y,X=X,Z=Z,G_eff=G_eff,
                          np_eff=npeff,bfixed=bfixed,
                          rho_fixed=rho_fixed,phi_fixed=phi_fixed)
  llik_reml <- foptim$llik_reml
  llik <- foptim$llik
  se_rho_an <- foptim$se_rho
  se_phi_an <- foptim$se_phi
  cov_rho_phi_an <- foptim$cov_rho_phi
  rm(foptim)
  if (var_num){
    hessian_optim <- numDeriv::hessian(llik_reml_fn,c(rho,phi),
                                       sig2u=sig2u,nsp=nsp,ntime=ntime,
                                       Wsp=Wsp,y=y,X=X,Z=Z,G_eff=G_eff,
                                       np_eff=np_eff,
                                       bfixed=bfixed,
                                       rho_fixed=rho_fixed,
                                       phi_fixed=phi_fixed)
    var_rho_phi_num <- solve(hessian_optim)
    rownames(var_rho_phi_num) <- colnames(var_rho_phi_num) <- c("rho","phi")
  } else var_rho_phi_num <- NULL
  ########################## CÁLCULO MATRIZ COVARIANZAS EFECTOS FIJOS Y ALEATORIOS
  # 		  # pp.375 Fahrmeir et al.
  #   # Se reescala A multiplicándola por sig2u toda la matriz
  #   # Matriz Covarianzas Bayesiana. OJO: Tiene en cuenta autocorrelación...
  if (phi_fixed & phi_init == 0) { # (Omega=I)
      Omega <- Omegainv <- It
  } else {
  # CONSTRUIMOS MATRIZ OMEGA (SI phi=0 ENTONCES OMEGA=I)
      call.Omega <- build_Omega_ar1(phi,ntime)
      Omega <- Matrix::Matrix(call.Omega$Omega)
      Omegainv <- Matrix::Matrix(call.Omega$Omegainv)
      rm(call.Omega)
  }
  Rinv <- kronecker(In,Omegainv)
  # 		  #chol.Rinv <- chol(Rinv)
  Uchol_Rinv <- Matrix::Matrix(lltA(as(Rinv,"matrix")))
  #chol_Rinv <- Matrix::Cholesky(Rinv)
  #factors_chol_Rinv <- Matrix::expand(chol_Rinv)
  #Lchol_Rinv <- Matrix::t(factors_chol_Rinv$P) %*% factors_chol_Rinv$L
  #Uchol_Rinv <- Matrix::t(factors_chol_Rinv$L) %*% factors_chol_Rinv$P
  #rm(chol_Rinv,factors_chol_Rinv)
  # Comprobación factor Choleski
  #Rinv2 <- Matrix::tcrossprod(Lchol_Rinv)
  #Rinv3 <- Matrix::crossprod(Uchol_Rinv)
  #range(Rinv - Rinv2); range(Rinv - Rinv3)
  #Xstar <- chol.Rinv %*% X
  #Zstar <- chol.Rinv %*% Z
  Xstar <- Uchol_Rinv %*% X
  Zstar <- Uchol_Rinv %*% Z
  # 		  #Astar <- as.matrix(chol.Rinv)%*%kronecker(A,It)
  start <- proc.time()[3]
  A_cov1 <- rbind(cbind(Matrix::crossprod(Xstar),
                        Matrix::t(Xstar) %*% Zstar),
                  cbind(Matrix::t(Zstar) %*% Xstar,
                        Matrix::crossprod(Zstar) +
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
  end <- proc.time()[3]
  cat("\n Time to compute covariances: ", (end-start), "seconds")
  # 		  # Matriz Covarianzas Frequentist tipo Sandwich (Algo menor las varianzas)
  # 		  # C.Rinv.C <- (1/sig2u)*rbind(cbind(t(X)%*%kronecker(Diagonal(nsp),Omegainv)%*%X,
  # 		  #                                   t(X)%*%kronecker(Diagonal(nsp),Omegainv)%*%Z),
  # 		  #                        cbind(t(Z)%*%kronecker(Diagonal(nsp),Omegainv)%*%X,
  # 		  #                              t(Z)%*%kronecker(Diagonal(nsp),Omegainv)%*%Z))
  # 		  C.Rinv.C <- (1/sig2u)*rbind(cbind(Matrix::crossprod(Xstar),
  # 		                                    Matrix::t(Xstar) %*% Zstar),
  # 		                              cbind(Matrix::t(Zstar) %*% Xstar,
  # 		                                    Matrix::crossprod(Zstar)))
  # 		  cov2.eff <- cov1_eff %*% C.Rinv.C %*% cov1_eff
  # CALCULAR VALORES ESTIMADOS Y RESIDUOS
  fit_Ay <- X %*% bfixed + Z %*% brandom # + offset
  fit <- try((Matrix::solve(A) %x% It) %*% fit_Ay)
  if (class(fit) == "try-error") {
           fit <- (MASS::ginv(as(A,"matrix")) %x% It) %*% eta  }
  ###### Commented Code: same result, much slower
  # var_fit_Ay <- cbind(X,Z) %*% cov1_eff %*% Matrix::t(cbind(X,Z))
  # var_fit <- try((Matrix::solve(A) %x% It) %*% var_fit_Ay %*%
  #  		                   (Matrix::t(Matrix::solve(A)) %x% It))
  # if (class(var_fit) == "try-error") { (MASS::ginv(as(A,"matrix")) %x% It) %*%
  #                       var_fit_Ay %*% (t(MASS::ginv(as(A,"matrix"))) %x% It) }
  # se_fit_Ay <- sqrt(diag(as.matrix(var_fit_Ay)))
  # se_fit <- sqrt(diag(as.matrix(var_fit)))
  se_fit_Ay <- Matrix::rowSums((cbind(X,Z) %*% cov1_eff) * cbind(X,Z))^0.5
  se_fit <- Matrix::rowSums(
                     (((Matrix::solve(A) %x% It) %*% cbind(X,Z))
                        %*% cov1_eff) *
                     ((Matrix::solve(A) %x% It) %*% cbind(X,Z)) )^0.5
  # resids N(0,sigma_u^2*Omega) and resids.norm N(0,sigma_eps*I)
  residuals <- as.vector((A %x% It) %*% y - fit_Ay)
  # REPASAR LOS RESIDS_NORM COMPARADO CON RESIDS.NORM
  #residuals_norm <- sqrt(1-phi^2)*(chol(as.matrix(Rinv)) %*% residuals)
  residuals_norm <- sqrt(1-phi^2)*(Uchol_Rinv %*% residuals)
  # Compute AIC y BIC based on loglik functions (Fahrmeir, pp. 664 and 677)
  aic <- -2*llik + 2*edftot
  bic <- -2*llik + log(length(y))*edftot

  res <- list(edfspt = edfspt, edfnopar = edfnopar, edftot = edftot,
              tauspt = tauspt, taunopar = taunopar,
              psanova = psanova, sar = sar, ar1 = ar1,
              fitted.values = as.vector(fit),
              fit_Ay = as.vector(fit_Ay),
              se_fitted.values = as.vector(se_fit),
              se_fit_Ay = as.vector(se_fit_Ay),
              residuals = as.vector(residuals),
              residuals_norm = as.vector(residuals_norm),
              sig2 = sig2u,
              rho = rho, phi = phi,
              se_rho = se_rho_an, se_phi = se_phi_an,
              cov_rho_phi_an = cov_rho_phi_an,
              var_rho_phi_num = var_rho_phi_num,
              bfixed = bfixed, brandom = brandom,
              se_bfixed = se_bfixed, se_brandom = se_brandom,
              llik = llik, llik_reml = llik_reml,
              aic = aic, bic = bic,
              vcov_b = cov1_eff,
              sp1 = sp1,sp2 = sp2, time = time)
} # end of function
