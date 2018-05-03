fit_ps3d_sar <- function(y,x1,x2,x3,decom=2,Xpar=NULL,Xnopar=NULL,
    psanova=FALSE,f1_main=TRUE,f2_main=TRUE,f3_main=TRUE,
    f12_int=TRUE,f13_int=TRUE,f23_int=TRUE,f123_int=TRUE,
    div_sp=c(1,1),div_tm=c(1,1),
    knots_sp=c(10,10),pord_sp=c(2,2),bdeg_sp=c(3,3),
    knots_time=10,pord_time=2,bdeg_time=3,
    x1lim=NULL,x2lim=NULL,x3lim=NULL,
    vary_init=1,rho_init=0,rho_fixed=TRUE,phi_init=0,phi_fixed=TRUE,
    tau_init = rep(1,12), tau_fixed = rep(FALSE,12),
    GLAM = FALSE,trace = FALSE,thr=1e-3,maxit=200,offset=0,
    bold = NULL,weights=NULL,Wsp=NULL,
    knots_nopar=NULL,bdeg_nopar=NULL,pord_nopar=NULL,
    tau_nopar_init=NULL,tau_nopar_fixed=FALSE,post_estim=TRUE,sc_rho_phi=NULL) {
# Function to estimate PS-ANOVA-SAR in 3d.
# x1 and x2 are spatial coordinates and x3 is time coordinate
# It is allowed to select the terms which enter in PS-ANOVA (main effects,
#  interactions 2d or interactions 3d).
# It is allowed to use nested basis for interactions terms to reduce
# the number of knots for interactions.
# The SAR specification is chosen if rho_fixed=FALSE (by default)
# The parametric covariates are included in Xpar matrix and
# non-parametric additive covariates are included in Xnopar matrix

  nsp <- length(x1);  ntime <- length(x3)
  In <- Matrix::Diagonal(nsp)
  It <- Matrix::Diagonal(ntime)
  knots <- c(knots_sp,knots_time);  pord <- c(pord_sp,pord_time)
  bdeg <- c(bdeg_sp,bdeg_time)
  if (!is.null(Xnopar)){
    if (is.null(knots_nopar)) knots_nopar <- rep(10,ncol(Xnopar))
    if (is.null(bdeg_nopar)) bdeg_nopar <- rep(3,ncol(Xnopar))
    if (is.null(pord_nopar)) pord_nopar <- rep(2,ncol(Xnopar))
  }
  if (trace) start <- proc.time()[3]
  # Build design matrices in mixed model
  mat.XZ <- XZ3d(x1=x1,x2=x2,x3=x3,decom=decom,psanova=psanova,
            f1_main=f1_main,f2_main=f2_main,f3_main=f3_main,
            f12_int=f12_int,f13_int=f13_int,f23_int=f23_int,f123_int=f123_int,
            div_sp=div_sp,div_tm=div_tm,
            x1lim=x1lim,x2lim=x2lim,x3lim=x3lim,
            knots=knots,pord=pord,bdeg=bdeg,
            Xpar=Xpar,Xnopar=Xnopar,
            knots_nopar=knots_nopar,bdeg_nopar=bdeg_nopar,
            pord_nopar=pord_nopar)
  X <- mat.XZ$X;  X_spt <- mat.XZ$X_spt;
  Xpar <- mat.XZ$Xpar;  Xnopar <- mat.XZ$Xnopar
  Z <- mat.XZ$Z;  Z_spt <- mat.XZ$Z_spt;  Znopar <- mat.XZ$Znopar
  d1 <- mat.XZ$d1; d2 <- mat.XZ$d2 ; d3 <- mat.XZ$d3
  d_nopar <- mat.XZ$d_nopar
  c1 <- mat.XZ$c1; c2 <- mat.XZ$c2;	c3 <- mat.XZ$c3;
  c_nopar <- mat.XZ$c_nopar
  nvar_par <- mat.XZ$nvar_par; nvar_nopar <- mat.XZ$nvar_nopar
  ncolZ_spt <- ncol(Z_spt)
  rm(mat.XZ)

  # Build vector and matrices for variance components in mixed model
  var_comp <- par_var_comp (la=c(vary_init), np_fixed = ncol(X), pord = pord,
                 c1 = c1, c2 = c2, c3 = c3, d1 = d1, d2 = d2, d3 = d3,
                 c_nopar = c_nopar, pord_nopar = pord_nopar,
                 nvar_nopar = nvar_nopar,
                 psanova = psanova, f1_main = f1_main, f2_main = f2_main,
                 f3_main = f3_main, f12_int = f12_int, f13_int = f13_int,
                 f23_int = f23_int, f123_int = f123_int,
                 tau_init = tau_init, tau_nopar_init = tau_nopar_init)

   # Number of parameters
   np <- var_comp$np; np.eff <- var_comp$np.eff;
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
  D <- diag(c(rep(0,np.eff[1]), rep(1,sum(np.eff[-1]))))
  ed <- rep(0,length(la)-1)
  # add rho and phi to la vector of parameters
  la <- c(la,rho_init,phi_init)
  # Initialise the parameters
  if (is.null(bold)) bold = rep(0,sum(np.eff))
  eta <- X%*%bold[1:np.eff[1]] + Z%*%bold[-(1:np.eff[1])] + offset

  for (iq in 1:maxit) { # Nested loops for SAP and Scores of phi and rho
      if (trace) start1 <- proc.time()[3]
      for (it in 1:maxit) {
	      if (trace) start2 <- proc.time()[3]
	      # CÁLCULO MATRIZ A PARA MODELO PS-SAR
	      rho <- la[length(la)-1]
	      phi <- la[length(la)]
	      if (is.null(Wsp)) {
	          A <- In
	          } else {
	          A <- as(In-rho*Wsp,"CsparseMatrix")
	      }

	      # fast index ntime, slow index nsp
	      A_It_y <- matrix(RH(A,RH(It, array(y,dim=c(ncol(It),ncol(A))))),ncol=1)
	      # Es lo mismo que esta operación (lo he comprobado):
	      #      A_It_y <- matrix(kronecker(A,I_t)%*%y,ncol=1)

	      # Build covariance G matrix for random effects: block diagonal matrix
	      lG <- build_G(la = la, lg = var_comp,
	             nvar_nopar = nvar_nopar, d_nopar = d_nopar,
	             psanova = psanova, f1_main = f1_main, f2_main = f2_main,
	             f3_main = f3_main, f12_int = f12_int, f13_int = f13_int,
	             f23_int = f23_int, f123_int = f123_int)
	      G <- lG$G; Ginv <- lG$Ginv;
	      G.eff <- lG$G.eff; Ginv.eff <- lG$Ginv.eff
	      rm(lG)

	      # X'WX X'WZG
	      # Z'WX Z'WZG
	      sig2u <- la[1]
	      if (phi_fixed & phi_init==0) { # (Omega=I)
	        Xstar <- X
	        Zstar <- Z
	        A_It_ystar <- A_It_y } else {
	        # CONSTRUIMOS MATRIZ OMEGA (SI phi=0 ENTONCES OMEGA=I)
	        call.Omega <- build_Omega_ar1(phi,ntime)
	        Omega <- call.Omega$Omega
	        Omegainv <- call.Omega$Omegainv
	        #Rinv <- as(kronecker(Matrix::Diagonal(nsp), Omegainv),
	        #           "symmetricMatrix")
	        #chol.Rinv <- chol(Rinv)
	        # SE PUEDEN UTILIZAR OPERACIONES RH PARA EVITAR
	        # KRONECKER PRODUCTS
	        #chol.Rinv <- kronecker(Diagonal(nsp),chol(Omegainv))
	        #chol.Omegainv <- chol(Omegainv)
	        chol.Omegainv <- lltA(Omegainv)
	        Xstar <- apply(aperm(RH(In,RH(chol.Omegainv,
	                 array(X,dim=c(ncol(chol.Omegainv),ncol(In),ncol(X))))),
	                 perm=c(2,1,3)),2,rbind)
	        Zstar <- apply(aperm(RH(In,RH(chol.Omegainv,
	                 array(Z,dim=c(ncol(chol.Omegainv),ncol(In),ncol(Z))))),
	                 perm=c(2,1,3)),2,rbind)
	        A_It_ystar <- matrix(RH(In,RH(chol.Omegainv,
	               array(A_It_y,dim=c(ncol(chol.Omegainv),ncol(In))))),ncol=1)

	        # Xstar <- NULL
	        # for(k in 1:ncol(X)){
	        #     Xstar <- cbind(Xstar,matrix(RH(In,RH(chol.Omegainv,
	        #        array(X[,k],dim=c(ncol(chol.Omegainv),ncol(In))))),ncol=1))
	        # }
	        # Zstar <- NULL
	        # for(k in 1:ncol(Z)){
	        #     Zstar <- cbind(Zstar,matrix(RH(In,RH(chol.Omegainv,
	        #        array(Z[,k],dim=c(ncol(chol.Omegainv),ncol(In))))),ncol=1))
	        # }
	        #Xstar <- as.matrix(chol.Rinv %*% X)
	        #Zstar <- as.matrix(chol.Rinv %*% Z)
	        #A_It_ystar <- as.matrix(chol.Rinv) %*% A_It_y
	        rm(call.Omega,Omegainv,chol.Omegainv)
	      }
	      if (is.null(weights)) {
	          w <- as.vector(matrix(1,nrow=length(y))) } else { w <- weights }
	      mat <- construct_matrices(as.matrix(Xstar),
	                 as.matrix(Zstar), as.matrix(A_It_ystar), w, GLAM)
	      C <- construct_block(mat$XtX, t(mat$ZtX*G.eff),
	                           mat$ZtX,
	                           t(mat$ZtZ*G.eff))
	      # ES LA MATRIZ C DE COEFICIENTES DEL SISTEMA (12) EN PAPER SAP
	      # NO SE PUEDE RESOLVER SI SE INCLUYEN LOS ZEROS EN MATRICES X,Z,G
	      Hinv <- try(solve((1/la[1])*C + D,tol = 1e-40))
	      if (class(Hinv) == "try-error")
	        Hinv <- MASS::ginv((1/la[1])*C + D)
	      b <- (1/la[1])*Hinv %*% mat$u
	      b_fixed <- b[1:np.eff[1]]
	      b_random <- G.eff*b[-(1:np.eff[1])]

	      # Compute effective dimensions and variances
	      # Only the diagonal of ZtPZ
	      dZtPZ <- 1/la[1]*apply((t(Hinv[-(1:np.eff[1]),])*mat$ZtXtZ),2,sum)
	      dZtPZ.wide <- rep(0,length(G))
	      b_random.wide <- rep(0,length(G))
	      index.zeros.G <- G==0
	      index <- 1
	      for (ide in 1:length(G)) {
	        if (!index.zeros.G[ide]) {
	          b_random.wide[ide] <- b_random[index]
	          dZtPZ.wide[ide] <- dZtPZ[index]
	          index <- index + 1
	        }
	      }
	      ltau_edf <- update_tau (la = la, lg = var_comp, G = G,
	                  dZtPZ.wide = dZtPZ.wide, b_random.wide = b_random.wide,
	                  psanova = psanova, f1_main = f1_main, f2_main = f2_main,
	                  f3_main = f3_main, f12_int = f12_int, f13_int = f13_int,
	                  f23_int = f23_int, f123_int = f123_int,
	                  nvar_nopar = nvar_nopar, d_nopar = d_nopar,
	                  tau_init = tau_init, tau_fixed = tau_fixed,
	                  tau_nopar_init = tau_nopar_init,
	                  tau_nopar_fixed = tau_nopar_fixed)


	      tau1 <- ltau_edf$tau1; tau2 <- ltau_edf$tau2; tau3 <- ltau_edf$tau3;
	      tau4 <- ltau_edf$tau4; tau5 <- ltau_edf$tau5; tau6 <- ltau_edf$tau6;
	      tau7 <- ltau_edf$tau7; tau8 <- ltau_edf$tau8; tau9 <- ltau_edf$tau9;
	      tau10 <- ltau_edf$tau10; tau11 <- ltau_edf$tau11;
	      tau12 <- ltau_edf$tau12; tau_nopar <- ltau_edf$tau_nopar
	      ed1 <- ltau_edf$ed1; ed2 <- ltau_edf$ed2; ed3 <- ltau_edf$ed3;
	      ed4 <- ltau_edf$ed4; ed5 <- ltau_edf$ed5; ed6 <- ltau_edf$ed6;
	      ed7 <- ltau_edf$ed7; ed8 <- ltau_edf$ed8; ed9 <- ltau_edf$ed9;
	      ed10 <- ltau_edf$ed10; ed11 <- ltau_edf$ed11;
	      ed12 <- ltau_edf$ed12; edf_nopar <- ltau_edf$edf_nopar

	      # OJO: USA LAS VARIABLES STAR PARA CALCULAR SSR Y EDF'S
	      # Regression (Fahrmeir et al.) pp. 180
	      ssr <- mat$yty - t(c(b_fixed, b_random)) %*% (2*mat$u - C %*% b)
	      # New variance components
	      if(!psanova) {
	          lanew <- c(sig2u, tau1, tau2, tau3)
	          edf_tot <- ed1 + ed2 + ed3 + ncol(Xstar)
	      } else {
	          lanew <- c(sig2u, tau1, tau2, tau3, tau4, tau5, tau6, tau7,
	                        tau8, tau9, tau10, tau11, tau12)
	          edf_tot <- ed1 + ed2 + ed3 + ed4 + ed5 + ed6 + ed7 + ed8 + ed9 +
	              ed10 + ed11 + ed12 + ncol(Xstar)
	      }
	      if (!is.null(Xnopar)){
	          lanew <- c(lanew,tau_nopar)
	          edf_tot <- edf_tot + sum(edf_nopar)
	      }
	      if ((!rho_fixed) || !(rho==0))  edf_tot <- edf_tot + 1
	      if ((!phi_fixed) || !(phi==0))  edf_tot <- edf_tot + 1
	      sig2u <- as.numeric((ssr/(length(y) - edf_tot)))
	      # Update first component of la
	      lanew[1] <- sig2u
	      # Update last two components of la
	      lanew <- c(lanew,rho,phi)
	      dla <- mean(abs(la - lanew))
	      la <- lanew

	      if (trace) {
	        #cat(sprintf("%1$3d %2$10.6f", it, dla), '\n')
	        #cat("la ",la,'\n')
	        cat('sig2u ',la[1],'\n')
	        cat('edf ')
	            if(!psanova){
	                cat(sprintf("%8.3f %8.3f %8.3f ",ed1,ed2,ed3),'\n')
	            } else {
	                cat(sprintf("%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f",
	                        ed1,ed2,ed3,ed4,ed5,ed6,ed7,ed8,ed9,
	                        ed10,ed11,ed12),'\n')
	            }
	        }
	      #  convergence check
	      if (dla < thr) break
	    } # end for (it in 1:maxit) #
	      if (trace) {
	      end2 <- proc.time()[3]
	      cat("Timings:\nSAP", (end2-start2), "seconds\n") }

	      # BEGINS RHO-PHI ITERATIONS
# 	      # Se igualan scores a 0 y se resuelven numéricamente las ecuaciones
        # Reparameterization of rho and phi to use unconstrained optimization
#         if(rho>=0) rho_unc <- rho/(1-rho) else rho_unc <- rho/(1+rho)
# 	    if(phi>=0) phi_unc <- phi/(1-phi) else phi_unc <- phi/(1+phi)
#         if (!rho_fixed & !phi_fixed) {
# 	      init.val <- c(rho_unc,phi_unc)
# 	      sc_rho_phi <- score_rho_phi }
# 	    if (!rho_fixed & phi_fixed) {
# 	      init.val <- c(rho_unc)
# 	      sc_rho_phi <- score_rho }
# 	    if (rho_fixed & !phi_fixed) {
# 	      init.val <- c(phi_unc)
# 	      sc_rho_phi <- score_phi }
# 	    if (rho_fixed & phi_fixed) {
# 	      sc_rho_phi <- NULL
# 	      rho <- rho_init
# 	      phi <- phi_init }
#	    if (!is.null(sc_rho_phi)) {

	      # rho_phi_roots <- nlminb(start=init.val,
	      #                         objective=llik_reml,
	      #                         gradient=sc_rho_phi,
	      #                      sig2u=sig2u,nsp=nsp,ntime=ntime,
	      #                     Wsp=Wsp,y=y,X=X,Z=Z,G.eff=G.eff,np.eff=np.eff,
	      #                                       b_fixed=b_fixed,
	      #                     rho_fixed=rho_fixed,
	      #                      phi_fixed=phi_fixed)
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
      if (!rho_fixed || !phi_fixed){
        #cat("rho: ",rho,"\n"); cat("phi: ",phi,"\n")
        rho_phi_optim <- nloptr::nloptr(x0=par_init,
                                   eval_f=llik_reml_fn,
                                   #eval_grad=score_phi,
                                   lb=lb, ub=ub,
                                   opts=list("algorithm"="NLOPT_LN_BOBYQA",
                                              "xtol_rel"=1.0e-1,
                                              "ftol_abs"=1.0e-2,
                                              "maxtime"=100,
                                              "print_level"=1),
                                    sig2u=sig2u,nsp=nsp,ntime=ntime,
                                    Wsp=Wsp,y=y,X=X,Z=Z,G.eff=G.eff,
                                    np.eff=np.eff,
                                    b_fixed=b_fixed,
                                    rho_fixed=rho_fixed,
                                    phi_fixed=phi_fixed)

        rho <- rho_phi_optim$solution[1]
        phi <- rho_phi_optim$solution[2]

      }

          # rho_phi_roots <- rootSolve::multiroot(sc_rho_phi,start=init.val,
	        #                    sig2u=sig2u,nsp=nsp,ntime=ntime,
	        #                    Wsp=Wsp,y=y,X=X,Z=Z,G.eff=G.eff,np.eff=np.eff,
	        #                    b_fixed=b_fixed,trace=trace)


# 	      if (class(rho_phi_roots)=="try-error"){
# 	          if (!rho_fixed & !phi_fixed) {
# 	              rho <- init.val[1]
# 	              phi <- init.val[2]
#               }
# 	          if (!rho_fixed & phi_fixed) rho <- init.val
# 	          if (rho_fixed & !phi_fixed) phi <- init.val
# 	          la[length(la)-1] <- rho
# 	          la[length(la)] <- phi
# 	          next
# 	      }
	      # Undo the reparameterization
	      # if (!rho_fixed && !phi_fixed) {
	      #     rho <- rho_phi_roots$root[1] / (1 + abs(rho_phi_roots$root[1]))
	      #     phi <- rho_phi_roots$root[2] / (1 + abs(rho_phi_roots$root[2]))
	      # } else if (!rho_fixed && phi_fixed) {
	      #     rho <- rho_phi_roots$root[1] / (1 + abs(rho_phi_roots$root[1]))
	      # } else if (rho_fixed && !phi_fixed) {
	      #     phi <- rho_phi_roots$root[1] / (1 + abs(rho_phi_roots$root[1]))
	      #
	      # }
#	     } # end of if (!is.null(sc_rho_phi))
	    # Check Convergence in rho and phi parameters
	    rho_old <- la[length(la)-1]
	    phi_old <- la[length(la)]
	    drhophi <- mean(abs(c(rho_old,phi_old)-c(rho,phi)))
	    la[length(la)-1] <- rho
	    la[length(la)] <- phi
	    if (trace) {
	        if (!rho_fixed) cat("rho ",rho,'\n')
	        if (!phi_fixed) cat("phi",phi,'\n')
	     }
	     if (drhophi < thr) break
	  } # End loop for SAP and Scores of phi and rho

  if (trace) { end1 <- proc.time()[3]
			cat("Timings:\nSAP-Scores rho-phi", (end1-start1), "seconds\n") }
  eta <- Xstar%*%b_fixed + Zstar%*%b_random + offset

  if (is.null(Xnopar)) { edf_nopar <- NULL; tau_nopar<-NULL }

#  FINAL ESTIMATES OF PARAMETERS
  sig2u <- la[1]
  if(!psanova){
      tau1 <- la[2]; tau2 <- la[3]; tau3 <- la[4]
  } else {
      tau_f1_main <- la[2]; tau_f2_main <- la[3]; tau_f3_main <- la[4]
      tau_f12 <- la[5:6]; tau_f13 <- la[7:8]; tau_f23 <- la[9:10];
      tau_f123 <- la[11:13]
  }
  rho <- la[length(la)-1]
  phi <- la[length(la)]

 # Valor de log.lik y log.lik.reml en el óptimo
  foptim <- llik_reml_var(c(rho,phi),
                         sig2u=sig2u,nsp=nsp,ntime=ntime,
                         Wsp=Wsp,y=y,X=X,Z=Z,G.eff=G.eff,
                         np.eff=np.eff,b_fixed=b_fixed,
                         rho_fixed=rho_fixed,phi_fixed=phi_fixed)
  llik_reml <- foptim$llik_reml
  llik <- foptim$llik
  sd_rho_an <- foptim$sd_rho
  sd_phi_an <- foptim$sd_phi
  cov_rho_phi_an <- foptim$cov_rho_phi
  rm(foptim)
  hessian_optim <- numDeriv::hessian(llik_reml_fn,c(rho,phi),
                                     sig2u=sig2u,nsp=nsp,ntime=ntime,
                                     Wsp=Wsp,y=y,X=X,Z=Z,G.eff=G.eff,
                                     np.eff=np.eff,
                                     b_fixed=b_fixed,
                                     rho_fixed=rho_fixed,
                                     phi_fixed=phi_fixed)
  var_rho_phi_num <- solve(hessian_optim)

  cat("llik_reml: ",llik_reml,"\n")
  cat("llik: ",llik,"\n")
  cat("sd_rho_an: ",sd_rho_an,"\n")
  cat("sd_phi_an: ",sd_phi_an,"\n")
  cat("cov_rho_phi_an: ",cov_rho_phi_an,"\n")
  cat("sd_rho_num: ",sqrt(var_rho_phi_num[1,1]),"\n")
  cat("sd_phi_num: ",sqrt(var_rho_phi_num[2,2]),"\n")
  cat("cov_rho_phi_num: ",var_rho_phi_num[1,2],"\n")

     # CONTINUAR AQUÍ
		# Asignación Efectos Fijos y Aleatorios según tipo variables
		b_fixed_par <- NULL; b_fixed_nopar <- NULL; b_random_nopar <- NULL
		if (nvar_par>0) b_fixed_par <- b_fixed[1:ncol(Xpar)]
    if (nvar_nopar>0){
      b_fixed_nopar <- b_fixed[(ncol(X)-ncol(Xnopar)+1):ncol(X)]
      b_random_spt <- b_random[1:(ncol(Z)-ncol(Znopar))]
      b_random_nopar <- b_random[(ncol(Z)-ncol(Znopar)+1):ncol(Z)]
      if (nvar_par>0)
      {
        b_fixed_spt <- b_fixed[(ncol(Xpar)+1):(ncol(X)-ncol(Xnopar))]
      } else { b_fixed_spt <- b_fixed[1:(ncol(X)-ncol(Xnopar))] }
    }	else {
      b_random_spt <- b_random
      if (nvar_par>0) { b_fixed_spt <- b_fixed[(ncol(Xpar)+1):ncol(X)]
      } else { b_fixed_spt <- b_fixed[1:ncol(X)]}
    }

		if (post_estim==TRUE)
		{
		  # CÁLCULO MATRIZ COVARIANZAS EFECTOS FIJOS Y ALEATORIOS
		  # pp.375 Fahrmeir et al.
		  # Se reescala A multiplicándola por sig2u toda la matriz
		  # Matriz Covarianzas Bayesiana. OJO: Tiene en cuenta autocorrelación...
		  if (phi_fixed && phi_init == 0) { # (Omega=I)
		    #Omega <- It
		    Omegainv <- It
		    } else {
		      # CONSTRUIMOS MATRIZ OMEGA (SI phi=0 ENTONCES OMEGA=I)
		      call.Omega <- build_Omega_ar1(phi,ntime)
		      #Omega <- call.Omega$Omega
		      Omegainv <- Matrix::Matrix(call.Omega$Omegainv)
		      rm(call.Omega)
		    }
		  # Rinv <- as(kronecker(In,Omegainv),"symmetricMatrix")
		  #chol.Rinv <- chol(Rinv)
		  Rinv <- kronecker(In,Omegainv) # OJO: SPARSE (CAMBIAR CÓDIGO)
		  chol.Rinv <- Matrix::Matrix(lltA(as(Rinv,"matrix")))
		  Xstar <- chol.Rinv %*% X
		  Zstar <- chol.Rinv %*% Z
		  #Astar <- as.matrix(chol.Rinv)%*%kronecker(A,It)
		  A_cov1 <- rbind(cbind(Matrix::crossprod(Xstar),
		                        Matrix::t(Xstar) %*% Zstar),
		                  cbind(Matrix::t(Zstar) %*% Xstar,
		                        Matrix::crossprod(Zstar) +
		                          sig2u*Matrix::Matrix(diag(Ginv.eff)) ))
		  cov1_eff <- try(sig2u*Matrix::solve(A_cov1))
		  if (class(cov1_eff) == "try-error") {
		    cov1_eff <- Matrix::Matrix(sig2u*MASS::ginv(as(A_cov1,"matrix")))
		    cov1_eff.ginv <- TRUE
		  } else {cov1_eff.ginv <- FALSE}
		  var1_b_fixed <- Matrix::diag(cov1_eff)[1:(ncol(X))]
		  var1_b_random <- Matrix::diag(cov1_eff)[(ncol(X)+1):(ncol(X)+ncol(Z))]
		  var1_b_fixed_par <- NULL; cov1_b_fixed_par <- NULL
		  var1_b_fixed_nopar <- NULL; var1_b_random_nopar <- NULL
		  if (nvar_par>0)
		  {
		    var1_b_fixed_par <- Matrix::diag(cov1_eff)[1:ncol(Xpar)]
		    cov1_b_fixed_par <- cov1_eff[1:ncol(Xpar),1:ncol(Xpar)]
		  }

		  if (nvar_nopar>0)
		  {
		    var1_b_random_spt <- Matrix::diag(cov1_eff)[(ncol(X)+1):
		                                          (ncol(X)+ncol(Z)-ncol(Znopar))]
		    var1_b_random_nopar <- Matrix::diag(cov1_eff)[(ncol(X)+length(var1_b_random_spt)+1):
		                                            (ncol(X)+ncol(Z))]
		    if (nvar_par>0){
		      var1_b_fixed_spt <- Matrix::diag(cov1_eff)[(ncol(Xpar)+1):
		                                           (ncol(X)-ncol(Xnopar))]
		      var1_b_fixed_nopar <- Matrix::diag(cov1_eff)[(ncol(Xpar) +
		                                          length(b_fixed_spt)+1):ncol(X)]
		    } else {
		      var1_b_fixed_spt <- Matrix::diag(cov1_eff)[1:ncol(X_spt)]
		      var1_b_fixed_nopar <- Matrix::diag(cov1_eff)[(ncol(X_spt)+1):ncol(X)]}
		  } else {
		    var1_b_random_spt <- Matrix::diag(cov1_eff)[(ncol(X)+1):
		                                          (ncol(X)+ncol(Z))]
		    if (nvar_par>0){
		      var1_b_fixed_spt <- Matrix::diag(cov1_eff)[(ncol(Xpar)+1):ncol(X)]
		    } else {
		      var1_b_fixed_spt <- Matrix::diag(cov1_eff)[1:ncol(X)]}
		  }
		  # Matriz Covarianzas Frequentist tipo Sandwich (Algo menor las varianzas)
		  # C.Rinv.C <- (1/sig2u)*rbind(cbind(t(X)%*%kronecker(Diagonal(nsp),Omegainv)%*%X,
		  #                                   t(X)%*%kronecker(Diagonal(nsp),Omegainv)%*%Z),
		  #                        cbind(t(Z)%*%kronecker(Diagonal(nsp),Omegainv)%*%X,
		  #                              t(Z)%*%kronecker(Diagonal(nsp),Omegainv)%*%Z))
		  C.Rinv.C <- (1/sig2u)*rbind(cbind(Matrix::crossprod(Xstar),
		                                    Matrix::t(Xstar) %*% Zstar),
		                              cbind(Matrix::t(Zstar) %*% Xstar,
		                                    Matrix::crossprod(Zstar)))
		  cov2.eff <- cov1_eff %*% C.Rinv.C %*% cov1_eff
		  var2_b_fixed <- Matrix::diag(cov2.eff)[1:ncol(X)]
		  var2_b_random <- Matrix::diag(cov2.eff)[(ncol(X)+1):ncol(X)+ncol(Z)]
		  var2_b_fixed_par <- NULL; cov2.b_fixed_par <- NULL
		  var2_b_fixed_nopar <- NULL; var2_b_random_nopar <- NULL

		  if (nvar_par>0)
		  {
		    var2_b_fixed_par <- Matrix::diag(cov2.eff)[1:ncol(Xpar)]
		    cov2.b_fixed_par <- cov2.eff[1:ncol(Xpar),1:ncol(Xpar)]
		  }

		  if (nvar_nopar>0)
		  {
		    var2_b_random_spt <- Matrix::diag(cov1_eff)[(ncol(X)+1):
		                                          (ncol(X)+ncol(Z)-ncol(Znopar))]
		    var2_b_random_nopar <- Matrix::diag(cov1_eff)[(ncol(X)+length(var2_b_random_spt)+1):
		                                            (ncol(X)+ncol(Z))]
		    if (nvar_par>0){
		      var2_b_fixed_spt <- Matrix::diag(cov1_eff)[(ncol(Xpar)+1):
		                                   (ncol(X)-ncol(Xnopar))]
		      var2_b_fixed_nopar <- Matrix::diag(cov1_eff)[(ncol(Xpar)+length(b_fixed_spt)+1):
		                                             ncol(X)]
		    } else {
		      var2_b_fixed_spt <- Matrix::diag(cov1_eff)[1:ncol(X_spt)]
		      var2_b_fixed_nopar <- Matrix::diag(cov1_eff)[(ncol(X_spt)+1):ncol(X)]
		    }
		  } else {
		          var2_b_random_spt <- Matrix::diag(cov1_eff)[(ncol(X)+1):
		                                          (ncol(X)+ncol(Z))]
		          if (nvar_par>0){
		              var2_b_fixed_spt <- Matrix::diag(cov1_eff)[(ncol(Xpar)+1):
		                                                   ncol(X)]
		          } else {
		              var2_b_fixed_spt <- Matrix::diag(cov1_eff)[1:ncol(X)]
		          }
		  }

		  # CALCULAR VALORES ESTIMADOS Y RESIDUOS
		  fit_Ay <- X %*% b_fixed + Z %*% b_random + offset
		  fit <- try((Matrix::solve(A) %x% It) %*% fit_Ay)
      if (class(fit) == "try-error") {
        fit <- (MASS::ginv(as(A,"matrix")) %x% It) %*% eta
      }
      # OJO: pp. 375 Fahrmeir et al.
		  # fit = mH*y
      # mH = kron(Ainv,It)*(X Z)*VAR(beta,alpha)*t(X Z)*Rinv*kron(A,It)
		  mH <- (1/sig2u)*(Matrix::solve(A) %x% It) %*% cbind(X,Z) %*%
		    (cov1_eff %*% (Matrix::t(cbind(X,Z))%*%Rinv %*% (A %x% It)))

		  var_fit_Ay <- cbind(X,Z) %*% cov1_eff %*% Matrix::t(cbind(X,Z))
		  var_fit <- try((solve(A) %x% It) %*% var_fit_Ay %*% (t(solve(A))%x%It))
		  if (class(var_fit) == "try-error") {
		    (MASS::ginv(as(A,"matrix")) %x% It) %*% var_fit_Ay %*%
		       (t(MASS::ginv(as(A,"matrix"))) %x% It)
      }
		  # resids N(0,sigma_u^2*Omega) and resids.norm N(0,sigma_eps*I)
		  resids <- as.vector((A %x% It)%*% y - fit_Ay)
		  resids.norm <- sqrt(1-phi^2)*(chol.Rinv %*% resids)
		  # Compute AIC y BIC based on loglik functions (Fahrmeir, pp. 664 and 677)
		  aic <- -2*llik + 2*edf_tot
		  bic <- -2*llik + log(length(y))*edf_tot
		  # CÁLCULO AIC and BIC based on ssr: LIBRO FAHRMEIR pp. 564 (STAR MODELS WITH GAUSSIAN ERRORS)
		  # VIP: SSR CON RESIDUOS NORMALIZADOS
		  ssr <- sum(resids^2)
		  ssr.norm <- sum(resids.norm^2)
		  aic.sig2 <- length(y)*log(ssr.norm/length(y))+2*(edf_tot+1)
		  bic.sig2 <- length(y)*log(ssr.norm/length(y))+log(length(y))*(edf_tot+1)

		  # COMPUTE FITS EFFECTS BY VARIABLE
		  fit_cov_par <- NULL;	fit_cov_nopar <- NULL
		  fit_cov_nopar_fixed <- NULL;	fit_cov_nopar_random <- NULL
		  var_fit_cov_nopar <- NULL; k.Znopar <- NULL; edf_nopar_tot <- NULL
		  if (nvar_par>0)
		  {
		    for (i in 1:nvar_par)
		    {
		      fit_cov_par <- cbind(fit_cov_par,Xpar[,i]*b_fixed_par[i])
		    }
		  }
		  if (nvar_nopar>0)
		  {
		    edf_nopar_tot <- rep(0,nvar_nopar)
		    k.Znopar <- rep(0,nvar_nopar)
		    for (i in 1:nvar_nopar)
		    {
		      fit_cov_nopar_fixed <- cbind(fit_cov_nopar_fixed,
		                                   Xnopar[,i]*b_fixed_nopar[i])
		      var_np<-Xnopar[,i]
		      # Matriz selección para varianza f(x_j) Libro Ruppert et al, pp. 175
		      E_var_np <- matrix(0,ncol=ncol(cbind(X,Z)),nrow=ncol(cbind(X,Z)))
		      vindex.var_np_fixed <- (ncol(X)-ncol(Xnopar)+i)
		      E_var_np[vindex.var_np_fixed,vindex.var_np_fixed] <- 1

		      MM.var_np <- MM_basis(var_np,min(var_np)-0.01,max(var_np)+0.01,
		                            knots_nopar[i],bdeg_nopar[i],pord_nopar[i])
		      X.var_np <- MM.var_np$X; Z.var_np <- MM.var_np$Z
		      k.Znopar[i] <- ncol(Z.var_np)
		      if (i==1){
		        fit_cov_nopar_random <- cbind(fit_cov_nopar_random,
		                                      Z.var_np%*%b_random_nopar[1:k.Znopar[1]])
		        vindex.var_np_random <- c((ncol(X)+ncolZ_spt+1):
		                                    (ncol(X)+ncolZ_spt+k.Znopar[1]))
		      } else {
		        fit_cov_nopar_random <- cbind(fit_cov_nopar_random,
		                                      Z.var_np%*%b_random_nopar[(sum(k.Znopar[1:(i-1)])+1):
		                                                                  sum(k.Znopar[1:i])])
		        vindex.var_np_random<- c((ncol(X)+ncolZ_spt+
		                                    sum(k.Znopar[1:(i-1)])+1):
		                                   (ncol(X)+ncolZ_spt+sum(k.Znopar[1:i])))
		      }
		      # Se cambia la diagonal de la matriz de Selección
		      for (j in 1:length(vindex.var_np_random))
		      {
		        E_var_np[vindex.var_np_random[j],vindex.var_np_random[j]] <- 1
		      }
		      # REPASAR
		      mH_var_np <- (1/sig2u)*(Matrix::solve(A) %x% It) %*% cbind(X,Z) %*%
		              E_var_np %*% (cov1_eff %*%
		                            (Matrix::t(cbind(X,Z)) %*% Rinv %*% (A %x% It)))
		      edf_nopar_tot[i] <- sum(Matrix::diag(mH_var_np))
		      var_fit_np <- rep(0,length(fit))
		      for (j in 1:length(fit))
		      {
		        var_fit_np[j] <- sig2u*sum(mH_var_np[j,]^2)
		      }
		      var_fit_cov_nopar <- cbind(var_fit_cov_nopar,var_fit_np)
		    }
		    fit_cov_nopar <- fit_cov_nopar_fixed+fit_cov_nopar_random
		  }

		  # COMPUTE SPATIO-TEMPORAL TREND.
		  spt_trend_fixed <- as.vector(X_spt %*% b_fixed_spt)
		  spt_trend_random <- as.vector(Z_spt %*% b_random_spt)
		  spt_trend <- spt_trend_fixed+spt_trend_random
		  E_spt_trend <- matrix(0,ncol=ncol(cbind(X,Z)),nrow=ncol(cbind(X,Z)))
		  vindex_spt_trend_fixed <- c((nvar_par+1):(ncol(X)-nvar_nopar))
		  for (i in 1:length(vindex_spt_trend_fixed))
		  {
		    E_spt_trend[vindex_spt_trend_fixed[i],
		                vindex_spt_trend_fixed[i]] <- 1
		  }
		  vindex_spt_trend_random <- c((ncol(X)+1):(ncol(X)+ncolZ_spt))
		  for (i in 1:length(vindex_spt_trend_random))
		  {
		    E_spt_trend[vindex_spt_trend_random[i],
		                vindex_spt_trend_random[i]] <- 1
		  }
		  mH_spt_trend <- (1/sig2u)*(Matrix::solve(A) %x% It) %*%
		                   cbind(X,Z) %*% E_spt_trend %*%
		    (cov1_eff %*% (Matrix::t(cbind(X,Z)) %*% Rinv %*% (A %x% It)))
		  edf_spt_tot <- sum(Matrix::diag(mH_spt_trend))
		  var_spt_trend <- rep(0,length(spt_trend))
		  for (j in 1:length(spt_trend))
		  {
		    var_spt_trend[j] <- sig2u*sum(mH_spt_trend[j,]^2)
		  }
		  spt_trend_fixed <- matrix(spt_trend_fixed,ncol=nsp,byrow=FALSE)
		  spt_trend_random <- matrix(spt_trend_random,ncol=nsp,byrow=FALSE)
		  spt_trend <- matrix(spt_trend,ncol=nsp,byrow=FALSE)
		  var_spt_trend <- matrix(var_spt_trend,ncol=nsp,byrow=FALSE)

		} else {
		     cov1_b_fixed_par=NULL;cov2.b_fixed_par=NULL
		        var1_b_fixed=NULL;var2_b_fixed=NULL
		        var1_b_random=NULL;var2_b_random=NULL
		        var1_b_fixed_spt=NULL;var2_b_fixed_spt=NULL
		        var1_b_random_spt=NULL;var2_b_random_spt=NULL
		        var1_b_fixed_par=NULL;var2_b_fixed_par=NULL
		        var1_b_fixed_nopar=NULL;var2_b_fixed_nopar=NULL
		        var1_b_random_nopar=NULL;var2_b_random_nopar=NULL
		        fit=NULL;resids=NULL;resids.norm=NULL;var_fit=NULL
		        fit_cov_par=NULL;fit_cov_nopar=NULL
		        fit_cov_nopar_fixed=NULL;fit_cov_nopar_random=NULL
		        var_fit_cov_nopar=NULL;spt_trend=NULL
		        spt_trend_fixed=NULL;spt_trend_random=NULL
		        var_spt_trend=NULL
		        k.Znopar=NULL;edf_nopar=NULL;edf_nopar_tot=NULL
		        edf_spt_tot=NULL;llik_reml=NULL
		      } # end if (post_estim==TRUE)
		  res <- list (call = match.call(),
               knots=knots,bdeg=bdeg,pord=pord,
	             knots_nopar=knots_nopar,bdeg_nopar=bdeg_nopar,pord_nopar=pord_nopar,
	             k.Znopar=k.Znopar,
	             # Parámetros estimados
	             sig2u=sig2u,rho=rho,phi=phi,
                 sd_rho=sd_rho_an,sd_phi=sd_phi_an,
                 cov_rho_phi=cov_rho_phi_an,
                 tau1=tau1,tau2=tau2,tau3=tau3,ed1=ed1,ed2=ed2,ed3=ed3,
	             tau_f1_main=tau_f1_main,tau_f2_main=tau_f2_main,
	             tau_f3_main=tau_f3_main,tau_f12=tau_f12,tau_f13=tau_f13,
	             tau_f23=tau_f23,tau_f123=tau_f123,tau_nopar=tau_nopar,
	             edf_f1_main=ed1,edf_f2_main=ed2,edf_f3_main=ed3,
	             edf_f12=c(ed4,ed5),edf_f13=c(ed6,ed7),edf_f23=c(ed8,ed9),
	             edf_f123=c(ed10,ed11,ed12),
	             edf_nopar=edf_nopar,edf_tot=edf_tot,
	             edf_nopar_tot=edf_nopar_tot,
	             edf_spt_tot=edf_spt_tot,
	             ssr=ssr,ssr.norm=ssr.norm,aic=aic,bic=bic,
               aic.sig2=aic.sig2,bic.sig2=bic.sig2,
               llik_reml=llik_reml,llik=llik,
               #V=V,Vinv=Vinv,P=P,Z=Z,G=G,X=X,
               #A=A,Omega=Omega,Omegainv=Omegainv,
               b_fixed=b_fixed,b_random=b_random,
	             b_fixed_spt=b_fixed_spt,b_random_spt=b_random_spt,
	             b_fixed_par=b_fixed_par,b_fixed_nopar=b_fixed_nopar,
	             b_random_nopar=b_random_nopar,
	             #cov1_eff=cov1_eff,cov2.eff=cov2.eff,
	             var1_b_fixed=var1_b_fixed,var2_b_fixed=var2_b_fixed,
	             var1_b_random=var1_b_random,var2_b_random=var2_b_random,
	             var1_b_fixed_spt=var1_b_fixed_spt,
	             var2_b_fixed_spt=var2_b_fixed_spt,
	             var1_b_random_spt=var1_b_random_spt,
	             var2_b_random_spt=var2_b_random_spt,
	             var1_b_fixed_par=var1_b_fixed_par,
	             cov1_b_fixed_par=cov1_b_fixed_par,
	             var2_b_fixed_par=var2_b_fixed_par,
	             cov2.b_fixed_par=cov2.b_fixed_par,
	             var1_b_fixed_nopar=var1_b_fixed_nopar,
	             var2_b_fixed_nopar=var2_b_fixed_nopar,
	             var1_b_random_nopar=var1_b_random_nopar,
	             var2_b_random_nopar=var2_b_random_nopar,
	             # fit, resids y var(fit)
               fit=fit,fit_Ay=fit_Ay,resids=resids,resids.norm=resids.norm,
               var_fit=var_fit,var_fit_Ay=var_fit_Ay,fit_cov_par=fit_cov_par,
	             fit_cov_nopar=fit_cov_nopar,
	             fit_cov_nopar_fixed=fit_cov_nopar_fixed,
	             fit_cov_nopar_random=fit_cov_nopar_random,
	             var_fit_cov_nopar=var_fit_cov_nopar,
	             # tendencias espacio-temporales
	             spt_trend=spt_trend,
	             spt_trend_fixed=spt_trend_fixed,spt_trend_random=spt_trend_random,
	             var_spt_trend=var_spt_trend,eta=eta,maxit=it)
	res
}

