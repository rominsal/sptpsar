fitpsar2d <-
function(y,Xsp,Xpar=NULL,Xnopar=NULL,
    knots=c(10,10),pord=c(2,2),bdeg=c(3,3),
    x1lim=NULL,x2lim=NULL,
    vary_init=var(y), tau_init = rep(1,2), tau_fixed = rep(FALSE,2),
    rho_init = 0, rho_fixed = FALSE,
    onlygam=FALSE,trace=FALSE,thr=1e-4,maxit=200,
    offset=0,bold=NULL,weights=NULL,Wsp=NULL,post_estim=TRUE,
    knots_nopar=NULL,bdeg_nopar=NULL,pord_nopar=NULL) {


#  if(!exists("spline.des")) require(splines)
#  if(!exists("ginv")) require(MASS)
#  if(!exists("multiroot")) require(rootSolve)

  nsp <- length(y)
  x1 <- Xsp[,1]
  x2 <- Xsp[,2]
  if(trace) start <- proc.time()[3]
  if(!onlygam)
  {
    if(is.null(x1lim))  x1lim <- c(min(x1)-0.01, max(x1)+0.01)
    if(is.null(x2lim))  x2lim <- c(min(x2)-0.01, max(x2)+0.01)

    MM1 <- MM_basis(x1, x1lim[1], x1lim[2], knots[1], bdeg[1], pord[1])
    MM2 <- MM_basis(x2, x2lim[1], x2lim[2], knots[2], bdeg[2], pord[2])
    X1 <- MM1$X; Z1<-MM1$Z; d1<-MM1$d; B1<-MM1$B
    X2 <- MM2$X; Z2<-MM2$Z;	d2<-MM2$d; B2<-MM2$B
    c1 <- ncol(B1); c2 <- ncol(B2); rm(MM1,MM2,B1,B2)
    g1u <- rep(1,pord[2])%x%d1
    g2u <- d2%x%rep(1,pord[1])
    g1b <- rep(1,c2-pord[2])%x%d1
    g2b <- d2%x%rep(1,c1-pord[1])

    X_spt <- Rten2(X2, X1)
    Z_spt <- cbind(Rten2(X2, Z1), Rten2(Z2, X1),Rten2(Z2, Z1))
    rm(X1,X2,Z1,Z2)
    ncolX_spt<-ncol(X_spt)
    ncolZ_spt<-ncol(Z_spt)
  } else {
    X_spt <- 1; Z_spt <- NULL
    ncolX_spt <- 1; ncolZ_spt <- 0
  }
	X <- X_spt
	Z <- Z_spt

	if (is.null(Xpar))
	{
	  X_par <- NULL
	  nvar_par <- 0
	} else {
	  X_par <- Xpar
	  nvar_par <- ncol(X_par)
	  X <- cbind(X_par,X)
	}

	# Number of parameters in each part
	if(!onlygam)
	{
	  np <-  c(prod(pord)+nvar_par, (c1-pord[1])*pord[2],
	           (c2-pord[2])*pord[1],(c1-pord[1])*(c2-pord[2]))
	} else { np <- ncol(X) }

	if(!is.null(Xnopar))
	{
	  nvar_nopar <- ncol(Xnopar)
	  X_nopar <- NULL; Z_nopar <- NULL; d_nopar <- NULL; c_nopar <- NULL
	  if(is.null(knots_nopar)) knots_nopar <- rep(10,nvar_nopar)
	  if(is.null(bdeg_nopar)) bdeg_nopar <- rep(3,nvar_nopar)
	  if(is.null(pord_nopar)) pord_nopar <- rep(2,nvar_nopar)
	  for(i in 1:nvar_nopar)
	  {
	    var_np<-Xnopar[,i]
	    MM.var <- MM_basis(var_np,min(var_np)-0.01,max(var_np)+0.01,
	                       knots_nopar[i],bdeg_nopar[i],pord_nopar[i])
	    X.var <- MM.var$X; Z.var <- MM.var$Z; d.var <- MM.var$d; B.var <- MM.var$B
	    X_nopar <- cbind(X_nopar,X.var[,c(-1)]) # Quita columna intercepto
	    Z_nopar <- cbind(Z_nopar,Z.var)
	    d_nopar <- cbind(d_nopar,d.var)
	    c_nopar <- c(c_nopar,ncol(B.var))
	    rm(MM.var,B.var)
	  }
	  X <- cbind(X,X_nopar)
	  Z <- cbind(Z,Z_nopar)
	  np[1] <- np[1]+ncol(X_nopar)
	  np <- c(np,c_nopar-pord_nopar)
	} else {
	  nvar_nopar <- 0
	  X_nopar <- NULL; Z_nopar <- NULL; d_nopar <- NULL; c_nopar <- NULL
	}
	# 0	0
	# 0 I
	D <- diag(c(rep(0,np[1]), rep(1,sum(np[-1]))))
	if(!onlygam)
	{
	  G1inv.n <- c(g1u, rep(0, np[3]), g1b)
	  G2inv.n <- c(rep(0,np[2]), g2u, g2b)
	# Initialize variance components
	  la <- c(vary_init,tau_init[1:2],rho_init)
	  ed <- c(0,0)
	} else {
	  la <- c(vary_init,rho_init)
	  ed <- NULL
	}

  # AÑADE PARÁMETROS SMOOTHING COVARIABLES NO PARAMÉTRICAS
  if(!is.null(Xnopar))
  {
    la <- c(la,rep(1,nvar_nopar))
    ed <- c(ed,rep(0,nvar_nopar))
  }

	# Initialise the parameters
	if(is.null(bold)) bold <- rep(0,sum(np))
	if(!is.null(Z)){ eta <- X%*%bold[1:np[1]] + Z%*%bold[-(1:np[1])] + offset } else {
	  eta <- X%*%bold[1:np[1]] + offset }

	for (iq in 1:maxit) { # Nested loops for SAP and Scores of rho
	  w <- rep(1,length(y))
		if(trace) start1 <- proc.time()[3]
		for (it in 1:maxit) { # Loop for SAP
		  if(trace) start2 <- proc.time()[3]
      # CÁLCULO MATRIZ A PARA MODELO PS-SAR
      if(!onlygam) rho <- la[4] else rho <- la[2]
      if(!is.null(Wsp))
      {
          A <- Matrix::Diagonal(nsp)-rho*Wsp
      } else {
          A <- Matrix::Diagonal(nsp)
      }
      Ay <- A%*%y
      # CASO GAUSSIANO
      mat <- construct_matrices(as.matrix(X),as.matrix(Z),
                                as.matrix(Ay),w=1,GLAM=FALSE)
			# Build penalty matrix: block diagonal matrix
      if(!onlygam)
      {
        Ginv <- c((1/la[2])*g1u, (1/la[3])*g2u, (1/la[2])*g1b + (1/la[3])*g2b)
      } else { Ginv <- NULL }
	 if(!is.null(Xnopar))
	 {
			  for(k in 1:nvar_nopar)
			  {
			    if(!onlygam){ Ginv <- c(Ginv,(1/la[4+k])*d_nopar[,k]) } else {
			                  Ginv <- c(Ginv,(1/la[2+k])*d_nopar[,k]) }
			  }
	}
	G <- 1/Ginv
    # X'WX X'WZG
	# Z'WX Z'WZG
	C <- construct_block(mat$XtX, t(mat$ZtX*G), mat$ZtX, t(mat$ZtZ*G))
    # ES LA MATRIZ C DE COEFICIENTES DEL SISTEMA (12) EN PAPER SAP
			Hinv <- try(solve((1/la[1])*C + D))
			if(class(Hinv) == "try-error") {
				Hinv <- MASS::ginv((1/la[1])*C + D)
			}
			b <- (1/la[1])*Hinv%*%mat$u
			b_fixed <- b[1:np[1]]
			b_random <- G*b[-(1:np[1])]
      if(nvar_par==0) {b_fixed_par <- NULL} else {b_fixed_par <- b_fixed[1:nvar_par]}
      sig2u<-la[1]

			# Compute effective dimensions and variances
			# Only the diagonal of ZtPZ
			dZtPZ <- 1/la[1]*apply((t(Hinv[-(1:np[1]),])*mat$ZtXtZ),2,sum)

			if(!onlygam)
			{
			  # Tau 1
			  if(is.null(Xnopar)) {G1inv.d <- (1/la[2])*G1inv.n} else{
			    G1inv.d <- c((1/la[2])*G1inv.n,rep(0,sum(np[-c(1:4)])))}
			  ed1 <- sum(dZtPZ*(G1inv.d*G^2))
			  ed1 <- ifelse(ed1 == 0, 1e-50,ed1)
			  if(!tau_fixed[1])
			  {
			    tau1 <- ifelse(is.null(Xnopar),sum(b_random^2*G1inv.n)/ed1,
			                   sum(b_random^2*c(G1inv.n,rep(0,sum(np[-c(1:4)]))))/ed1)
			    tau1 <- ifelse(tau1 == 0, 1e-50, tau1)
			  } else {tau1 <- tau_init[1]}

			  # Tau 2
			  if(is.null(Xnopar)) {G2inv.d <- (1/la[3])*G2inv.n} else{
			    G2inv.d <- c((1/la[3])*G2inv.n,rep(0,sum(np[-c(1:4)])))}
			  ed2 <- sum(dZtPZ*(G2inv.d*G^2))
			  ed2 <- ifelse(ed2 == 0, 1e-50, ed2)
			  if(!tau_fixed[2])
			  {
			    tau2 <- ifelse(is.null(Xnopar),sum(b_random^2*G2inv.n)/ed2,
			                   sum(b_random^2*c(G2inv.n,rep(0,sum(np[-c(1:4)]))))/ed2)
			    tau2 <- ifelse(tau2 == 0, 1e-50, tau2)
			  } else {tau2 <- tau_init[2]}
			}

			if(!is.null(Xnopar))
			{
			  Ginv.d_nopar <- matrix(0,nrow=sum(np[2:length(np)]),ncol=ncol(Xnopar))
			  tau_nopar <- rep(0,nvar_nopar);  edf_nopar <- rep(0,nvar_nopar)
			  for (k in 1:nvar_nopar)
			  {
			    if(!onlygam)
			    {
			      Ginv.d_nopar[(sum(np[2:(4+k-1)])+1):sum(np[2:(4+k)]),k] <- (1/la[4+k])*d_nopar[,k]
			    } else {
			      cumnp_nopar <- cumsum(np[-1])
			      if(k==1)
			      {Ginv.d_nopar[(1:cumnp_nopar[1]),k] <- (1/la[2+k])*d_nopar[,k]} else {
			       Ginv.d_nopar[((cumnp_nopar[k-1]+1):cumnp_nopar[k]),k] <- (1/la[2+k])*d_nopar[,k]}
			    }
			    edf_nopar[k] <- sum(dZtPZ*(Ginv.d_nopar[,k]*G^2))
			    edf_nopar[k] <- ifelse(edf_nopar[k] == 0, 1e-50, edf_nopar[k])
			    tau_nopar[k] <- sum(b_random^2*Ginv.d_nopar[,k])/edf_nopar[k]
			    tau_nopar[k] <- ifelse(tau_nopar[k] == 0, 1e-50, tau_nopar[k])
			  }
			}
			P.Hv <- (1/sig2u)*(Matrix::Diagonal(nsp)-
			           (1/sig2u)*( cbind(X, Z %*% Matrix::Matrix(diag(G))) ) %*%
			           ( Hinv %*% Matrix::t(cbind(X,Z))) )


				# Sigma^2
			  ssr <- mat$yty - t(c(b_fixed, b_random)) %*% (2*mat$u - C %*% b)
        # OJO: SE DISMINUYE EN UN GRADO DE LIBERTAD POR LA ESTIMACIÓN DE RHO
			  if(!onlygam) edf_tot <- ed1+ed2+length(b_fixed) else edf_tot <- length(b_fixed)
			  if(!rho_fixed) edf_tot <- edf_tot+1
			  if(!is.null(Xnopar)){edf_tot <- edf_tot+sum(edf_nopar)}
        sig2u <- as.numeric((ssr / (length(y) - edf_tot)))

			# New variance components and convergence check
        if(!onlygam){lanew = c(sig2u, tau1, tau2,rho)} else {lanew = c(sig2u,rho)}
			  if(!is.null(Xnopar)) {lanew = c(lanew,tau_nopar)}
				dla = mean(abs(la - lanew))
				la = lanew
			if(trace) {
				cat(sprintf("%1$3d %2$10.6f", it, dla))
				if(!onlygam) cat(sprintf("%8.3f", c(ed1, ed2)), '\n')
				cat('la ',la,"\n")
			}
			if (dla < thr) break
		} # end for (it in 1:maxit) SAP loop
		if (trace) {
		  end2 <- proc.time()[3]
		  cat("Timings:\nSAP", (end2-start2), "seconds\n")
		}

		# Se calcula el score analítico de L_{R} respecto a rho.
		if(!rho_fixed)
		{
		  score.REML_rho <- function(rho)
		  {
		    A <- as(Matrix::Diagonal(nsp)-rho*Wsp,"CsparseMatrix")
		    score.REML_rho <- t(as.matrix(P.Hv)%*%as.matrix(A%*%y))%*%as.matrix(Wsp%*%y) -
		      sum(diag(as.matrix(solve(as.matrix(A),as.matrix(Wsp)))))
		    as.numeric(score.REML_rho)
		  }

		  # Se iguala el score a 0 y se resuelve numéricamente la ecuación unidimensional
		  rho_roots <- rootSolve::multiroot(score.REML_rho,start=rho)
		  if(abs(rho_roots$f.root[1])<1.0E-6 &
		     rho_roots$root[1]<1 & rho_roots$root[1]> (-1))
		  {rho <- rho_roots$root[1]}
		} else {rho <- rho_init}

		if(!onlygam) rho_old <- la[4] else rho_old <- la[2]
		drho <- abs(rho_old-rho)
		if(!onlygam) la[4] <- rho else la[2] <- rho
		if(trace)
		  {
		    cat(sprintf("%1$3d %2$10.6f", iq, drho), '\n')
		    cat("la ",la,'\n')
		  }
		if (drho < thr) break
		} # End loop for SAP and Scores of phi and rho

  if (trace)
    {
      end <- proc.time()[3]
      cat("All process", end - start, "seconds\n")
    }

  eta <- X%*%b_fixed + Z%*%b_random + offset
  if(is.null(Xnopar)){ edf_nopar <- NULL; tau_nopar<-NULL }

  #  ESTIMATES OF PARAMETERS
  sig2u <- la[1]
  if(!onlygam)
  { tau_sp <- la[2:3]; edf_sp <- c(ed1,ed2); rho <- la[4]
  } else {
    tau_sp<-NULL; edf_sp <- NULL; rho <- la[2]
  }
  G.eff <- G # NO ESTAMOS TRABAJANDO CON PS-ANOVA
  Ginv.eff <- Ginv
  if(!is.null(Wsp)) {
      A <- as(Matrix::Diagonal(nsp)-rho*Wsp,"CsparseMatrix") } else {
      A <- Matrix::Diagonal(nsp) }

  # Valor del log.lik.reml en el óptimo
  V <- Z %*% (Matrix::Matrix(diag(G.eff)) %*% Matrix::t(Z)) +
                      sig2u*Matrix::Diagonal(nsp)
  Vinv <- try(Matrix::solve(V))
  if(class(V) == "try-error")
  {
    Vinv <- Matrix::Matrix(MASS::ginv(as.matrix(V)))
  }
  Vinv.X <- Vinv %*% X
  Xt.Vinv.X.inv <- try(Matrix::solve(Matrix::t(X) %*% (Vinv.X)))
  if(class(Xt.Vinv.X.inv)=="try-error")
  {
      Xt.Vinv.X.inv <- Matrix::Matrix(
        MASS::ginv(as.matrix(t(as.matrix(X)) %*% (Vinv.X))) )
  }
  P <- try(Vinv-Vinv.X %*% Xt.Vinv.X.inv %*% Matrix::t(Vinv.X))
  if(class(P) == "try-error")
  {
    P <- Matrix::Matrix( Vinv - (Vinv.X %*%
                  MASS::ginv( as.matrix((Matrix::t(X) %*% Vinv.X) %*%
                              Matrix::t(Vinv.X)) ) ) )
  }

  log.lik.reml <- as.numeric( -0.5*(Matrix::determinant(V)$modulus+
          Matrix::determinant(Matrix::t(X) %*% Vinv.X)$modulus +
          Matrix::t(A %*% y) %*% (P %*% (A %*%y))) +
          Matrix::determinant(A)$modulus )

  log.lik <- as.numeric( -0.5*(Matrix::determinant(V)$modulus +
                   Matrix::t( A%*% y) %*% (P %*% (A %*% y))) +
                   Matrix::determinant(A)$modulus )

  # Asignación Efectos Fijos y Aleatorios según tipo variables
  b_fixed_par <- NULL; b_fixed_nopar <- NULL; b_random_nopar <- NULL
  if(nvar_par>0) b_fixed_par <- b_fixed[1:ncol(X_par)]
  if(nvar_nopar>0)
  {
    b_fixed_nopar <- b_fixed[(ncol(X)-ncol(X_nopar)+1):ncol(X)]
    if(!onlygam) b_random_spt <- b_random[1:(ncol(Z)-ncol(Z_nopar))] else b_random_spt <- NULL
    b_random_nopar <- b_random[(ncol(Z)-ncol(Z_nopar)+1):ncol(Z)]
    if(nvar_par>0)
    {
      b_fixed_spt <- b_fixed[(ncol(X_par)+1):(ncol(X)-ncol(X_nopar))]
    } else { b_fixed_spt <- b_fixed[1:(ncol(X)-ncol(X_nopar))] }
  }	else {
    if(!onlygam) b_random_spt <- b_random else b_random_spt <- NULL
    if(nvar_par>0) { b_fixed_spt <- b_fixed[(ncol(X_par)+1):ncol(X)]
    } else { b_fixed_spt <- b_fixed[1:ncol(X)]}
  }

   if(post_estim==TRUE)
   {
  #   # CÁLCULO MATRIZ COVARIANZAS EFECTOS FIJOS Y ALEATORIOS
  #   # pp.375 Fahrmeir et al.
  #   # Se reescala A multiplicándola por sig2u toda la matriz
  #   # Matriz Covarianzas Bayesiana.
     A_cov1 <- rbind( cbind(Matrix::crossprod(X),Matrix::t(X) %*% Z),
                      cbind(Matrix::t(Z) %*% X,
                            Matrix::crossprod(Z) +
                             sig2u*Matrix::Matrix(diag(Ginv.eff))) )
     cov1_eff <- try(sig2u*Matrix::solve(A_cov1))
     if(class(cov1_eff) == "try-error") {
       cov1_eff <- sig2u*MASS::ginv(as.matrix(A_cov1))
       cov1_eff.ginv <- TRUE
     } else { cov1_eff.ginv <- FALSE }
     cov1_eff <- as.matrix(cov1_eff)

     var1_b_fixed <- diag(cov1_eff)[1:(ncol(X))]
     var1_b_random <- diag(cov1_eff)[(ncol(X)+1):(ncol(X)+ncol(Z))]
     var1_b_fixed_par <- NULL; cov1_b_fixed_par <- NULL
     var1_b_fixed_nopar <- NULL; var1_b_random_nopar <- NULL
     if(nvar_par>0)
     {
       var1_b_fixed_par <- diag(cov1_eff)[1:ncol(X_par)]
       cov1_b_fixed_par <- cov1_eff[1:ncol(X_par),1:ncol(X_par)]
     }

     if(nvar_nopar>0)
     {
       if(!onlygam)
       {
         var1_b_random_spt <- diag(cov1_eff)[(ncol(X)+1):
                                         (ncol(X)+ncol(Z)-ncol(Z_nopar))]
         var1_b_random_nopar <- diag(cov1_eff)[(ncol(X) +
                                      length(var1_b_random_spt)+1):
                                                 (ncol(X)+ncol(Z))]
       } else {
         var1_b_random_spt <- NULL
         var1_b_random_nopar <- diag(cov1_eff)[(ncol(X)+1):(ncol(X)+ncol(Z))]
       }
       if(nvar_par>0){
         var1_b_fixed_spt <- diag(cov1_eff)[(ncol(X_par)+1):
                                              (ncol(X)-ncol(X_nopar))]
         var1_b_fixed_nopar <-
                  diag(cov1_eff)[(ncol(X_par)+length(b_fixed_spt)+1):
                            ncol(X)]
       } else {
         var1_b_fixed_spt <- diag(cov1_eff)[1:ncol(X_spt)]
         var1_b_fixed_nopar <- diag(cov1_eff)[(ncol(X_spt)+1):ncol(X)]}
     } else {
       if(!onlygam)
       {
         var1_b_random_spt <- diag(cov1_eff)[(ncol(X)+1):
                                               (ncol(X)+ncol(Z))]
       } else var1_b_random_spt <- NULL
       if(nvar_par>0){
         var1_b_fixed_spt <- diag(cov1_eff)[(ncol(X_par)+1):ncol(X)]
       } else {
         var1_b_fixed_spt <- diag(cov1_eff)[1:ncol(X)]}
     }

#     # CALCULAR VALORES ESTIMADOS Y RESIDUOS
     fit_Ay <- as.vector(X %*% b_fixed + Z %*% b_random + offset)
     fit <- try(Matrix::solve(A) %*% fit_Ay)
     if(class(fit) == "try-error")
     {
       fit <- MASS::ginv(as.matrix(A)) %*% eta
     }
     fit <- as.vector(fit)

     # mH = kron(Ainv,It)*(X Z)*VAR(beta,alpha)*t(X Z)*Rinv*kron(A,It)
     mH <- (1/sig2u)*(Matrix::solve(A) %*% cbind(X,Z)) %*%
                     (cov1_eff %*% (Matrix::t(cbind(X,Z)) %*% A))
     var_fit_Ay <- as.matrix( cbind(X,Z) %*% cov1_eff %*%
                                Matrix::t(cbind(X,Z)) )
     var_fit <- try(Matrix::solve(A) %*% (var_fit_Ay %*%
                      Matrix::t(Matrix::solve(A))) )

     if(class(var_fit) == "try-error")
     {
         MASS::ginv(as.matrix(A)) %*% (var_fit_Ay %*%
                                  t(MASS::ginv(as.matrix(A))))
     }
     var_fit <- as.vector(var_fit)
     resids <- as.vector( (A %*% y) - fit_Ay )
     # Compute AIC y BIC based on loglik functions (Fahrmeir, pp. 664 and 677)
     aic <- -2*log.lik + 2*edf_tot
     bic <- -2*log.lik + log(length(y))*edf_tot
     # CÁLCULO AIC and BIC based on ssr: LIBRO FAHRMEIR pp. 564 (STAR MODELS WITH GAUSSIAN ERRORS)
     # VIP: SSR CON RESIDUOS NORMALIZADOS
     ssr <- sum(resids^2)
     aic.sig2 <- length(y)*log(ssr/length(y))+2*(edf_tot+1)
     bic.sig2 <- length(y)*log(ssr/length(y))+log(length(y))*(edf_tot+1)

#     # COMPUTE FITS EFFECTS BY VARIABLE
     fit_cov_par <- NULL;	fit_cov_nopar <- NULL
     fit_cov_nopar_fixed <- NULL;	fit_cov_nopar_random <- NULL
     var_fit_cov_nopar <- NULL; k.Z_nopar <- NULL; edf_nopar_tot <- NULL
     if(nvar_par>0)
     {
       for(i in 1:nvar_par)
       {
         fit_cov_par <- cbind(fit_cov_par,X_par[,i]*b_fixed_par[i])
       }
     }
     if(nvar_nopar>0)
     {
       edf_nopar_tot <- rep(0,nvar_nopar)
       k.Z_nopar <- rep(0,nvar_nopar)
       for(i in 1:nvar_nopar)
       {
         fit_cov_nopar_fixed <- cbind(fit_cov_nopar_fixed,
                                      X_nopar[,i]*b_fixed_nopar[i])
         var_np<-Xnopar[,i]
         # Matriz selección para varianza f(x_j) Libro Ruppert et al, pp. 175
         E_var_np <- matrix(0,ncol=ncol(cbind(X,Z)),nrow=ncol(cbind(X,Z)))
         vindex.var_np_fixed <- (ncol(X)-ncol(X_nopar)+i)
         E_var_np[vindex.var_np_fixed,vindex.var_np_fixed] <- 1

         MM.var_np <- MM_basis(var_np,min(var_np)-0.01,max(var_np)+0.01,
                               knots_nopar[i],bdeg_nopar[i],pord_nopar[i])
         X.var_np <- MM.var_np$X; Z.var_np <- MM.var_np$Z
         k.Z_nopar[i] <- ncol(Z.var_np)
         rm(MM.var_np)
         if(i==1){
           fit_cov_nopar_random <- cbind(fit_cov_nopar_random,
                                     Z.var_np%*%b_random_nopar[1:k.Z_nopar[1]])
           vindex.var_np_random <- c((ncol(X)+ncolZ_spt+1):
                                     (ncol(X)+ncolZ_spt+k.Z_nopar[1]))
         } else {
           fit_cov_nopar_random <- cbind(fit_cov_nopar_random,
                       Z.var_np%*%b_random_nopar[(sum(k.Z_nopar[1:(i-1)])+1):
                                                         sum(k.Z_nopar[1:i])])
           vindex.var_np_random<- c((ncol(X)+ncolZ_spt+
                                       sum(k.Z_nopar[1:(i-1)])+1):
                                      (ncol(X)+ncolZ_spt+sum(k.Z_nopar[1:i])))
         }
         # Se cambia la diagonal de la matriz de Selección
         for(j in 1:length(vindex.var_np_random))
         {
           E_var_np[vindex.var_np_random[j],vindex.var_np_random[j]] <- 1
         }
         mH_var_np <- as.matrix( (1/sig2u)*(Matrix::solve(A) %*% cbind(X,Z) %*%
                                   E_var_np) %*% (cov1_eff %*%
                                    (Matrix::t(cbind(X,Z))%*%A)) )

         edf_nopar_tot[i] <- sum(diag(mH_var_np))
         var_fit_np <- rep(0,length(fit))
         for(j in 1:length(fit))
         {
           var_fit_np[j] <- sig2u*sum(mH_var_np[j,]^2)
         }
         var_fit_cov_nopar <- cbind(var_fit_cov_nopar,var_fit_np)
       }
       fit_cov_nopar <- fit_cov_nopar_fixed+fit_cov_nopar_random
     }

#     # COMPUTE SPATIO-TEMPORAL TREND.
     if(!onlygam)
     {
       spt_trend_fixed <- as.vector(X_spt%*%b_fixed_spt)
       spt_trend_random <- as.vector(Z_spt%*%b_random_spt)
       spt_trend <- spt_trend_fixed+spt_trend_random
       E_spt_trend <- matrix(0,ncol=ncol(cbind(X,Z)),nrow=ncol(cbind(X,Z)))
       vindex_spt_trend_fixed <- c((nvar_par+1):(ncol(X)-nvar_nopar))
       for(i in 1:length(vindex_spt_trend_fixed))
       {
         E_spt_trend[vindex_spt_trend_fixed[i],
                     vindex_spt_trend_fixed[i]] <- 1
       }
       vindex_spt_trend_random <- c((ncol(X)+1):(ncol(X)+ncolZ_spt))
       for(i in 1:length(vindex_spt_trend_random))
       {
         E_spt_trend[vindex_spt_trend_random[i],
                     vindex_spt_trend_random[i]] <- 1
       }
       mH_spt_trend <- as.matrix( (1/sig2u)*(Matrix::solve(A) %*% cbind(X,Z) %*%
                                    E_spt_trend) %*%
                      (cov1_eff %*% (Matrix::t(cbind(X,Z)) %*% A)) )
       edf_spt_tot <- sum(diag(mH_spt_trend))
       var_spt_trend <- rep(0,length(spt_trend))
       for(j in 1:length(spt_trend))
       {
         var_spt_trend[j] <- sig2u*sum(mH_spt_trend[j,]^2)
       }
       spt_trend_fixed <- matrix(spt_trend_fixed,ncol=nsp,byrow=FALSE)
       spt_trend_random <- matrix(spt_trend_random,ncol=nsp,byrow=FALSE)
       spt_trend <- matrix(spt_trend,ncol=nsp,byrow=FALSE)
       var_spt_trend <- matrix(var_spt_trend,ncol=nsp,byrow=FALSE)
     } else {
       spt_trend_fixed <- as.vector(X_spt%*%b_fixed_spt)
       spt_trend <- spt_trend_fixed
       spt_trend_random <- NULL
       var_spt_trend <- NULL
       edf_spt_tot <- 1
     }
   } # end if(post_estim==TRUE)


	res <- list (call = match.call(),
	             # Parámetros estimados
               knots = knots, bdeg = bdeg, pord = pord,
	             knots_nopar=knots_nopar,bdeg_nopar=bdeg_nopar,pord_nopar=pord_nopar,
               tau_sp = tau_sp,tau_nopar=tau_nopar,
	             edf_sp = edf_sp,edf_nopar=edf_nopar,edf_tot=edf_tot,
	             edf_spt_tot=edf_spt_tot,
	             sig2u=sig2u,rho=rho,
	             ssr=ssr,aic=aic,bic=bic,
	             aic.sig2=aic.sig2,bic.sig2=bic.sig2,
	             log.lik.reml=log.lik.reml,log.lik=log.lik,
	             b_fixed=b_fixed,b_random=b_random,
	             b_fixed_spt=b_fixed_spt,b_random_spt=b_random_spt,
	             b_fixed_par=b_fixed_par,b_fixed_nopar=b_fixed_nopar,
	             b_random_nopar=b_random_nopar,
	             var1_b_fixed=var1_b_fixed,
	             var1_b_random=var1_b_random,
	             var1_b_fixed_spt=var1_b_fixed_spt,
	             var1_b_random_spt=var1_b_random_spt,
	             var1_b_fixed_par=var1_b_fixed_par,
	             cov1_b_fixed_par=cov1_b_fixed_par,
	             var1_b_fixed_nopar=var1_b_fixed_nopar,
	             var1_b_random_nopar=var1_b_random_nopar,
	             # fit, resids y var(fit)
	             fit=fit,fit_Ay=fit_Ay,resids=resids,
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
