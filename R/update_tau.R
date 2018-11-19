###############################################################################
# Function to update tau's in each iteration of the loop.
update_tau2d <- function (la, lg, G, dZtPZ_wide, brandom_wide,
                          nvarnopar = 0, dnoparlist = NULL, psanova = TRUE,
                          f1_main = TRUE, f2_main = TRUE, f12_int = TRUE,
                          tau_init = NULL, tau_fixed = NULL,
                          taunopar_fixed = NULL, taunopar_init = NULL){
  np <- lg$np
  G1inv.n <- lg$G1inv.n; G2inv.n <- lg$G2inv.n
  g1u <- lg$g1u; g2u <- lg$g2u; g1b <- lg$g1b; g2b <- lg$g2b
  if (!is.null(g1u) & !is.null(g2u) & !is.null(g1b) & !is.null(g2b)) {
    # With Spatial Trend
    if (!psanova) {
      if (is.null(tau_init)) tau_init <- rep(1,2)
      if (is.null(tau_fixed)) tau_fixed <- rep(FALSE,2)
      # Tau 1
      G1inv.d <- (1/la[2])*G1inv.n
      ed1 <- sum(dZtPZ_wide*(G1inv.d*G^2))
      ed1 <- ifelse(ed1 == 0, 1e-50,ed1)
      tau1 <- ifelse(!tau_fixed[1],
                     sum(brandom_wide^2*G1inv.n)/ed1,tau_init[1])
      # Tau 2
      G2inv.d <- (1/la[3])*G2inv.n
      ed2 <- sum(dZtPZ_wide*(G2inv.d*G^2))
      ed2 <- ifelse(ed2 == 0, 1e-50, ed2)
      tau2 <- ifelse(!tau_fixed[2],
                     sum(brandom_wide^2*G2inv.n)/ed2,tau_init[2])
    } else { # PS-ANOVA=TRUE
      if (is.null(tau_init)) tau_init <- rep(1,4)
      if (is.null(tau_fixed)) tau_fixed <- rep(FALSE,4)
      G3inv.n <- lg$G3inv.n; G4inv.n <- lg$G4inv.n
      # Tau 1, Tau 2, Tau 3 (Main Effects)
      if(f1_main){
        G1inv.d <- (1/la[2])*G1inv.n
        ed1 <- sum(dZtPZ_wide*(G1inv.d*G^2))
        ed1 <- ifelse(ed1 == 0, 1e-50,ed1)
        tau1 <- ifelse(!tau_fixed[1],
                       sum(brandom_wide^2*G1inv.n)/ed1,tau_init[1])
      } else { tau1 <- 0; ed1 <- 0 }
      if(f2_main){
        G2inv.d <- (1/la[3])*G2inv.n
        ed2 <- sum(dZtPZ_wide*(G2inv.d*G^2))
        ed2 <- ifelse(ed2 == 0, 1e-50, ed2)
        tau2 <- ifelse(!tau_fixed[2],
                       sum(brandom_wide^2*G2inv.n)/ed2,tau_init[2])
      } else { tau2 <- 0; ed2 <- 0 }
      # Interactions 2nd order
      if (f12_int) {
        G3inv.d <- (1/la[4])*G3inv.n
        ed3 <- sum(dZtPZ_wide*(G3inv.d*G^2))
        ed3 <- ifelse(ed3 == 0, 1e-50,ed3)
        tau3 <- ifelse(!tau_fixed[3],
                       sum(brandom_wide^2*G3inv.n)/ed3,tau_init[3])
        G4inv.d <- (1/la[5])*G4inv.n
        ed4 <- sum(dZtPZ_wide*(G4inv.d*G^2))
        ed4 <- ifelse(ed4 == 0, 1e-50,ed4)
        tau4 <- ifelse(!tau_fixed[4],
                       sum(brandom_wide^2*G4inv.n)/ed4,tau_init[4])
      } else { tau3 <- 0; tau4 <- 0; ed3 <- 0; ed4 <- 0 }
    }
  } else { # Without Spatial Trend
    tau1 <- tau2 <- tau3 <- tau4 <- NULL
    ed1 <- ed2 <- ed3 <- ed4 <- NULL
  }
  if (nvarnopar>0) {
    if (is.null(taunopar_fixed)) taunopar_fixed <- rep(FALSE,nvarnopar)
    Ginv_dnopar <- matrix(0, nrow=sum(np[2:length(np)]),ncol=nvarnopar)
    if (!is.null(taunopar_init)) {
      taunopar <- taunopar_init  } else taunopar <- rep(1,nvarnopar)
    edfnopar <- rep(0,nvarnopar)
    for (k in 1:nvarnopar) {
      if (!is.null(g1u) & !is.null(g2u) & !is.null(g1b) & !is.null(g2b)) {
        # With Spatial Trend
        if(!psanova){
          Ginv_dnopar[(sum(np[2:(4+k-1)])+1):sum(np[2:(4+k)]),k] <-
            (1/la[length(la)-(nvarnopar+1)+k])*dnoparlist[[k]]
        } else {
          Ginv_dnopar[(sum(np[2:(6+k-1)])+1):sum(np[2:(6+k)]),k] <-
            (1/la[length(la)-(nvarnopar+1)+k])*dnoparlist[[k]]
        }
      } else { # Without Spatial Trend
        if (k==1){
          Ginv_dnopar[1:np[2],k] <-
             (1/la[length(la)-(nvarnopar+1)+k])*dnoparlist[[k]]
        } else {
          Ginv_dnopar[(sum(np[2:k])+1):(sum(np[2:(k+1)])),k] <-
            (1/la[length(la)-(nvarnopar+1)+k])*dnoparlist[[k]]
        }
      }
      edfnopar[k] <- sum(dZtPZ_wide*(Ginv_dnopar[,k]*G^2))
      edfnopar[k] <- ifelse(edfnopar[k] == 0, 1e-50, edfnopar[k])
      taunopar[k] <- ifelse(!taunopar_fixed[k],
                            sum(brandom_wide^2*Ginv_dnopar[,k]) /
                              edfnopar[k], taunopar_init[k])
      taunopar[k] <- ifelse(taunopar[k] == 0, 1e-50, taunopar[k])
    } # end for (k in 1:nvarnopar)
  } else { taunopar <- NULL; edfnopar <- NULL }
  if (psanova) {
    res <- list( tau1 = tau1, tau2 = tau2, tau3 = tau3, tau4 = tau4,
                 ed1 = ed1, ed2 = ed2, ed3 = ed3, ed4 = ed4,
                 taunopar = taunopar, edfnopar = edfnopar )
  } else {
    res <- list( tau1 = tau1, tau2 = tau2, ed1 = ed1, ed2 = ed2,
                 taunopar = taunopar, edfnopar = edfnopar )
  }
  res
}



###############################################################################
# Function to update tau's in each iteration of the loop.
update_tau3d <- function (la, lg, G, dZtPZ_wide, brandom_wide,
                     nvarnopar = 0, dnoparlist = NULL, psanova = TRUE,
                     f1_main = TRUE, f2_main = TRUE, ft_main = TRUE,
                     f12_int = TRUE,f1t_int = TRUE, f2t_int = TRUE,
                     f12t_int=TRUE,
                     tau_init = NULL, tau_fixed = NULL,
                     taunopar_fixed = NULL, taunopar_init = NULL){
    np <- lg$np
    G1inv.n <- lg$G1inv.n; G2inv.n <- lg$G2inv.n; G3inv.n <- lg$G3inv.n;
    if (!psanova) {
      if (is.null(tau_init)) tau_init <- rep(1,3)
      if (is.null(tau_fixed)) tau_fixed <- rep(FALSE,3)
        # Tau 1
        G1inv.d <- (1/la[2])*G1inv.n
        ed1 <- sum(dZtPZ_wide*(G1inv.d*G^2))
        ed1 <- ifelse(ed1 == 0, 1e-50,ed1)
        tau1 <- ifelse(!tau_fixed[1],
                       sum(brandom_wide^2*G1inv.n)/ed1,tau_init[1])
        # Tau 2
        G2inv.d <- (1/la[3])*G2inv.n
        ed2 <- sum(dZtPZ_wide*(G2inv.d*G^2))
        ed2 <- ifelse(ed2 == 0, 1e-50, ed2)
        tau2 <- ifelse(!tau_fixed[2],
                       sum(brandom_wide^2*G2inv.n)/ed2,tau_init[2])
        # Tau 3
        G3inv.d <- (1/la[4])*G3inv.n
        ed3 <- sum(dZtPZ_wide*(G3inv.d*G^2))
        ed3 <- ifelse(ed3 == 0, 1e-50,ed3)
        tau3 <- ifelse(!tau_fixed[3],
                       sum(brandom_wide^2*G3inv.n)/ed3,tau_init[3])
    } else { # PS-ANOVA=TRUE
      if (is.null(tau_init)) tau_init <- rep(1,12)
      if (is.null(tau_fixed)) tau_fixed <- rep(FALSE,12)
        G4inv.n <- lg$G4inv.n; G5inv.n <- lg$G5inv.n; G6inv.n <- lg$G6inv.n
        G7inv.n <- lg$G7inv.n; G8inv.n <- lg$G8inv.n; G9inv.n <- lg$G9inv.n
        G10inv.n <- lg$G10inv.n; G11inv.n <- lg$G11inv.n
        G12inv.n <- lg$G12inv.n
        # Tau 1, Tau 2, Tau 3 (Main Effects)
        if(f1_main){
            G1inv.d <- (1/la[2])*G1inv.n
            ed1 <- sum(dZtPZ_wide*(G1inv.d*G^2))
            ed1 <- ifelse(ed1 == 0, 1e-50,ed1)
            tau1 <- ifelse(!tau_fixed[1],
                           sum(brandom_wide^2*G1inv.n)/ed1,tau_init[1])
        } else { tau1 <- 0; ed1 <- 0 }
        if(f2_main){
            G2inv.d <- (1/la[3])*G2inv.n
            ed2 <- sum(dZtPZ_wide*(G2inv.d*G^2))
            ed2 <- ifelse(ed2 == 0, 1e-50, ed2)
            tau2 <- ifelse(!tau_fixed[2],
                           sum(brandom_wide^2*G2inv.n)/ed2,tau_init[2])
        } else { tau2 <- 0; ed2 <- 0 }
        if (ft_main) {
            G3inv.d <- (1/la[4])*G3inv.n
            ed3 <- sum(dZtPZ_wide*(G3inv.d*G^2))
            ed3 <- ifelse(ed3 == 0, 1e-50,ed3)
            tau3 <- ifelse(!tau_fixed[3],
                           sum(brandom_wide^2*G3inv.n)/ed3,tau_init[3])
        } else {  tau3 <- 0; ed3 <- 0 }
        # Interactions 2nd order
        if (f12_int) {
            G4inv.d <- (1/la[5])*G4inv.n
            ed4 <- sum(dZtPZ_wide*(G4inv.d*G^2))
            ed4 <- ifelse(ed4 == 0, 1e-50,ed4)
            tau4 <- ifelse(!tau_fixed[4],
                           sum(brandom_wide^2*G4inv.n)/ed4,tau_init[4])
            G5inv.d <- (1/la[6])*G5inv.n
            ed5 <- sum(dZtPZ_wide*(G5inv.d*G^2))
            ed5 <- ifelse(ed5 == 0, 1e-50,ed5)
            tau5 <- ifelse(!tau_fixed[5],
                           sum(brandom_wide^2*G5inv.n)/ed5,tau_init[5])
        } else { tau4 <- 0; tau5 <- 0; ed4 <- 0; ed5 <- 0 }
        if (f1t_int) {
            G6inv.d <- (1/la[7])*G6inv.n
            ed6 <- sum(dZtPZ_wide*(G6inv.d*G^2))
            ed6 <- ifelse(ed6 == 0, 1e-50,ed6)
            tau6 <- ifelse(!tau_fixed[6],
                           sum(brandom_wide^2*G6inv.n)/ed6,tau_init[6])
            G7inv.d <- (1/la[8])*G7inv.n
            ed7 <- sum(dZtPZ_wide*(G7inv.d*G^2))
            ed7 <- ifelse(ed7 == 0, 1e-50,ed7)
            tau7 <- ifelse(!tau_fixed[7],
                           sum(brandom_wide^2*G7inv.n)/ed7,tau_init[7])
        } else { tau6 <- 0; tau7 <- 0; ed6 <- 0; ed7 <- 0 }
        if(f2t_int){
            G8inv.d <- (1/la[9])*G8inv.n
            ed8 <- sum(dZtPZ_wide*(G8inv.d*G^2))
            ed8 <- ifelse(ed8 == 0, 1e-50,ed8)
            tau8 <- ifelse(!tau_fixed[8],
                           sum(brandom_wide^2*G8inv.n)/ed8,tau_init[8])
            G9inv.d <- (1/la[10])*G9inv.n
            ed9 <- sum(dZtPZ_wide*(G9inv.d*G^2))
            ed9 <- ifelse(ed9 == 0, 1e-50,ed9)
            tau9 <- ifelse(!tau_fixed[9],
                           sum(brandom_wide^2*G9inv.n)/ed9,tau_init[9])
        } else { tau8 <- 0; tau9 <- 0; ed8 <- 0; ed9 <- 0 }
        # Interactions 3rd order
        if(f12t_int){
            G10inv.d <- (1/la[11])*G10inv.n
            ed10 <- sum(dZtPZ_wide*(G10inv.d*G^2))
            ed10 <- ifelse(ed10 == 0, 1e-50,ed10)
            tau10 <- ifelse(!tau_fixed[10],
                            sum(brandom_wide^2*G10inv.n)/ed10,tau_init[10])
            G11inv.d <- (1/la[12])*G11inv.n
            ed11 <- sum(dZtPZ_wide*(G11inv.d*G^2))
            ed11 <- ifelse(ed11 == 0, 1e-50,ed11)
            tau11 <- ifelse(!tau_fixed[11],
                            sum(brandom_wide^2*G11inv.n)/ed11,tau_init[11])
            G12inv.d <- (1/la[13])*G12inv.n
            ed12 <- sum(dZtPZ_wide*(G12inv.d*G^2))
            ed12 <- ifelse(ed12 == 0, 1e-50,ed12)
            tau12 <- ifelse(!tau_fixed[12],
                            sum(brandom_wide^2*G12inv.n)/ed12,tau_init[12])
        } else { tau10 <- 0; tau11 <- 0; tau12 <- 0
        ed10 <- 0; ed11 <- 0; ed12 <- 0 }
    }

    if (nvarnopar>0) {
        if (is.null(taunopar_fixed)) taunopar_fixed <- rep(FALSE,nvarnopar)
        Ginv_dnopar <- matrix(0, nrow=sum(np[2:length(np)]),
                               ncol=nvarnopar)
        if (!is.null(taunopar_init)) {
          taunopar <- taunopar_init  } else taunopar <- rep(1,nvarnopar)
        #if (!taunopar_fixed) { taunopar <- rep(1,nvarnopar) } else {
        #    taunopar <- taunopar_init  }
        edfnopar <- rep(0,nvarnopar)
        for (k in 1:nvarnopar) {
            if(!psanova){
                Ginv_dnopar[(sum(np[2:(8+k-1)])+1):sum(np[2:(8+k)]),k] <-
                    (1/la[length(la)-(nvarnopar+2)+k])*dnoparlist[[k]]
            } else {
                Ginv_dnopar[(sum(np[2:(20+k-1)])+1):sum(np[2:(20+k)]),k] <-
                    (1/la[length(la)-(nvarnopar+2)+k])*dnoparlist[[k]]
            }
            edfnopar[k] <- sum(dZtPZ_wide*(Ginv_dnopar[,k]*G^2))
            edfnopar[k] <- ifelse(edfnopar[k] == 0, 1e-50, edfnopar[k])
            taunopar[k] <- ifelse(!taunopar_fixed[k],
                                   sum(brandom_wide^2*Ginv_dnopar[,k]) /
                                       edfnopar[k], taunopar_init[k])
            taunopar[k] <- ifelse(taunopar[k] == 0, 1e-50, taunopar[k]) }
        # end for (k in 1:nvarnopar)
    } else {
        taunopar <- NULL; edfnopar <- NULL
    }
    if (psanova) {
      res <- list( tau1 = tau1, tau2 = tau2, tau3 = tau3, tau4 = tau4,
                   tau5 = tau5, tau6 = tau6, tau7 = tau7, tau8 = tau8,
                   tau9 = tau9, tau10 = tau10, tau11 = tau11, tau12 = tau12,
                   ed1 = ed1, ed2 = ed2, ed3 = ed3, ed4 = ed4,
                   ed5 = ed5, ed6 = ed6, ed7 = ed7, ed8 = ed8,
                   ed9 = ed9, ed10 = ed10, ed11 = ed11, ed12 = ed12,
                   taunopar = taunopar, edfnopar = edfnopar )
    } else {
      res <- list( tau1 = tau1, tau2 = tau2, tau3 = tau3,
                   ed1 = ed1, ed2 = ed2, ed3 = ed3,
                   taunopar = taunopar, edfnopar = edfnopar )
    }
    res
}
