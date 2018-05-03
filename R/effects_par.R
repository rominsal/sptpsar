#' Compute Direct, Indirect and Total effects for parametric covariates
#' in spatio-temporal semiparametric regression models.

#' @param sptsarfit An object of class \emph{psar} usually fitted using command
#'                  \code{\link{psar}}


eff_par <- function(sptsarfit,variables,
                    nrep=1000,seed=1111,
                    m=100,p=50,tol=0.01){
# Function to estimate Direct, Indirect and Total Effects
#  for Parametric Covariates

  bfixed <- sptsarfit$bfixed
  match_varpar <- unique(grep(paste(variables,
                                    collapse="|"),
                              names(bfixed), value=TRUE))
  bfixed_par <- bfixed[match_varpar]

  cov_b <- sptsarfit$var_betas_alphas
  row_cov_fixed <- c(grepl("fixed",rownames(cov_b)))
  col_cov_fixed <- c(grepl("fixed",colnames(cov_b)))
  cov_bfixed <- cov_b[row_cov_fixed,col_cov_fixed]
  match_varpar <- unique(grep(paste(variables,
                                    collapse="|"),
                              colnames(cov_bfixed), value=TRUE))
  cov_bfixed_par <- cov_bfixed[match_varpar,match_varpar]

  sar <- sptsarfit$sar
  if(sar){
    rho <- sptsarfit$rho
    se_rho <- sptsarfit$se_rho
    Wsp <- Matrix::Matrix(sptsarfit$Wsp)
    nfull <- length(sptsarfit$y)
    time <- sptsarfit$time
    if(is.null(time)){ # Case 2d
      nsp_short <- nrow(Wsp)
      Wsp_full <- Wsp %x% Matrix::Diagonal(nfull/nsp_short)
      In <- Matrix::Diagonal(nfull)
      A <- In - rho*Wsp_full
    } else { # Case 3d
      nsp <- length(sptsarfit$sp1)
      ntime <- length(time)
      In <- Matrix::Diagonal(nsp)
      It <- Matrix::Diagonal(ntime)
      A <- In - rho*Wsp
      A <- A %x% It
    }
  } else {
    A <- Matrix::Diagonal(nfull)
    stop ("Model fitted must be sar")
  } # end if (sar)

  tot_eff <- dir_eff <- ind_eff <- NULL
  if(!is.null(seed)) set.seed(seed)
  if(!is.null(Wsp)){
    trWsp <- spdep::trW(Wsp,m=m,p=p,type="MC")
  } else stop("W matrix is null")

  rho_sim <- rnorm(nrep,rho,se_rho)
  beta_par_sim <- MASS::mvrnorm(nrep,bfixed_par,cov_bfixed_par,tol=tol)
  beta_par_sim <- t(beta_par_sim)
  rownames(beta_par_sim) <- variables
  q <- length(trWsp)
  nsp <- nrow(Wsp)
  mT <- matrix(c(1,trWsp[1:q]/nsp),nrow=1)
  a <- matrix(1,nrow=q+1,ncol=1)
  #   G <- diag(g)
  mP <- list()
  k <- length(variables)
  dir <- matrix(NA,nrow=k,ncol=nrep)
  tot <- matrix(NA,nrow=k,ncol=nrep)
  for (i in 1:nrep)
  {
    mP[[i]] <- matrix(beta_par_sim[1:k,i],ncol=1)
    g <- 1
    for(j in 1:q) { g <- c(g,rho_sim[i]^j) }
    G <- diag(g)
    dir[,i] <- mP[[i]] %*% mT %*% G %*% a
    tot[,i] <- matrix(rowSums(mP[[i]]),ncol=1) %*%
                      (matrix(g,nrow=1) %*% a)
  }
  ind <- tot - dir
  rownames(dir) <- rownames(ind) <- rownames(tot) <- variables
  res <- list( tot_eff = tot, dir_eff = dir, ind_eff = ind)
  class(res) <- ("par.eff.psar")
  res
}

