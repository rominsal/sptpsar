#' @name eff_par
#' @rdname eff_par
#'
#' @title Compute direct, indirect and total effects (or impacts)
#'   for continous parametric covariates in spatial or 
#'   spatio-temporal semiparametric PS-SAR regression models.
#'        
#' @description Compute direct, indirect and total effects (also named
#'   impacts) for non-parametric covariates included in a semiparametric spatial
#'   or spatio-temporal SAR model. This model must include a spatial
#'   lag of the dependent variable (SAR) to have indirect effects 
#'   different from 0, otherwise, total and direct effects 
#'   are the same.         
#'
#' @param sptsarfit A \emph{psar} object fitted using \code{\link{psar}} function 
#' @param variables vector including names of non-parametric covariates.
#' @param nrep number of repetitions for the simulation. Default 1000.
#' @param seed initial seed to get random numbers. Must be set to a 
#'             specific value to make reproducible results. Default 1111.
#' @param m number of powers to compute a vector of traces of powers
#'          of a spatial weight matrix (see \code{\link[spdep]{trW}} in 
#'          \pkg{spdep} package). Default 100.
#' @param p number of samples used in MC simulation of traces
#'          of a spatial weight matrix (see \code{\link[spdep]{trW}} in
#'          \pkg{spdep} package).Default 50.
#' @param tol tolerance (relative to largest variance) for numerical lack 
#'            of positive-definiteness in Sigma when simulate \eqn{\betas}
#'            from the maximum likelihood estimates (see \code{\link[MASS]{mvrnorm}}
#'            in \pkg{MASS} package). Default 0.01.
#'                     
#' @details DESCRIBE ALGORITHM TO SIMULATE PARAMETRIC EFFECTS
#' 
#' @return An object of class \emph{par.eff.psar}. Can be printed
#'         with \code{summary}.
#'         
#'         The object returned is a list with 3 matrices including
#'         the results of simulated effects:
#'          \tabular{ll}{
#'            \emph{tot_eff} \tab Matrix including simulated total effects 
#'                                for each variable in rows. \cr
#'            \emph{dir_eff} \tab Matrix including simulated direct effects 
#'                                for each variable in rows. \cr
#'            \emph{ind_eff} \tab Matrix including simulated indirect effects 
#'                                for each variable in rows. \cr
#'          }                      
#' 
#' @seealso
#' \itemize{
#'   \item \code{\link{psar}} estimate spatial or spatio-temporal 
#'           semiparametric PS-SAR regression models.
#'   \item \code{\link{eff_nopar}} compute total, direct and indirect effect
#'           functions for non-parametric continuous covariates.
#'   \item \code{\link{fit_terms}} compute smooth functions for non-parametric
#'           continuous covariates.
#'   \item \code{\link[spdep]{impacts}} similar function in \pkg{spdep} 
#'           package to compute impacts in spatial parametric econometric 
#'           models.                              
#' }
#' 
#' @references \itemize{ 
#'   \item LeSage, J. and Pace, K. (2009). \emph{Introduction to 
#'         Spatial Econometrics}. CRC Press, Boca Raton.
#'  }       
#'
#' @examples
#' ################################################
#'  ###################### Examples using a panel data of rate of
#'  ###################### unemployment for 103 Italian provinces in period 1996-2014.
#' library(sptpsar)
#' data(unemp_it); Wsp <- Wsp_it
#' 
#' ######################  No Spatial Trend: PSAR including a spatial 
#' ######################  lag of the dependent variable
#' form1 <- unrate ~ partrate + agri + cons +
#'                  pspl(serv,nknots=15) +
#'                  pspl(empgrowth,nknots=20) 
#'  gamsar <- psar(form1,data=unemp_it,sar=TRUE,Wsp=Wsp_it)
#'  summary(gamsar)
#'  ###### Parametric Total, Direct and Indirect Effects
#'  list_varpar <- c("partrate","agri","cons")
#'  eff_parvar <- eff_par(gamsar,list_varpar)
#'  summary(eff_parvar)
#'  
#' ######################   PSAR-ANOVA with spatial trend
#' form2 <- unrate ~ partrate + agri + cons +
#'                   pspl(serv,nknots=15) + pspl(empgrowth,nknots=20) +
#'                   pspt(long,lat,nknots=c(20,20),psanova=TRUE,
#'                   nest_sp1=c(1,2),nest_sp2=c(1,2))
#' ##### Spatial trend fixed for period 1996-2014
#' geospanova_sar <- psar(form2,data=unemp_it,Wsp=Wsp_it,sar=TRUE,
#'                    control=list(thr=1e-1,maxit=200,trace=FALSE))
#' summary(geospanova_sar)
#'  ###### Parametric Total, Direct and Indirect Effects
#'  list_varpar <- c("partrate","agri","cons")
#'  eff_parvar <- eff_par(geospanova_sar,list_varpar)
#'  summary(eff_parvar)
#'  
#' ######################   PSAR-ANOVA with spatio-temporal trend and 
#' ######################   temporal autorregresive noise
#'  form3 <- unrate ~ partrate + agri + cons +
#'                    pspl(serv,nknots=15) + pspl(empgrowth,nknots=20) +
#'                    pspt(long,lat,year,nknots=c(18,18,8),psanova=TRUE,
#'                    nest_sp1=c(1,2,3),nest_sp2=c(1,2,3),
#'                    nest_time=c(1,2,2),ntime=19)
#' sptanova_sar_ar1 <- psar(form3,data=unemp_it,Wsp=Wsp_it,sar=TRUE,ar1=TRUE,
#'                     control=list(thr=1e-1,maxit=200,trace=FALSE))
#' summary(sptanova_sar_ar1)
#'  ###### Parametric Total, Direct and Indirect Effects
#'  list_varpar <- c("partrate","agri","cons")
#'  eff_parvar <- eff_par(sptanova_sar_ar1,list_varpar)
#'  summary(eff_parvar)
#'
#' @keywords Indirect effects, Direct effects, SAR, parametric covariates.
#'
#' @export
#'
eff_par <- function(sptsarfit,variables,
                    nrep=1000,seed=1111,
                    m=100,p=50,tol=0.01){
  bfixed <- sptsarfit$bfixed
  match_varpar <- unique(grep(paste(variables,
                                    collapse="|"),
                              names(bfixed), value=TRUE))
  bfixed_par <- bfixed[match_varpar]

  cov_b <- sptsarfit$vcov_b
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
    trWsp <- spdep::trW(as(Wsp,"CsparseMatrix"),m=m,p=p,type="MC")
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

