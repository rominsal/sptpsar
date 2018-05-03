#' Compute Total, Direct and Indirect effects functions for  geoadditive
#' spatial or spatio-temporal semiparametric PS-SAR regression models.
#'
#' Compute effects function for non-parametric covariates in semiparametric
#' models.
#'
#' @param fitted_nopar A matrix including the non-parametric fitted
#'    functions (GAM) for each covariate.
#' @param sd_fitted_nopar A matrix including the standard deviations of
#'    non-parametric functions(GAM) for each convariate.
#' @param Wsp A neighbour spatial matrix.
#'    The dimension of the matrix is always \emph{nxn}, where \emph{n} is the
#'    dimension of each spatial coordinate. Default NULL.
#' @param rho Spatial parameter of the spatial lag of the dependent
#'    variable for SAR models. Default 0.
#' @param durbin A logical value indicating if the model include a spatial lag
#'    of each covariate (either parametric or non-parametric). Default FALSE.
#' @param conflevel Numerical value for the confidence interval of the
#'    effect functions.
#'
eff_nopar <- function(sptsarfit,variables,conflevel=0.95)
{
  sar <- sptsarfit$sar
  if(sar){
     rho <- sptsarfit$rho
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
  } # end if (sar)

  crval <- qnorm((1-conflevel)/2,mean=0,sd=1,lower.tail=FALSE)
  tot_eff <- dir_eff <- ind_eff <- NULL
  uptot_eff <- updir_eff <- upind_eff <- NULL
  lowtot_eff <- lowdir_eff <- lowind_eff <- NULL
  for (i in 1:length(variables)) {
    var_name <- variables[i]
    fitted_terms_i <- fit_terms(sptsarfit,var_name)
    fit_i <- fitted_terms_i$fitted_terms
    colnames(fit_i) <- var_name

    tot_eff_i <- as.matrix(Matrix::solve(A,fit_i))
    colnames(tot_eff_i) <- var_name

    dir_eff_i <- as.matrix( diag(
            as.matrix(Matrix::solve(A,diag(as.vector(fit_i))))) )
    colnames(dir_eff_i) <- var_name

    ind_eff_i <- tot_eff_i - dir_eff_i
    colnames(ind_eff_i) <- var_name

    tot_eff <- cbind(tot_eff,tot_eff_i)
    dir_eff <- cbind(dir_eff,dir_eff_i)
    ind_eff <- cbind(ind_eff,ind_eff_i)

    se_fit_i <- fitted_terms_i$se_fitted_terms
    colnames(se_fit_i) <- var_name

    upfit_i <- fit_i + crval*se_fit_i
    uptot_eff_i <- as.matrix(Matrix::solve(A,upfit_i))
    colnames(uptot_eff_i) <- var_name
    updir_eff_i <- as.matrix( diag(
      as.matrix(Matrix::solve(A,diag(as.vector(upfit_i))))) )
    colnames(updir_eff_i) <- var_name
    upind_eff_i <- uptot_eff_i - updir_eff_i
    colnames(upind_eff_i) <- var_name
    uptot_eff <- cbind(uptot_eff,uptot_eff_i)
    updir_eff <- cbind(updir_eff,updir_eff_i)
    upind_eff <- cbind(upind_eff,upind_eff_i)

    lowfit_i <- fit_i - crval*se_fit_i
    lowtot_eff_i <- as.matrix(Matrix::solve(A,lowfit_i))
    colnames(lowtot_eff_i) <- var_name
    lowdir_eff_i <- as.matrix( diag(
      as.matrix(Matrix::solve(A,diag(as.vector(lowfit_i))))) )
    colnames(lowdir_eff_i) <- var_name
    lowind_eff_i <- lowtot_eff_i - lowdir_eff_i
    colnames(lowind_eff_i) <- var_name
    lowtot_eff <- cbind(lowtot_eff,lowtot_eff_i)
    lowdir_eff <- cbind(lowdir_eff,lowdir_eff_i)
    lowind_eff <- cbind(lowind_eff,lowind_eff_i)
  }
   res <- list(effnopar_tot = tot_eff,
               effnopar_dir = dir_eff,
               effnopar_ind = ind_eff,
               effnopar_tot_up = uptot_eff,
               effnopar_dir_up = updir_eff,
               effnopar_ind_up = upind_eff,
               effnopar_tot_low = lowtot_eff,
               effnopar_dir_low = lowdir_eff,
               effnopar_ind_low = lowind_eff)
   res
}



