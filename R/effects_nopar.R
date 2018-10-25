#' @name eff_nopar
#' @rdname eff_nopar
#'
#' @title Compute direct, indirect and total effect functions 
#'   for continous non-parametric covariates in spatial or 
#'   spatio-temporal semiparametric PS-SAR regression models.
#'        
#' @description Compute direct, indirect and total effect functions for 
#'   non-parametric covariates included in a semiparametric spatial
#'   or spatio-temporal SAR model. This model must include a spatial
#'   lag of the dependent variable (SAR) to have indirect effects 
#'   different from 0, otherwise, total and direct function effects 
#'   are the same.         
#'
#' @param sptsarfit \emph{psar} object fitted using \code{\link{psar}} function 
#' @param variables vector including names of non-parametric covariates.
#' @param conflevel mumerical value for the confidence interval of the
#'    effect functions. Default 0.95.
#'
#' @details DESCRIBE ALGORITHM TO COMPUTE EFFECT FUNCTIONS AND THE 
#'          SMOOTHING TO PLOT        
#'
#' @return A list including
#'   \tabular{ll}{
#'     \emph{effnopar_tot} \tab Matrix including total effects in columns. \cr
#'     \emph{effnopar_dir} \tab Matrix including direct effects in columns. \cr
#'     \emph{effnopar_ind} \tab Matrix including indirect effects in columns. \cr
#'     \emph{effnopar_tot_up} \tab Matrix including upper bounds of total effects in columns. \cr
#'     \emph{effnopar_dir_up} \tab Matrix including upper bounds of direct effects in columns. \cr
#'     \emph{effnopar_ind_up} \tab Matrix including upper bounds of indirect effects in columns. \cr
#'     \emph{effnopar_tot_low} \tab Matrix including lower bounds of total effects in columns. \cr
#'     \emph{effnopar_dir_low} \tab Matrix including lower bounds of direct effects in columns. \cr
#'     \emph{effnopar_ind_low} \tab Matrix including lower bounds of indirect effects in columns. \cr
#'  }
#'         
#' @author Roman Minguez \email{roman.minguez@@uclm.es}
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{psar}} estimate spatial or spatio-temporal semiparametric PS-SAR
#'   regression models.
#'   \item \code{\link{eff_par}} compute and simulate total, direct and indirect effect
#'                                (or impacts) for parametric continuous covariates.
#'   \item \code{\link{fit_terms}} compute smooth functions for non-parametric
#'                                 continuous covariates.
#'   \item \code{\link{plot_effects_nopar}} plot the non-parametric effects functions
#'                                          allowing for previous smoothing.
#' }
#' 
#' @references \itemize{ 
#'   \item Basile, R., Durbán, M., Mínguez, R., Montero, J.
#'         M., and Mur, J. (2014). Modeling regional economic 
#'         dynamics: Spatial dependence, spatial heterogeneity and 
#'         nonlinearities. \emph{Journal of Economic Dynamics and 
#'         Control}, (48), 229-245.
#'         
#'   \item LeSage, J. and Pace, K. (2009). \emph{Introduction to 
#'         Spatial Econometrics}. CRC Press, Boca Raton.
#'         
#'   \item Mínguez, R.; Basile, R. and Durbán, M. (2018). An Alternative Semiparametric Model
#'         for Spatial Panel Data. Under evaluation in \emph{Statistical
#'         Methods and Applications}.
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
#'  ###### Non-Parametric Total, Direct and Indirect Effects
#'  list_varnopar <- c("serv","empgrowth")
#'  eff_nparvar <- eff_nopar(gamsar,list_varnopar)
#'  plot_effects_nopar(eff_nparvar,unemp_it,smooth=TRUE)
#'  plot_effects_nopar(eff_nparvar,unemp_it,smooth=FALSE)
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
#'  ###### Non-Parametric Total, Direct and Indirect Effects
#'  list_varnopar <- c("serv","empgrowth")
#'  eff_nparvar <- eff_nopar(geospanova_sar,list_varnopar)
#'  plot_effects_nopar(eff_nparvar,unemp_it,smooth=TRUE)
#'  plot_effects_nopar(eff_nparvar,unemp_it,smooth=FALSE)
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
#'  ###### Non-Parametric Total, Direct and Indirect Effects
#'  list_varnopar <- c("serv","empgrowth")
#'  eff_nparvar <- eff_nopar(sptanova_sar_ar1,list_varnopar)
#'  plot_effects_nopar(eff_nparvar,unemp_it,smooth=TRUE)
#'  plot_effects_nopar(eff_nparvar,unemp_it,smooth=FALSE)
#'
#' @keywords Indirect effects, Direct effects, SAR, non-parametric covariates.
#'
#' @export
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



