#' @name plot_eff_nopar
#' @rdname plot_eff_nopar
#'
#' @title Plot direct, indirect and total effects functions 
#'   for continous non-parametric covariates in PS-SAR regression models.
#'        
#' @description Plot direct, indirect and total effect functions for 
#'   non-parametric covariates included in a semiparametric spatial
#'   or spatio-temporal SAR model. This model must include a spatial
#'   lag of the dependent variable (SAR) to have indirect effects 
#'   different from 0, otherwise, total and direct function effects 
#'   are the same. The effect functions can be smoothed to overcome 
#'   the instabilities created by the premultiplication of matrix
#'   \eqn{(I - \rho W)^{-1}} 
#'
#' @param effnopar object returned from \code{\link{eff_nopar}} function.
#' @param data dataframe with the data. 
#' @param smooth logical value to choose smoothing of the effects function
#'               prior to plot. Default TRUE.
#' @param span span for the kernel of the smoothing (see \code{\link{loess}} 
#'             for details). Default c(0.1,0.1,0.2). 
#'
#' @return plot of the direct, indirect and total effects function for each non-parametric
#'   covariate included in the object returned from \code{\link{effects_nopar}}.
#'                                 
#' @author Roman Minguez \email{roman.minguez@@uclm.es}
#'
#' @family Direct, Indirect and Total Effects.
#' 
#' @seealso
#' \itemize{
#'   \item \code{\link{eff_nopar}} compute total, direct and indirect effect
#'           functions for non-parametric continuous covariates.
#'   \item \code{\link{fit_terms}} compute smooth functions for non-parametric
#'                                 continuous covariates.
#'   \item \code{\link{plot_terms}} plot the terms of non-parametric covariates.
#' }
#' 
#' @references \itemize{ 
#'   \item Basile, R., Durbán, M., Mínguez, R., Montero, J.
#'         M., and Mur, J. (2014). Modeling regional economic 
#'         dynamics: Spatial dependence, spatial heterogeneity and 
#'         nonlinearities. \emph{Journal of Economic Dynamics and 
#'         Control}, (48), 229-245.                  
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
#'  #' ######################   PSAR-ANOVA with spatial trend
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
#'  @export

plot_eff_nopar <- function(effnopar,data,smooth=TRUE,
                               span=c(0.1,0.1,0.2)){
  tot <- effnopar$effnopar_tot
  uptot <- effnopar$effnopar_tot_up
  lowtot <- effnopar$effnopar_tot_low
  dir <- effnopar$effnopar_dir
  updir <- effnopar$effnopar_dir_up
  lowdir <- effnopar$effnopar_dir_low
  ind <- effnopar$effnopar_ind
  upind <- effnopar$effnopar_ind_up
  lowind <- effnopar$effnopar_ind_low

  for(i in 1:ncol(tot))
  {
    name_var <- colnames(tot)[i]
    var <- matrix(data[,c(name_var)],ncol=1)
    colnames(var) <- name_var
    ord <- order(var)

    tot_i <- matrix(tot[,c(name_var)],ncol=1)
    colnames(tot_i) <- name_var
    uptot_i <- matrix(uptot[,c(name_var)],ncol=1)
    colnames(uptot_i) <- name_var
    lowtot_i <- matrix(lowtot[,c(name_var)],ncol=1)
    colnames(lowtot_i) <- name_var
    dir_i <- matrix(dir[,c(name_var)],ncol=1)
    colnames(dir_i) <- name_var
    updir_i <- matrix(updir[,c(name_var)],ncol=1)
    colnames(updir_i) <- name_var
    lowdir_i <- matrix(lowdir[,c(name_var)],ncol=1)
    colnames(lowdir_i) <- name_var
    ind_i <- matrix(ind[,c(name_var)],ncol=1)
    colnames(ind_i) <- name_var
    upind_i <- matrix(upind[,c(name_var)],ncol=1)
    colnames(upind_i) <- name_var
    lowind_i <- matrix(lowind[,c(name_var)],ncol=1)
    colnames(lowind_i) <- name_var
    if(smooth){
      span_tot <- span[1]; span_dir <- span[2]; span_ind <- span[3]
      tot_i_smooth <- predict(loess(tot_i~var, span = span_tot),
                            method = "loess()")
      uptot_i_smooth <- predict(loess(uptot_i~var, span = span_tot),
                              method = "loess()")
      lowtot_i_smooth <- predict(loess(lowtot_i~var, span = span_tot),
                                method = "loess()")
      tot_i <- tot_i_smooth
      uptot_i <- uptot_i_smooth
      lowtot_i <- lowtot_i_smooth

      dir_i_smooth <- predict(loess(dir_i~var, span = span_dir),
                              method = "loess()")
      updir_i_smooth <- predict(loess(updir_i~var, span = span_dir),
                                method = "loess()")
      lowdir_i_smooth <- predict(loess(lowdir_i~var, span = span_dir),
                                 method = "loess()")
      dir_i <- dir_i_smooth
      updir_i <- updir_i_smooth
      lowdir_i <- lowdir_i_smooth

      ind_i_smooth <- predict(loess(ind_i~var, span = span_ind),
                              method = "loess()")
      upind_i_smooth <- predict(loess(upind_i~var, span = span_ind),
                                method = "loess()")
      lowind_i_smooth <- predict(loess(lowind_i~var, span = span_ind),
                                 method = "loess()")
      ind_i <- ind_i_smooth
      upind_i <- upind_i_smooth
      lowind_i <- lowind_i_smooth
      ##### SET TOTAL SMOOTH AS THE ADDTION OF SMOOTH DIRECT AND INDIRECT
      # tot_i <- dir_i + ind_i
      # uptot_i <- updir_i + upind_i
      # lowtot_i <- lowdir_i + lowind_i
    }
    par(mfrow=c(3,1))
    plot(var[ord],tot_i[ord],type="l",
         ylab=paste("f(",name_var,")"),xlab=name_var,
         ylim=c(min(lowtot_i),max(uptot_i)),cex.lab=1.0,col=2,lty=1,lwd=2,
         cex.main=1.0,main="Total Effects")
    lines(var[ord],uptot_i[ord],xlab="",ylab="",type="l",col=2,lty=2,lwd=1.5)
    lines(var[ord],lowtot_i[ord],xlab="",ylab="",type="l",col=2,lty=2,lwd=1.5)
    abline(a=0,b=0)

    plot(var[ord],dir_i[ord],type="l",
         ylab=paste("f(",name_var,")"),xlab=name_var,
         ylim=c(min(lowdir_i),max(updir_i)),cex.lab=1.0,col=3,lty=1,lwd=2,
         cex.main=1.0,main="Direct Effects")
    lines(var[ord],updir_i[ord],xlab="",ylab="",type="l",col=3,lty=2,lwd=1.5)
    lines(var[ord],lowdir_i[ord],xlab="",ylab="",type="l",col=3,lty=2,lwd=1.5)
    abline(a=0,b=0)

    plot(var[ord],ind_i[ord],type="l",
         ylab=paste("f(",name_var,")"),xlab=name_var,
         ylim=c(min(lowind_i),max(upind_i)),cex.lab=1.0,col=4,lty=1,lwd=2,
         cex.main=1.0,main="Indirect Effects")
    lines(var[ord],upind_i[ord],xlab="",ylab="",type="l",col=4,lty=2,lwd=1.5)
    lines(var[ord],lowind_i[ord],xlab="",ylab="",type="l",col=4,lty=2,lwd=1.5)
    abline(a=0,b=0)
    readline(prompt="Press [enter] to continue")
  }
  par(mfrow=c(1,1))
}
