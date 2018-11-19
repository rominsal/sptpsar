#' @name plot_main_spt
#' @rdname plot_main_spt
#'
#' @title Plot main terms in ANOVA spatial or spatio-temporal trends.
#'        
#' @description Plot main terms for spatial and temporal coordinates 
#'   in spatial (2d) or spatio-temporal (3d) trends decomposed in ANOVA
#'   way.   
#'
#' @param spttrend object returned from \code{\link{fit_terms}} function
#'   including \emph{spttrend} in the \emph{variables} argument.
#' @param sp1 vector of first spatial coordinate. 
#' @param sp2 vector of second spatial coordinate.
#' @param time vector of temporal coordinate. It is NULL in spatial (2d)
#'   trends. Default NULL.
#' @param nT Number of time periods (1 for non-temporal data). Default 1.    
#' @param conflevel numerical value for the confidence interval of the
#'   trend functions. Default 0.95.
#'   
#' @return plot of each main trend (spatial and temporal) in ANOVA models.
#'                                 
#' @author Roman Minguez \email{roman.minguez@@uclm.es}
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{fit_terms}} compute the terms for non-parametric 
#'     trend and smooth functions for non-parametric continuous covariates .
#'   \item \code{\link[mgcv]{plot.gam}} plot the terms fitted by 
#'     \code{\link[mgcv]{gam}} function in \pkg{mgcv} package.   
#' }
#' 
#' @references 
#' \itemize{ 
#'   \item Lee, D. and Durbán, M. (2011). P-Spline ANOVA Type Interaction 
#'         Models for Spatio-Temporal Smoothing. \emph{Statistical Modelling}, (11), 49-69.
#'  }
#'         
#' @examples
#' ################################################
#'  ###################### Examples using a panel data of rate of
#'  ###################### unemployment for 103 Italian provinces in period 1996-2014.
#' library(sptpsar)
#' data(unemp_it); Wsp <- Wsp_it
#' ###############################################
#'  # Spatial (2d) semiparametric ANOVA model without spatial lag
#'  # Interaction term f12 with nested basis
#' form3 <- unrate ~ partrate + agri + cons +
#'                   pspl(serv,nknots=15) + pspl(empgrowth,nknots=20) +
#'                   pspt(long,lat,nknots=c(20,20),psanova=TRUE,
#'                   nest_sp1=c(1,2),nest_sp2=c(1,2))
#' # Spatial trend fixed for period 1996-2014
#' geospanova <- psar(form3,data=unemp_it)
#' summary(geospanova)
#' ### Plot spatial trend (ANOVA)
#' spttrend <- fit_terms(geospanova,"spttrend")
#' lon <- scale(unemp_it$long); lat <- scale(unemp_it$lat)
#' ### Plot main effects
#' plot_main_spt(spttrend,sp1=lon,sp2=lat,nT=19)
#' 
#' #'  ###############################################
#'  # Spatio-temporal (3d) semiparametric ANOVA model without spatial lag
#'  # Interaction terms f12,f1t,f2t and f12t with nested basis
#'  # Remark: It is necessary to include ntime as argument
#'  # Remark: nest_sp1, nest_sp2 and nest_time must be divisors of nknots
#'  form4 <- unrate ~ partrate + agri + cons +
#'                    pspl(serv,nknots=15) + pspl(empgrowth,nknots=20) +
#'                    pspt(long,lat,year,nknots=c(18,18,8),psanova=TRUE,
#'                    nest_sp1=c(1,2,3),nest_sp2=c(1,2,3),
#'                    nest_time=c(1,2,2),ntime=19)
#'  sptanova <- psar(form4,data=unemp_it,
#'                   control=list(thr=1e-2,maxit=200,trace=FALSE))
#'  summary(sptanova)
#'  ### Plot spatial trend (ANOVA)
#'  spttrend <- fit_terms(sptanova,"spttrend")
#'  lon <- scale(unemp_it$long); lat <- scale(unemp_it$lat)
#'  time <- unemp_it$year
#'  ### Plot main effects
#'  plot_main_spt(spttrend,sp1=lon,sp2=lat,time=time,nT=19)
#' 
#' @export

plot_main_spt <- function(spttrend,sp1,sp2,nT=1,time=NULL,conflevel=0.95){
  fit_spt <- spttrend$fitted_terms
  se_fit_spt <- spttrend$se_fitted_terms
  n <- nrow(fit_spt)
  if (is.null(time)){ # 2d
    main_effects <- c("f1_main","f2_main")
    par(mfrow=c(2,1))
  } else {
    main_effects <- c("f1_main","f2_main","ft_main")
    par(mfrow=c(3,1))
  }
  crval <- qnorm((1-conflevel)/2,mean=0,sd=1,lower.tail=FALSE)
  for (i in 1:length(main_effects)){
    main <- main_effects[i]
    if (main!="ft_main"){
      seq_sel <- seq(from=1,to=n,by=nT)
    } else {
      seq_sel <- seq(from=1,to=nT,by=1)
    }
    fmain <- fit_spt[seq_sel,main]
    se_fmain <- se_fit_spt[seq_sel,main]
    up_fmain <- fmain + crval*se_fmain
    low_fmain <- fmain - crval*se_fmain
    if (main=="f1_main") {
      coord <- sp1[seq_sel]
      names(coord) <- "sp1"
    } else if (main=="f2_main") {
      coord <- sp2[seq_sel]
      names(coord) <- "sp2"
    } else {
      coord <- time[seq_sel]
      names(coord) <- "time"
    }
    ord <- order(coord)
    plot(coord[ord],fmain[ord],type="l",
         ylab=main,xlab=names(coord),
         ylim=c(min(low_fmain),max(up_fmain)),
         cex.lab=1.0,col=2,lty=1,lwd=2,cex.main=1.0,
         main=paste("Spatio-Temporal ANOVA Effect: ",main))
    lines(coord[ord],up_fmain[ord],xlab="",ylab="",type="l",col=2,lty=2,lwd=1.5)
    lines(coord[ord],low_fmain[ord],xlab="",ylab="",type="l",col=2,lty=2,lwd=1.5)
    abline(a=0,b=0)
  }
  par(mfrow=c(1,1))
}
