#' @name plot_timetrend
#' @rdname plot_timetrend
#'
#' @title Plot time trend for each spatial unit.
#' 
#' @description In models including a spatio-temporal trend plot a time
#'   trend for each spatial unit in the dataframe. 
#'   It can be useful to detect heterogeneous pattern of time trends for each
#'   spatial unit.   
#'
#' @param sptsarfit \emph{psar} object fitted using \code{\link{psar}} function.
#' @param data dataframe with the data.
#' @param time name of variable including the temporal coordinate. Default "year".
#' @param spat name of variable including the spatial units. Default "region".
#' @param xlab label of x-axis. Default "Time".
#' @param ylab label of y-axis. Default "Time Trend".
#' @param title title of graphics. Default "Non-Parametric Time Trend by Spatial Unit"
#'                                 
#' @author Roman Minguez \email{roman.minguez@@uclm.es}
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{plot_main_spt}} plot main terms of non-parametric 
#'     spatial (2d) or spatio-temporal (3d) trends in ANOVA models.
#'   \item \code{\link[mgcv]{plot.gam}} plot the terms fitted by 
#'     \code{\link[mgcv]{gam}} function in \pkg{mgcv} package.   
#' }
#'   
#' @examples
#' ################################################
#'  ###################### Examples using a panel data of rate of
#'  ###################### unemployment for 103 Italian provinces in period 1996-2014.
#' library(sptpsar)
#' data(unemp_it); Wsp <- Wsp_it
#'  ###############################################
#'  # Spatio-temporal semiparametric ANOVA model without spatial lag
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
#'  plot_timetrend(sptanova,data=unemp_it,time="year",
#'                spat="name")
#'                
#' @export

plot_timetrend <- function(sptsarfit,data,
                           time="year",spat="region",
                           xlab="Time",ylab="Time Trend",
                           title="Non-Parametric Time Trend by Spatial Unit"){
  terms_spt_trend <- fit_terms(sptsarfit,"spttrend")
  spt_trend <- terms_spt_trend$fitted_terms
  se_spt_trend <- terms_spt_trend$se_fitted_terms
  data$spttrend <- spt_trend[,c("spttrend")]
  data$se_spttrend <- se_spt_trend[,c("spttrend")]
  data$spat_var <- as.factor(data[,colnames(data)==spat])
  data$time_var <- as.factor(data[,colnames(data)==time])
  spttrend_by_year <- NULL
  for (i in 1:length(levels(data$spat_var)))
  {
    name_reg <- levels(data$spat_var)[i]
    spttrend_reg <- data$spttrend[data$spat_var==name_reg]
    spttrend_by_year <- cbind(spttrend_by_year,spttrend_reg)
  }
  colnames(spttrend_by_year) <- levels(data$spat_var)
  matplot(levels(data$time_var),spttrend_by_year,type="l",
          xlab=xlab,ylab=ylab,
          main=title)
}
