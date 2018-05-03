#' Estimate geoadditive spatial or spatio-temporal semiparametric PS-SAR
#' regression models.
#'
#' Estimate geoadditive spatial or spatio-temporal semiparametric PS-SAR
#' regression models including parametric and non-parametric covariates, spatial
#' or spatio-temporal non-parametric trends and spatial lags. The non-parametric
#' terms (either trends or covariates) are modelled using P-Splines. The
#' spatio-temporal non-parametric trend can be decomposed in an ANOVA way
#' including main effects and interactions of 2nd and 3rd order. For the spatial
#' lag a SAR model is specified. The estimation method
#' is REML (citar...)
#'
#' @param formula A formula similar to GAM specification. Parametric covariates
#'   are like linear models and non-parametric p-spline smooth terms are
#'   specified using \code{pspl(.)} and \code{pspt(.)} for the non-parametric covariates
#'   and spatial   or spatio-temporal trend, respectively.
#'   More details in \emph{Details} and \emph{Examples}.
#' @param data  A data frame containing the parametric and non-parametric
#'   covariates included in the model. Also, if a \code{pspt(.)} term is included
#'   in formula, the data frame must include the spatial and temporal
#'   coordinates specified in \code{pspt(.)}.
#' @param Wsp A neighbour spatial matrix. The dimension of the matrix is always
#'   \emph{NxN}, where \emph{N} is the dimension of each spatial coordinate.
#'   Default NULL.
#' @param sar A logical value used to include a spatial lag of the dependent
#'   variable in the model, that is, to estimate PS-SAR models. Default TRUE.
#' @param ar1 A logical value indicating if the noise of spatio-temporal model
#'   follows a first order autoregressive process. Default FALSE.
#' @param control List of extra control arguments - see section below
#'
#' @details Function to estimate the model
#'   \deqn{ y = (\rho*W_N \otimes I_T) y + X \beta +
#'     f(s_1,s_2,\tau_t) + \sum_{i=1}^k g(z_i) + \epsilon }
#' where \eqn{f(s_1,s_2,\tau_t)} is a smooth spatio-temporal trend
#' of the spatial coordinates \eqn{s1,s_2} and of the temporal
#' coordinate \eqn{\tau_t}.
#'
#' On the other hand, \eqn{g(z_i)} are also non-parametric smooth
#' functions of the covariates \eqn{z_i} and X includes the parametric
#' covariates; \eqn{W_N} is the spatial weights matrix,
#' \eqn{\rho} the spatial spillover parameter, \eqn{I_T} is identity matrix
#' or order \eqn{T} (equals to 1 for pure spatial data) and
#' \eqn{\epsilon ~ N(0,R)} where \eqn{R} can be multiple
#' of the identity (if errors are uncorrelated), or include an AR(1)
#' temporal autoregressive.
#' COMMENT HOW TO SPECIFY EACH TERM IN THE FORMULA
#'
#' In many situations the spatio-temporal trend to be estimated by
#' \eqn{f(s_1,s_2,\tau_t)} can be complex, and the use of a
#' multidimensional smooth function may not be flexible enough to capture the structure in the data.
#' To solve this problem, Lee and Durbán (2011) proposed an ANOVA-type
#' decomposition of this spatio-temporal trend where
#' spatial and temporal main effects, and second- and third-order
#' interactions between them can be identified:
#' \deqn{
#'   f(s_1,s_2,\tau_t) =
#'   f_1(s_1) + f_2(s_2) + f_t(\tau_t) +
#'   f_{1,2}(s_1,s_2) +  f_{1,t}(s_1,\tau_t) +
#'   f_{2,t}(s_2,\tau_t) + f_{1,2,t}(s_1,s_2,\tau_t) }
#'
#' @return  A list object of class \emph{psar}
#' @return \emph{call} Call of the function.
#' @return \emph{terms} The terms object used.
#' @return \emph{contrasts} (only where relevant) the contrasts used
#'              for parametric covariates.
#' @return \emph{xlevels}(only where relevant) a record of the levels
#'              of the parametric factors used in fitting.
#' @return \emph{edftot} Equivalent degrees of freedom for the whole model.
#' @return \emph{edfspt} Equivalent degrees of freedom for smooth
#'              spatio-temporal or spatial trend.
#' @return \emph{psanova} TRUE if spatio-temporal or spatial trend is PS-ANOVA.
#' @return \emph{edfnopar} Equivalent degrees of freedom for
#'              non-parametric covariates.
#' @return \emph{bfixed} Estimated betas corresponding to fixed effects in
#'              mixed model. These betas comes from either parametric
#'              covariates or fixed coefficients of smooth terms
#'              reparameterized as mixed models.
#' @return \emph{se_bfixed} Standard errors of fixed betas.
#' @return \emph{brandom} Estimated betas corresponding to random effects
#'              in mixed model. These betas comes from random coefficients of smooth
#'              terms reparameterized as mixed models.
#' @return \emph{cov_bfixed_brandom} Covariance matrix of fixed and random
#'              effects.
#' @return \emph{se_brandom} Standard errors of random betas.
#' @return \emph{sar} TRUE if model is PS-SAR.
#' @return \emph{rho} Estimated rho. If sar=FALSE always rho=0.
#' @return \emph{se_rho} Standard error of rho.
#' @return \emph{ar1} TRUE if model has a temporal autoregressive of
#'              first order in residuals.
#' @return \emph{phi} Estimated phi. If ar1=FALSE always phi=0.
#' @return \emph{se_rho} Standard error of phi.
#' @return \emph{fitted.values} Vector of fitted values of the dependent
#'              variable.
#' @return \emph{se_fitted.values} Vector of standard errors of fitted values
#'              values.
#' @return \emph{fit_Ay} Vector of fitted values of the spatial lag of
#'              dependent variable: \eqn{(\rho*W_N \otimes I_T) y}.
#' @return \emph{se_fitted.values} Vector of standard errors of fit_Ay.
#' @return \emph{residuals} Vector of residuals.
#' @return \emph{df.residual} Equivalent degrees of freedom for residuals.
#' @return \emph{sig2} Residual variance computed as SSR/df.residual.
#' @return \emph{residuals_norm} Uncorrelated residuals. For non-temporal
#'              data they are the same than \emph{residuals}.
#' @return \emph{llik} Log-likelihood value.
#' @return \emph{llik_reml} Restricted log-likelihood value.
#' @return \emph{aic} Akaike information criterion.
#' @return \emph{bic} Bayesian information criterion.
#' @return \emph{sp1} First spatial coordinate.
#' @return \emph{sp2} Second spatial coordinate.
#' @return \emph{time} Time coordinate.
#' @return \emph{y} Dependent variable.
#' @return \emph{X} Model matrix for fixed effects.
#' @return \emph{Z} Model matrix for random effects.
#'
#' @section Control arguments:
#' \itemize{
#' \item \emph{vary_init} A numerical
#'   value for the initial value of the noise variance in the model. Default
#'   var(y). Default NULL.
#'   \item \emph{trace} A logical value set to TRUE to provide information about
#'   the convergence of the estimation process. Default FALSE.
#'   \item \emph{thr}   Numerical value for the threshold of convergence in the estimation process.
#'   Default 1e-3.
#'   \item \emph{maxit} An integer value for the limit of the
#'   number of iterations during estimation process before reach convergence.
#'   Default 200.
#'   \item \emph{rho_init} An initial value for rho parameter if sar=TRUE. Default NULL.
#'   \item \emph{phi_init} An initial value for phi parameter if ar1=TRUE.
#'   }
#'
#' @author Roman Minguez \email{roman.minguez@@uclm.es}
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{eff__par}} compute total, direct and indirect effect
#'                                functions for parametric continuous covariates.
#'   \item \code{\link{eff_nopar}} compute total, direct and indirect effect
#'                                 functions for non-parametric continuous covariates.
#'   \item \code{\link{fit_terms}} compute smooth functions for non-parametric
#'                                 continuous covariates.
#'    \item\code{\link[mgcv]{gam}} well-known alternative of estimation
#'                                 of semiparametric models.
#' }
#'
#' @references \itemize{ \item Basile, R., Durbán, M., Mínguez, R., Montero, J.
#' M., and Mur, J. (2014). Modeling regional economic dynamics: Spatial
#' dependence, spatial heterogeneity and nonlinearities. \emph{Journal of
#' Economic Dynamics and Control}, (48), 229-245.
#'
#' \item Eilers, P. and Marx, B. (1996). Flexible Smoothing with B-Splines and
#' Penalties. \emph{Statistical Science}, (11), 89-121.
#'
#' \item Lee, D. and Durbán, M. (2011). P-Spline ANOVA Type Interaction Models
#' for Spatio-Temporal Smoothing. \emph{Statistical Modelling}, (11), 49-69.
#'
#' \item Lee, D. J., Durban, M., and Eilers, P. (2013). Efficient
#' two-dimensional smoothing with P-spline ANOVA mixed models and nested bases.
#' \emph{Computational Statistics & Data Analysis}, (61), 22-37.
#'
#' \item LeSage, J. and Pace, K. (2009). \emph{Introduction to Spatial
#' Econometrics}. CRC Press, Boca Raton.
#'
#' \item Montero, J., Mínguez, R., and Durbán, M. (2012). SAR models with
#' nonparametric spatial trends: A P-Spline approach. \emph{Estadística
#' Española}, (54:177), 89-111. }
#'
#' @examples
#' ################################################
#'  ###################### Examples using a panel data of rate of
#'  ###################### unemployment for 103 Italian provinces in period 1996-2014.
#'  library(sptpsar)
#'  data(unemp_it); Wsp <- Wsp_it
#'  ######################  GAM pure
#'  form1 <- unrate ~ partrate + agri + cons +
#'                    pspl(serv,nknots=15) +
#'                    pspl(empgrowth,nknots=20)
#'  gampure <- psar(form1,data=unemp_it)
#'  summary(gampure)
#'  library(mgcv)
#'
#'  form1_mgcv <- unrate ~ partrate + agri + cons +
#'                         s(serv,bs="ps",m=2,k=15) +
#'                         s(empgrowth,bs="ps",m=2,k=20)
#'  gampure_mgcv <- gam(form1_mgcv,data=unemp_it,method="REML")
#'  summary(gampure_mgcv)
#'  plot.gam(gampure_mgcv)
#'  ######################  Fit non-parametric terms (spatial trend must be name "spttrend")
#'  list_varnopar <- c("serv","empgrowth")
#'  terms_nopar <- fit_terms(gampure,list_varnopar)
#'  ######################  Plot non-parametric terms
#'  plot_terms(terms_nopar,unemp_it)
#' ######################  GAM SAR
#'  gamsar <- psar(form1,data=unemp_it,sar=TRUE,Wsp=Wsp_it)
#'  summary(gamsar)
#'  ###### Non-Parametric Total, Direct and Indirect Effects
#'  list_varnopar <- c("serv","empgrowth")
#'  eff_nparvar <- eff_nopar(gamsar,list_varnopar)
#'  plot_effects_nopar(eff_nparvar,unemp_it,smooth=TRUE)
#'  plot_effects_nopar(eff_nparvar,unemp_it,smooth=FALSE)
#'  ###################### Parametric Total, Direct and Indirect Effects
#'  list_varpar <- c("partrate","agri","cons")
#'  eff_parvar <- eff_par(gamsar,list_varpar)
#'  summary(eff_parvar)
#'
#'  ###############################################
#'  ### Spatial semiparametric model without spatial lag
#' form2 <- unrate ~ partrate + agri + cons +
#'                  pspl(serv,nknots=15) + pspl(empgrowth,nknots=20) +
#'                  pspt(long,lat,nknots=c(20,20),psanova=FALSE)
#'  ### Spatial trend fixed for 2014
#'   unemp_it_2d <- unemp_it[unemp_it$year==2014,]
#'   geosp1 <- psar(form2,data=unemp_it_2d)
#'   summary(geosp1)
#'  ### Spatial trend fixed for period 1996- 2014
#'  geosp2 <- psar(form2,data=unemp_it)
#'  summary(geosp2)
#'  ###############################################
#'  ### Spatial semiparametric model with spatial lag
#'  ### Spatial trend fixed for 2014
#'  geospar1 <- psar(form2,data=unemp_it_2d,Wsp=Wsp_it,sar=TRUE)
#'  summary(geospar1)
#'  ### Spatial trend fixed for period 1996-2014
#'  geospar2 <- psar(form2,data=unemp_it,Wsp=Wsp_it,sar=TRUE)
#'  summary(geospar2)
#'  #### Fitted Values and Residuals
#'  plot(geospar2$fit,unemp_it$unrate,xlab='fitted values',ylab="unrate",
#'       type="p",cex.lab=1.3,cex.main=1.3,
#'       main="Spatial semiparametric model with spatial lag",
#'       sub="Spatial trend fixed for period 1996-2014")

#'       plot(geospar2$fit,geospar2$residuals,xlab='fitted values',
#'            ylab="residuals",type="p",cex.lab=1.3,cex.main=1.3,
#'            main="Spatial semiparametric model with spatial lag",
#'            sub="Spatial trend fixed for period 1996-2014")
#'
#'  ### Fit non-parametric terms (spatial trend must be name "spttrend")
#' list_varnopar <- c("serv","empgrowth")
#' terms_nopar <- fit_terms(geospar2,list_varnopar)
#' ### Plot non-parametric terms
#'  plot_terms(terms_nopar,unemp_it)
#' ### Plot spatial trend (no ANOVA)
#' sptrend <- fit_terms(geospar2,"spttrend")
#' n <- dim(unemp_it)[1]
#' nT <- 19
#' sp_seq <- seq(from=1,to=n,by=nT)
#' lon <- scale(unemp_it[sp_seq,c("long")])
#' lat <- scale(unemp_it[sp_seq,c("lat")])
#' spt <- sptrend$fitted_terms[sp_seq]
#' require(rgl)
#' open3d(); plot3d(x=lon,y=lat,z=spt)
#' require(akima)
#' akima_spt <- akima::interp(x=lon, y=lat, z=spt,
#'              xo=seq(min(lon), max(lon), length = 1000),
#'              yo=seq(min(lat), max(lat), length = 1000))
#' persp3d(akima_spt$x,akima_spt$y,akima_spt$z,
#'         xlim=range(lon),ylim=range(lat),
#'         zlim=range(spt),col="green3",aspect="iso",add=TRUE)
#'  ###### Non-Parametric Total, Direct and Indirect Effects
#'  eff_nparvar <- eff_nopar(geospar2,list_varnopar)
#'  plot_effects_nopar(eff_nparvar,unemp_it,smooth=TRUE)
#'  plot_effects_nopar(eff_nparvar,unemp_it,smooth=FALSE)
#'  ###### Parametric Total, Direct and Indirect Effects
#'  list_varpar <- c("partrate","agri","cons")
#'  eff_parvar <- eff_par(geospar2,list_varpar)
#'  summary(eff_parvar)
#'
#'  ###############################################
#'  # Spatial semiparametric ANOVA model without spatial lag
#'  # Interaction term f12 with nested basis
 # form3 <- unrate ~ partrate + agri + cons +
 #                   pspl(serv,nknots=15) + pspl(empgrowth,nknots=20) +
 #                   pspt(long,lat,nknots=c(20,20),psanova=TRUE,
 #                   nest_sp1=c(1,2),nest_sp2=c(1,2))
#' # Spatial trend fixed for period 1996-2014
#' geospanova <- psar(form3,data=unemp_it)
#' summary(geospanova)
#' ### Plot spatial trend (ANOVA)
#' spttrend <- fit_terms(geospanova,"spttrend")
#' lon <- scale(unemp_it$long); lat <- scale(unemp_it$lat)
#' ### Plot main effects
#' plot_main_spt(spttrend,sp1=lon,sp2=lat,nT=19)
#' ### Plot Interaction
#' nT <- 19
#' sp_seq <- seq(from=1,to=n,by=nT)
#' f12_int <- spttrend$fitted_terms[sp_seq,c("f12_int")]
#' lon <- scale(unemp_it[sp_seq,c("long")])
#' lat <- scale(unemp_it[sp_seq,c("lat")])
#' require(rgl)
#' open3d(); plot3d(x=lon,y=lat,z=f12_int)
#' require(akima)
#' akima_f12 <- akima::interp(x=lon, y=lat, z=f12_int,
#'                       xo=seq(min(lon), max(lon), length = 1000),
#'                       yo=seq(min(lat), max(lat), length = 1000))
#' persp3d(akima_f12$x,akima_f12$y,akima_f12$z,
#'         xlim=range(lon),ylim=range(lat),
#'         zlim=range(f12_int),col="green3",aspect="iso",add=TRUE)

#' # Spatial semiparametric ANOVA model with spatial lag
#' # Interaction term f12 with nested basis
#' geospanova_sar <- psar(form3,data=unemp_it,Wsp=Wsp_it,sar=TRUE)
#' summary(geospanova_sar)

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
#'                   control=list(thr=1e-2,maxit=200,trace=TRUE))
#'  summary(sptanova)
#'  ### Plot spatial trend (ANOVA)
#'  spttrend <- fit_terms(sptanova,"spttrend")
#'  lon <- scale(unemp_it$long); lat <- scale(unemp_it$lat)
#'  time <- unemp_it$year
#'  ### Plot main effects
#'  plot_main_spt(spttrend,sp1=lon,sp2=lat,time=time,nT=19)
#'
#'  ### Spatio-temporal semiparametric ANOVA model with spatial lag
#'  sptanova_sar <- psar(form4,data=unemp_it,Wsp=Wsp_it,sar=TRUE,
#'                       control=list(thr=1e-1,maxit=100,trace=TRUE))
#'  summary(sptanova_sar)
#'  ### Spatio-temporal semiparametric ANOVA model with spatial lag
#'  ### and temporal autorregresive noise
#'  sptanova_sar_ar1 <- psar(form4,data=unemp_it,Wsp=Wsp_it,sar=TRUE,ar1=TRUE,
#'                     control=list(thr=1e-1,maxit=200,trace=TRUE))
#'  summary(sptanova_sar_ar1)
#'  ###### Non-Parametric Total, Direct and Indirect Effects
#'  list_varnopar <- c("serv","empgrowth")
#'  eff_nparvar <- eff_nopar(sptanova_sar_ar1,list_varnopar)
#'  plot_effects_nopar(eff_nparvar,unemp_it,smooth=TRUE)
#'  plot_effects_nopar(eff_nparvar,unemp_it,smooth=FALSE)
#'  ###### Parametric Total, Direct and Indirect Effects
#'  list_varpar <- c("partrate","agri","cons")
#'  eff_parvar <- eff_par(sptanova_sar_ar1,list_varpar)
#'  summary(eff_parvar)


#'  ###############################################
#'  # Spatio-temporal semiparametric ANOVA model without spatial lag
#'  # Interaction terms f12,f1t and f12t with nested basis
#'  # Interaction term f2t restricted to 0
#'   form5 <- unrate ~ partrate + agri + cons +
#'                   pspl(serv,nknots=15) + pspl(empgrowth,nknots=20) +
#'                   pspt(long,lat,year,nknots=c(18,18,6),psanova=TRUE,
#'                   nest_sp1=c(1,2,3),nest_sp2=c(1,2,3),
#'                   nest_time=c(1,2,2),f2t_int=FALSE,ntime=19)
#'  sptanova2 <- psar(form5,data=unemp_it,
#'                   control=list(thr=1e-2,maxit=200,trace=TRUE))
#'  summary(sptanova2)
#'  # Spatio-temporal semiparametric ANOVA model with spatial lag
#'  sptanova_sar2 <- psar(form5,data=unemp_it,Wsp=Wsp_it,sar=TRUE,
#'                       control=list(thr=1e-2,maxit=100,trace=TRUE))
#'  summary(sptanova_sar2)
#'  # Spatio-temporal semiparametric ANOVA model with spatial lag
#'  # and temporal autorregresive noise
#'  sptanova_sar_ar1_2 <- psar(form5,data=unemp_it,Wsp=Wsp_it,sar=TRUE,ar1=TRUE,
#'                     control=list(thr=1e-2,maxit=200,trace=TRUE))
#'  summary(sptanova_sar_ar1_2)
#'  # Spatio-temporal semiparametric ANOVA model without spatial lag
#'  # and temporal autorregresive noise
#'  sptanova_ar1_2 <- psar(form5,data=unemp_it,Wsp=Wsp_it,ar1=TRUE,
#'                     control=list(thr=1e-2,maxit=200,trace=TRUE))
#'  summary(sptanova_ar1_2)
#'  ###############################################
#'  # Spatio-temporal semiparametric model (no ANOVA) without spatial lag
#'  # Remark: Be careful with the number of knots
#'  #    (spend a lot of time fitting...)
#'   form6 <- unrate ~ partrate + agri + cons +
#'                   pspl(serv,nknots=15) + pspl(empgrowth,nknots=20) +
#'                   pspt(long,lat,year,nknots=c(10,10,6),ntime=19)
#'  spt <- psar(form6,data=unemp_it,Wsp=Wsp_it,
#'                     control=list(thr=1e-1,maxit=200,trace=TRUE))
#'  summary(spt)
#'
#'
#' @keywords P-Spline, SAR, Spatio-Temporal Trends, Semiparametric Models
#'
#' @export

psar <- function(formula,data,Wsp=NULL,ar=FALSE,ar1=FALSE,
                 control=list(vary_init=NULL,thr=1e-3,maxit=200,trace=TRUE,
                              rho_init=NULL,phi_init=NULL),...)
{
    if (any(grepl("pspt",formula))) formula <- update(formula, . ~ . - 1)
    cl <- match.call()
    mf <- match.call(expand.dots=TRUE)
    m <- match(c("formula","data",
                 #"subset", "weights", "na.action",
                 "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    mt <- terms(formula,specials=c("pspl","pspt"))
    names_var <- labels(mt)
    names_varspt <- names_var[grepl("pspt",names_var)]
    nvarspt <- length(names_varspt)
    names_varnopar <- names_var[grepl("pspl",names_var)]
    nvarnopar <- length(names_varnopar)
    names_varpar <- names_var[!grepl("pspl",names_var) & !grepl("pspt",names_var)]
    nvarpar <- length(names_varpar)
    #spec_bsp <- attr(mt,"specials")
    y <- model.response(mf, "numeric")
    nfull <- length(y)
    #X <-  model.matrix(mt, mf, contrasts)
    #Xpar <- X[,labels_varpar]
    Xfull <- Zfull <- cfull <- NULL
    pordfull <- bdegfull <- nknotsfull <- decomfull <- NULL
    dfull_list <- list()
    if (nvarpar>0) Xpar <- as.matrix(mf[,c(names_varpar)]) else Xpar <- NULL
    if (attr(mt,"intercept")==1){
      Intercept <- matrix(1,nrow=nfull,ncol=1)
      colnames(Intercept) <- c("Intercept")
      Xpar <- cbind(Intercept,Xpar)
      names_varpar <- c("Intercept",names_varpar)
      nvarpar <- length(names_varpar)
    }
    if (!is.null(Xpar)) Xfull <- cbind(Xfull,Xpar)

    if (nvarspt>0){
      #for (i in 1:length(names_varspt)) {
        varspt <- names_varspt
        Bi <- mf[,c(varspt)]
        sp1 <- attr(Bi,"sp1")
        sp2 <- attr(Bi,"sp2")
        nsp <- length(sp1)
        time <- attr(Bi,"time")
        ntime <- attr(Bi,"ntime")
        nknotsspt <- attr(Bi,"nknots")
        if (length(nknotsspt)==2) names(nknotsspt) <- c("sp1","sp2")
        if (length(nknotsspt)==3) names(nknotsspt) <- c("sp1","sp2","time")
        nknotsfull <- c(nknotsfull,nknotsspt)
        bdegspt <- attr(Bi,"bdeg")
        if (length(bdegspt)==2) names(bdegspt) <- c("sp1","sp2")
        if (length(bdegspt)==3) names(bdegspt) <- c("sp1","sp2","time")
        bdegfull <- c(bdegfull,bdegspt)
        pordspt <- attr(Bi,"pord")
        if (length(pordspt)==2) names(pordspt) <- c("sp1","sp2")
        if (length(pordspt)==3) names(pordspt) <- c("sp1","sp2","time")
        pordfull <- c(pordfull,pordspt)
        decomspt <- attr(Bi,"decom")
        names(decomspt) <- c("spt")
        decomfull <- c(decomfull,decomspt)
        psanova <- attr(Bi,"psanova")
        nest_sp1 <- attr(Bi,"nest_sp1")
        nest_sp2 <- attr(Bi,"nest_sp2")
        nest_time <- attr(Bi,"nest_time")
        f1_main <- attr(Bi,"f1_main")
        f2_main <- attr(Bi,"f2_main")
        ft_main <- attr(Bi,"ft_main")
        f12_int <- attr(Bi,"f12_int")
        f1t_int <- attr(Bi,"f1t_int")
        f2t_int <- attr(Bi,"f2t_int")
        f12t_int <- attr(Bi,"f12t_int")

        Bsptfull <- Bspt(sp1=sp1,sp2=sp2,
                          time=time,nfull=nfull,
                          ntime=ntime,psanova=psanova,Bi=Bi,bdegspt=bdegspt)
        sp1 <- Bsptfull$sp1; sp2 <- Bsptfull$sp2
        nsp <- length(sp1)
        time <- Bsptfull$time
        Bsptlist <- Bsptfull$Bsptlist
        rm(Bsptfull,Bi)
        XZsptlist <- B_XZ_spt(sp1=sp1,sp2=sp2,time=time,pordspt=pordspt,
                        psanova=psanova,decomspt=decomspt,
                        f1_main=f1_main,f2_main=f2_main,
                        ft_main=ft_main,f12_int=f12_int,f1t_int=f1t_int,
                        f2t_int=f2t_int,f12t_int=f12t_int,
                        Bsptlist=Bsptlist)
        Xspt <- XZsptlist$Xspt
        Zspt <- XZsptlist$Zspt
        dsptlist <- XZsptlist$dsptlist
        cspt <- unlist(XZsptlist$csptlist)
        Xfull <- cbind(Xfull,Xspt)
        Zfull <- cbind(Zfull,Zspt)
        dfull_list <- c(dfull_list,dsptlist)
        cfull <- c(cfull,cspt)
        #rm(Bsptlist,XZsptlist)
    } else {
      sp1 <- sp2 <- time <- NULL
      Xspt <- Zspt <- dsptlist <- cspt <- NULL
      nknotsspt <- pordspt <- bdegspt <- decomspt <- NULL
    } # if (!is.null(names_varspt))

    if (nvarnopar>0){
      Xnopar <- Znopar <-  cnopar <- NULL
      nknotsnopar <- bdegnopar <- pordnopar <- decomnopar <- NULL
      dnoparlist <- list()
      dnoparlistnames <- vector()
      for (i in 1:length(names_varnopar)) {
        varnopar <- names_varnopar[i]
        Bi <- mf[,c(varnopar)]
        colnames(Bi) <- paste(varnopar,1:ncol(Bi),sep=".")
        nknots_i <- attr(Bi,"nknots")
        pord_i <- attr(Bi,"pord")
        bdeg_i <- attr(Bi,"bdeg")
        decom_i <- attr(Bi,"decom")
        names(nknots_i) <- names(pord_i) <- varnopar
        names(bdeg_i) <- names(decom_i) <- varnopar
        nknotsnopar <- c(nknotsnopar,nknots_i)
        pordnopar <- c(pordnopar,pord_i)
        bdegnopar <- c(bdegnopar,bdeg_i)
        decomnopar <- c(decomnopar,decom_i)
        #attr(Bi,'varnopar') <- varnopar
        #Bnoparlist[[i]] <- Bi
        BtoXZ <- B_XZ(Bi)
        Xi <- as.matrix(BtoXZ$X[,-c(1)]) # Elimina intercepto
        colnames(Xi) <- paste(varnopar,1:ncol(Xi),sep=".")
        Zi <- BtoXZ$Z
        colnames(Zi) <- paste(varnopar,1:ncol(Zi),sep=".")
        dnoparlist[[i]] <- BtoXZ$d
        dnoparlistnames <- c(dnoparlistnames,
                             paste("d_",varnopar,sep=""))
        ci <- ncol(Bi)
        names(ci) <- paste(varnopar,sep="")
        Xnopar <- cbind(Xnopar,Xi)
        Znopar <- cbind(Znopar,Zi)
        cnopar <- c(cnopar,ci)
        rm(BtoXZ,Bi,Xi,Zi,ci,nknots_i,pord_i,bdeg_i,decom_i,varnopar)
      }
      names(dnoparlist) <- dnoparlistnames
      rm(dnoparlistnames)
      Xfull <- cbind(Xfull,Xnopar)
      Zfull <- cbind(Zfull,Znopar)
      dfull_list <- c(dfull_list,dnoparlist)
      cfull <- c(cfull,cnopar)
      nknotsfull <- c(nknotsfull,nknotsnopar)
      pordfull <- c(pordfull,pordnopar)
      bdegfull <- c(bdegfull,bdegnopar)
      decomfull <- c(decomfull,decomnopar)
    } else {
      Xnopar <- Znopar <- dnoparlist <- cnopar <- NULL
      nknotsnopar <- pordnopar <- bdegnopar <- decomnopar <- NULL
    }    # end if(!is.null(names_varnopar))
    thr <- control$thr
    maxit <- control$maxit
    trace <- control$trace
    vary_init <- control$vary_init
    if (is.null(vary_init)) vary_init <- var(y)
    cat("\nFitting Model...\n")
    if (!any(grepl("pspt",formula))) { # Without Spatial Trend
      model_fit <- fit_pspl2d_sar(y = y, vary_init = vary_init,
                                  sp1=NULL, sp2=NULL,
                                  Xfull = Xfull, Zfull = Zfull,
                                  Wsp = Wsp, nvarpar = nvarpar,
                                  nvarnopar = nvarnopar, nvarspt = nvarspt,
                                  #weights = NULL, GLAM = FALSE,
                                  cspt = NULL, dsptlist = NULL,
                                  bdegspt = NULL, pordspt = NULL,
                                  nknotsspt = NULL,
                                  cnopar = cnopar, dnoparlist = dnoparlist,
                                  bdegnopar = bdegnopar, pordnopar = pordnopar,
                                  nknotsnopar = nknotsnopar,
                                  names_varnopar = names_varnopar,
                                  names_varpar = names_varpar,
                                  sar = sar, rho_init = rho_init, rho_fixed = FALSE,
                                  bold = NULL, maxit = maxit, thr = thr,
                                  trace=trace)
    } else { # With Spatial or Spatio-Temporal Trend
      if (is.null(time)){ # With Spatial Trend
        model_fit <- fit_pspt2d_sar(y = y, vary_init = vary_init,
                                    sp1 = sp1, sp2 = sp2,
                                    Xfull = Xfull, Zfull = Zfull,Wsp = Wsp,
                                    nvarpar = nvarpar, nvarnopar = nvarnopar,
                                    nvarspt = nvarspt,
                                    #weights = NULL, GLAM = FALSE,
                                    cspt = cspt, dsptlist = dsptlist,
                                    bdegspt = bdegspt, pordspt = pordspt,
                                    nknotsspt = nknotsspt,
                                    cnopar = cnopar, dnoparlist = dnoparlist,
                                    bdegnopar = bdegnopar, pordnopar = pordnopar,
                                    nknotsnopar = nknotsnopar,
                                    names_varnopar = names_varnopar,
                                    names_varpar = names_varpar,
                                    psanova = psanova,
                                    f1_main = f1_main, f2_main = f2_main,
                                    f12_int = f12_int,
                                    sar = sar, rho_init = rho_init, rho_fixed = FALSE,
                                    bold = NULL, maxit = maxit, thr = thr,
                                    trace=trace)

      } else { # With Spatio-Temporal Trend
        model_fit <- fit_pspt3d_sar(y = y, vary_init = vary_init,
                                    sp1 = sp1, sp2 = sp2, time = time,
                                    Xfull = Xfull, Zfull = Zfull,Wsp = Wsp,
                                    nvarpar = nvarpar, nvarnopar = nvarnopar,
                                    nvarspt = nvarspt,
                                    #weights = NULL, GLAM = FALSE,
                                    cspt = cspt, dsptlist = dsptlist,
                                    bdegspt = bdegspt, pordspt = pordspt,
                                    nknotsspt = nknotsspt,
                                    cnopar = cnopar, dnoparlist = dnoparlist,
                                    bdegnopar = bdegnopar, pordnopar = pordnopar,
                                    nknotsnopar = nknotsnopar,
                                    names_varnopar = names_varnopar,
                                    names_varpar = names_varpar,
                                    psanova = psanova, f1_main = f1_main,
                                    f2_main = f2_main, ft_main = ft_main,
                                    f12_int = f12_int, f1t_int = f1t_int,
                                    f2t_int = f2t_int, f12t_int = f12t_int,
                                    sar = sar, rho_init = rho_init, rho_fixed = FALSE,
                                    ar1 = ar1, phi_init = phi_init, phi_fixed = FALSE,
                                    bold = NULL, maxit = maxit, thr = thr, trace=trace)
      }

    }
    mt_terms <- attr(mt,"term.labels")
    model_fit$contrasts <- attr(Xpar, "contrasts")
    model_fit$xlevels <- .getXlevels(mt, mf)
    model_fit$call <- cl
    model_fit$terms <- mt
    model_fit$model <- mf
    model_fit$X <- Xfull
    model_fit$Z <- Zfull
    model_fit$y <- y
    model_fit$df.residual <- length(y) - model_fit$edftot
    if (!is.null(Wsp)) model_fit$Wsp <- Wsp else  model_fit$Wsp <- NULL
    class(model_fit) <- c("psar","lm")
    model_fit
    # #z$durbin <- durbin
    # if(durbin)
    # {
    #     if(!is.null(Xpar)) Xpar <- cbind(Xpar,as.matrix(Wsp %*% Xpar))
    #     if(!is.null(Xnopar)) Xnopar <- cbind(Xnopar,as.matrix(Wsp %*% Xnopar))
    # }
}
