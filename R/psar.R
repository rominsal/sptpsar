#' @name psar
#' @rdname psar
#'
#' @title Estimate geoadditive spatial or spatio-temporal semiparametric PS-SAR
#' regression models.
#'
#' @description Estimate geoadditive spatial or spatio-temporal semiparametric PS-SAR
#' regression models including parametric and non-parametric covariates, spatial
#' or spatio-temporal non-parametric trends and spatial lags. The non-parametric
#' terms (either trends or covariates) are modelled using P-Splines. The
#' non-parametric trend can be decomposed in an ANOVA way including main
#' and interactions effects of 2nd and 3rd order. For the spatial
#' lag a SAR model is specified. The estimation method
#' is Restricted Maximum Likelihood (REML).
#'
#' @param formula A formula similar to GAM specification including
#'   parametric and non-parametric terms. Parametric covariates
#'   are included in the usual way and non-parametric p-spline smooth terms are
#'   specified using \code{pspl(.)} and \code{pspt(.)} for the non-parametric covariates
#'   and spatial   or spatio-temporal trend, respectively.
#'   More details in \emph{Details} and \emph{Examples}.
#' @param data  A data frame containing the parametric and non-parametric
#'   covariates included in the model. Also, if a \code{pspt(.)} term is included
#'   in formula, the data frame must include the spatial and temporal
#'   coordinates specified in \code{pspt(.)}. In this case coordinates
#'   must be ordered choosing time as fast index and spatial coordinates
#'   as low indexes. See \code{head(unemp_it)} as an example.
#' @param Wsp A neighbour spatial matrix. The dimension of the matrix is always
#'   \emph{NxN}, where \emph{N} is the dimension of each spatial coordinate.
#'   Default NULL.
#' @param sar A logical value used to include a spatial lag of the dependent
#'   variable in the model, that is, to estimate PS-SAR models. Default TRUE.
#' @param ar1 A logical value indicating if the noise of spatio-temporal model
#'   follows a first order autoregressive process. Default FALSE.
#' @param control List of extra control arguments - see section below
#'
#' @details
#' Function to estimate the model:
#'
#'   \deqn{ y = (\rho*W_N \otimes I_T) y + X \beta +
#'     f(s_1,s_2,\tau_t) + \sum_{i=1}^k g(z_i) + \epsilon }
#' where:
#' \itemize{
#'   \item \eqn{f(s_1,s_2,\tau_t)} is a smooth spatio-temporal trend
#' of the spatial coordinates \eqn{s1,s_2} and of the temporal
#' coordinate \eqn{\tau_t}.
#'   \item \eqn{X} is a matrix including values of parametric covariates.
#'   \item \eqn{g(z_i)} are non-parametric smooth functions of the
#'   covariates \eqn{z_i}.
#'   \item \eqn{W_N} is the spatial weights matrix.
#'   \item \eqn{\rho} is the spatial spillover parameter.
#'   \item \eqn{I_T} is an identity matrix of order \eqn{T} (\emph{T=1}
#'   for pure spatial data).
#'   \item \eqn{\epsilon ~ N(0,R)} where \eqn{ R = \sigma^2 I_T} if errors
#'    are uncorrelated or it follows an AR(1) temporal autoregressive structure
#'    for serially correlated errors.
#' }
#'
#' The non-parametric terms are included in \code{formula} using
#' \code{pspt(.)} for spatial or spatio-temporal trends and \code{pspl(.)}
#' for other non-parametric smooth additive terms.
#'
#' For example, if a model includes a spatio-temporal trend with \emph{long},
#' and \emph{lat} as spatial coordinates and \emph{year} as temporal coordinate;
#' two  non-parametric covariates (\emph{empgrowth} and \emph{serv}) and
#' three parametric covariates (\emph{partrate}, \emph{agri} and
#' emph{cons}), the formula (choosing default values
#' for each term) should be written as:
#'
#' \code{ unrate ~ partrate + agri + cons +
#'                    pspl(serv) + pspl(empgrowth) +
#'                    pspt(long,lat,year) }
#'
#'For spatial trend case the term \code{pspt(.)} does not include a temporal coordinate,
#'that is, in the previous example would be specified as \code{pspt(long,lat)}.
#'
#' In many situations the  spatio-temporal trend given by \eqn{f(s_1,s_2,\tau_t)}
#' can be very complex and the use of a multidimensional smooth function may not
#' be flexible enough to capture
#' the structure in the data. Furthermore, the estimation of this trend
#' can become computationally intensive especially for large databases.
#' To solve this problem, Lee and Durbán (2011) proposed an ANOVA-type
#' decomposition of this spatio-temporal trend where
#' spatial and temporal main effects, and second- and third-order
#' interaction effects can be identified as:
#'
#' \deqn{
#'   f(s_1,s_2,\tau_t) =
#'   f_1(s_1) + f_2(s_2) + f_t(\tau_t) +
#'   f_{1,2}(s_1,s_2) +  f_{1,t}(s_1,\tau_t) +
#'   f_{2,t}(s_2,\tau_t) + f_{1,2,t}(s_1,s_2,\tau_t) }
#'
#'   In this equation \eqn{f_1(s_1), f_2(s_2)} and
#'   \eqn{f_t(\tau_t)} are the main effects,
#'   \eqn{f_{1,2}(s_1,s_2), f_{1,t}(s_1,\tau_t)} and
#'   \eqn{f_{2,t}(s_2,\tau_t)} are the
#'   second-order interaction effects and
#'   \eqn{f_{1,2,t}(s_1,s_2,\tau_t)}
#'   is the third-order interaction effect. In this case, each effect
#'   can have its own degree of smoothing allowing a greater flexibility
#'   for the spatio-temporal trend. The ANOVA decomposition of the trend
#'   can be set as an argument in \code{pspl(.)} terms choosing
#'   \code{psanova=TRUE}. For example to choose an ANOVA decomposition in
#'   previous case we can set
#'   \code{pspt(long,lat,year,nknots=c(18,18,8),psanova=TRUE,ntime=19)}.
#'   Note that for spatio-temporal data it is needed to specify the number
#'   of temporal periods in \code{ntime} but the number of spatial
#'   coordinates are not needed. See \emph{Examples} for more details.
#'
#'   In most empirical cases main effects are more flexible than interaction
#'   effects and therefore, the number of knots in B-Spline basis for
#'   interaction effects do not need to be as big as the number of knots for
#'   main effects. \emph{Lee et al.}, (2013) proposed a nested basis procedure
#'   in which the number of knots for the interaction effects are reduced using
#'   \emph{divisors} such that the space spanned by B-spline bases used for
#'   interaction effects are a subset of the space spanned by B-spline bases
#'   used for main effects. The \emph{divisors} can be specified as an argument
#'   in \code{pspt(.)} terms. See \emph{Examples} for details.
#'
#'   All the terms are jointly fitted using Separation of Anisotropic Penalties
#'   (SAP) algorithm (see \emph{Rodríguez-Álvarez et al., (2015)}) which allows
#'   to the mixed model reparameterization of the model. When \emph{sar} or
#'   \emph{ar1} arguments are \emph{TRUE}, \eqn{\rho} and \eqn{\phi} parameters
#'   are numerically estimated using function \code{\link[nloptr]{nloptr}}
#'   implemented in pakage \pkg{nloptr} . In these cases, an iterative process between
#'   SAP and numerical optimization of \eqn{\rho} and \eqn{\phi} is applied until
#'   convergence. See details in \emph{Mínguez et al.}, (2018).
#'
#' @return
#' A list object of class \emph{psar}
#' \tabular{ll}{
#'  \code{call} \tab Call of the function. \cr
#'  \code{terms} \tab The terms object used. \cr
#'  \code{contrasts} \tab (only where relevant) the contrasts used
#'              for parametric covariates. \cr
#'  \code{xlevels} \tab (only where relevant) a record of the levels
#'              of the parametric factors used in fitting. \cr
#'  \code{edftot} \tab Equivalent degrees of freedom for the whole model. \cr
#'  \code{edfspt} \tab Equivalent degrees of freedom for smooth
#'              spatio-temporal or spatial trend. \cr
#'  \code{edfnopar} \tab Equivalent degrees of freedom for
#'              non-parametric covariates. \cr
#'  \code{psanova} \tab \emph{TRUE} if spatio-temporal or spatial trend is PS-ANOVA. \cr
#'  \code{bfixed} \tab Estimated betas corresponding to fixed effects in
#'              mixed model. These betas comes from either parametric
#'              covariates or fixed coefficients of smooth terms
#'              reparameterized as mixed models. \cr
#'  \code{se_bfixed} \tab Standard errors of fixed betas. \cr
#'  \code{brandom} \tab Estimated betas corresponding to random effects
#'              in mixed model. These betas comes from random coefficients of smooth
#'              terms reparameterized as mixed models. \cr
#'  \code{se_brandom}\tab Standard errors of random betas. \cr
#'  \code{vcov_b} \tab Covariance matrix of fixed and random
#'              effects. \cr
#'  \code{sar} \tab \emph{TRUE} if model is PS-SAR. \cr
#'  \code{rho} \tab Estimated rho. If \code{sar=FALSE} always \eqn{rho=0}. \cr \cr
#'  \code{se_rho} \tab Standard error of \eqn{rho}. \cr
#'  \code{ar1} \tab \emph{TRUE} if model has a temporal autoregressive of
#'              first order in residuals. \cr
#'  \code{phi} \tab Estimated phi. If \code{ar1=FALSE} always \eqn{phi=0}. \cr
#'  \code{se_rho} \tab Standard error of \eqn{phi}. \cr
#'  \code{fitted.values} \tab Vector of fitted values of the dependent
#'              variable. \cr
#'  \code{se_fitted.values} \tab Vector of standard errors of
#'       \code{fitted.values}. \cr
#'  \code{fitted.values_Ay} \tab Vector of fitted values of the spatial lag of
#'      dependent variable: \eqn{(\rho*W_N \otimes I_T) y}. \cr
#'  \code{se_fitted.values_Ay} \tab Vector of standard errors of
#'       \code{fitted.values_Ay}. \cr
#'  \code{residuals} \tab Vector of residuals. \cr
#'  \code{residuals_norm} \tab Uncorrelated residuals. For non-temporal
#'              data they are the same than \code{residuals}. \cr
#'  \code{df.residual} \tab Equivalent degrees of freedom for \code{residuals}. \cr
#'  \code{sig2}  \tab Residual variance computed as SSR/df.residual. \cr
#'  \code{llik} \tab Log-likelihood value. \cr
#'  \code{llik_reml} \tab Restricted log-likelihood value. \cr
#'  \code{aic} \tab Akaike information criterion. \cr
#'  \code{bic} \tab Bayesian information criterion. \cr
#'  \code{sp1} \tab First spatial coordinate. \cr
#'  \code{sp2} \tab Second spatial coordinate. \cr
#'  \code{time} \tab Time coordinate. \cr
#'  \code{y} \tab Dependent variable. \cr
#'  \code{X} \tab Model matrix for fixed effects. \cr
#'  \code{Z} \tab Model matrix for random effects. \cr
#' }
#'
#' @section Control arguments:
#' \tabular{ll}{
#'   \code{vary_init} \tab Initial value of the noise variance in the model.
#'     Default NULL.\cr
#'   \code{trace} \tab A logical value set to \emph{TRUE} to provide information
#'     about the convergence of the estimation process. Default \emph{FALSE}. \cr
#'   \code{thr} \tab Numerical value for the threshold of convergence
#'     in the estimation process. Default 1e-3. \cr
#'   \code{maxit} \tab An integer value for the limit of the
#'     number of iterations during estimation process before reach convergence.
#'     Default 200. \cr
#'   \code{rho_init} \tab An initial value for \eqn{rho} parameter if
#'     \code{sar=TRUE}. Default 0. \cr
#'   \code{phi_init} \tab An initial value for \eqn{phi} parameter if
#'     \code{ar1=TRUE}. Default 0. \cr
#'   \code{rho_fixed} \tab A logical value to set fixed \eqn{rho}
#'          parameter during estimation process. Default \emph{FALSE}. \cr
#'   \code{phi_fixed} \tab A logical value to set fixed \eqn{phi}
#'          parameter during estimation process. Default \emph{FALSE}. \cr
#' }
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
#'                                 of semiparametric models in \pkg{mgcv} package.
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
#' Mínguez, R.; Basile, R. and Durbán, M. (2018). An Alternative Semiparametric Model
#' for Spatial Panel Data. Under evaluation
#' in \emph{Regional Science and Urban Economics}.
#'
#' \item Montero, J., Mínguez, R., and Durbán, M. (2012). SAR models with
#' nonparametric spatial trends: A P-Spline approach. \emph{Estadística
#' Española}, (54:177), 89-111.
#'
#' \item Rodríguez-Alvarez, M. X.; Kneib, T.; Durban, M.; Lee, D.J.
#' and Eilers, P. (2015). Fast smoothing parameter separation in
#' multidimensional generalized P-splines: the SAP algorithm.
#' \emph{Statistics and Computing} 25 (5), 941-957.
#' }
#'
#' @examples
#' ################################################
#'  ###################### Examples using a panel data of rate of
#'  ###################### unemployment for 103 Italian provinces in period 1996-2014.
#' library(sptpsar)
#' data(unemp_it); Wsp <- Wsp_it
#'  ######################  GAM pure
#' form1 <- unrate ~ partrate + agri + cons +
#'                  pspl(serv,nknots=15) +
#'                  pspl(empgrowth,nknots=20)
#' gampure <- psar(form1,data=unemp_it)
#' summary(gampure)
#' library(mgcv)
#' form1_mgcv <- unrate ~ partrate + agri + cons +
#'                        s(serv,bs="ps",m=2,k=15) +
#'                        s(empgrowth,bs="ps",m=2,k=20)
#' gampure_mgcv <- gam(form1_mgcv,data=unemp_it,method="REML")
#' summary(gampure_mgcv)
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
#'                 pspl(serv,nknots=15) + pspl(empgrowth,nknots=20) +
#'                 pspt(long,lat,nknots=c(20,20),psanova=FALSE)
#'  ### Spatial trend fixed for 2014
#'  unemp_it_2d <- unemp_it[unemp_it$year==2014,]
#'  geosp1 <- psar(form2,data=unemp_it_2d)
#'  summary(geosp1)
#'  ### Spatial trend fixed for period 1996- 2014
#' geosp2 <- psar(form2,data=unemp_it)
#' summary(geosp2)
#'  ###############################################
#'  ### Spatial semiparametric model with spatial lag
#'  ### Spatial trend fixed for 2014
#' geospar1 <- psar(form2,data=unemp_it_2d,Wsp=Wsp_it,sar=TRUE)
#' summary(geospar1)
#'  ### Spatial trend fixed for period 1996-2014
#' geospar2 <- psar(form2,data=unemp_it,Wsp=Wsp_it,sar=TRUE)
#' summary(geospar2)
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
#'                    control=list(thr=1e-2,maxit=200,trace=TRUE))
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

psar <- function(formula, data, Wsp = NULL, sar = FALSE,
                 ar1 = FALSE,
                 control=list( vary_init = NULL,
                               thr = 1e-3, maxit = 200,
                               trace = FALSE, rho_init = 0,
                               phi_init = 0, rho_fixed = FALSE,
                               phi_fixed = FALSE), ...)
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
    rho_init <- control$rho_init
    phi_init <- control$phi_init
    rho_fixed <- control$rho_fixed
    phi_fixed <- control$phi_fixed
    if (is.null(trace)) trace <- FALSE
    if (is.null(maxit)) maxit <- 200
    if (is.null(thr)) thr <- 1e-3
    if (is.null(vary_init)) vary_init <- var(y)
    if (is.null(rho_init)) rho_init <- 0
    if (is.null(phi_init)) phi_init <- 0
    if (is.null(rho_fixed)) rho_fixed <- FALSE
    if (is.null(phi_fixed)) phi_fixed <- FALSE
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
                                  sar = sar, rho_init = rho_init,
                                  rho_fixed = rho_fixed,
                                  bold = NULL, maxit = maxit, thr = thr,
                                  trace = trace)
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
                                    sar = sar, rho_init = rho_init,
                                    rho_fixed = rho_fixed,
                                    bold = NULL, maxit = maxit, thr = thr,
                                    trace = trace)

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
                                    sar = sar, rho_init = rho_init,
                                    rho_fixed = rho_fixed,
                                    ar1 = ar1, phi_init = phi_init,
                                    phi_fixed = phi_fixed,
                                    bold = NULL, maxit = maxit, thr = thr,
                                    trace = trace)
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
