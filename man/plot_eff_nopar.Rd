% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_effects_nopar.R
\name{plot_eff_nopar}
\alias{plot_eff_nopar}
\title{Plot direct, indirect and total effects functions 
  for continous non-parametric covariates in PS-SAR regression models.}
\usage{
plot_eff_nopar(effnopar, data, smooth = TRUE, span = c(0.1, 0.1, 0.2))
}
\arguments{
\item{effnopar}{object returned from \code{\link{eff_nopar}} function.}

\item{data}{dataframe with the data.}

\item{smooth}{logical value to choose smoothing of the effects function
prior to plot. Default TRUE.}

\item{span}{span for the kernel of the smoothing (see \code{\link{loess}} 
for details). Default c(0.1,0.1,0.2).}
}
\value{
plot of the direct, indirect and total effects function for each non-parametric
  covariate included in the object returned from \code{\link{effects_nopar}}.
}
\description{
Plot direct, indirect and total effect functions for 
  non-parametric covariates included in a semiparametric spatial
  or spatio-temporal SAR model. This model must include a spatial
  lag of the dependent variable (SAR) to have indirect effects 
  different from 0, otherwise, total and direct function effects 
  are the same. The effect functions can be smoothed to overcome 
  the instabilities created by the premultiplication of matrix
  \eqn{(I - \rho W)^{-1}}
}
\examples{
################################################
 ###################### Examples using a panel data of rate of
 ###################### unemployment for 103 Italian provinces in period 1996-2014.
library(sptpsar)
data(unemp_it); Wsp <- Wsp_it

######################  No Spatial Trend: PSAR including a spatial 
######################  lag of the dependent variable
form1 <- unrate ~ partrate + agri + cons +
                 pspl(serv,nknots=15) +
                 pspl(empgrowth,nknots=20) 
 gamsar <- psar(form1,data=unemp_it,sar=TRUE,Wsp=Wsp_it)
 summary(gamsar)
 ###### Non-Parametric Total, Direct and Indirect Effects
 list_varnopar <- c("serv","empgrowth")
 eff_nparvar <- eff_nopar(gamsar,list_varnopar)
 plot_effects_nopar(eff_nparvar,unemp_it,smooth=TRUE)
 plot_effects_nopar(eff_nparvar,unemp_it,smooth=FALSE)
 
 #' ######################   PSAR-ANOVA with spatial trend
form2 <- unrate ~ partrate + agri + cons +
                  pspl(serv,nknots=15) + pspl(empgrowth,nknots=20) +
                  pspt(long,lat,nknots=c(20,20),psanova=TRUE,
                  nest_sp1=c(1,2),nest_sp2=c(1,2))
##### Spatial trend fixed for period 1996-2014
geospanova_sar <- psar(form2,data=unemp_it,Wsp=Wsp_it,sar=TRUE,
                   control=list(thr=1e-1,maxit=200,trace=FALSE))
summary(geospanova_sar)
 ###### Non-Parametric Total, Direct and Indirect Effects
 list_varnopar <- c("serv","empgrowth")
 eff_nparvar <- eff_nopar(geospanova_sar,list_varnopar)
 plot_effects_nopar(eff_nparvar,unemp_it,smooth=TRUE)
 plot_effects_nopar(eff_nparvar,unemp_it,smooth=FALSE)
 
######################   PSAR-ANOVA with spatio-temporal trend and 
######################   temporal autorregresive noise
 form3 <- unrate ~ partrate + agri + cons +
                   pspl(serv,nknots=15) + pspl(empgrowth,nknots=20) +
                   pspt(long,lat,year,nknots=c(18,18,8),psanova=TRUE,
                   nest_sp1=c(1,2,3),nest_sp2=c(1,2,3),
                   nest_time=c(1,2,2),ntime=19)
sptanova_sar_ar1 <- psar(form3,data=unemp_it,Wsp=Wsp_it,sar=TRUE,ar1=TRUE,
                    control=list(thr=1e-1,maxit=200,trace=FALSE))
summary(sptanova_sar_ar1)
 ###### Non-Parametric Total, Direct and Indirect Effects
 list_varnopar <- c("serv","empgrowth")
 eff_nparvar <- eff_nopar(sptanova_sar_ar1,list_varnopar)
 plot_effects_nopar(eff_nparvar,unemp_it,smooth=TRUE)
 plot_effects_nopar(eff_nparvar,unemp_it,smooth=FALSE)

}
\references{
\itemize{ 
  \item Basile, R., Durbán, M., Mínguez, R., Montero, J.
        M., and Mur, J. (2014). Modeling regional economic 
        dynamics: Spatial dependence, spatial heterogeneity and 
        nonlinearities. \emph{Journal of Economic Dynamics and 
        Control}, (48), 229-245.                  
 }
}
\seealso{
\itemize{
  \item \code{\link{eff_nopar}} compute total, direct and indirect effect
          functions for non-parametric continuous covariates.
  \item \code{\link{fit_terms}} compute smooth functions for non-parametric
                                continuous covariates.
  \item \code{\link{plot_terms}} plot the terms of non-parametric covariates.
}

Other Direct, Indirect and Total Effects.: \code{\link{eff_nopar}},
  \code{\link{eff_par}}
}
\author{
Roman Minguez \email{roman.minguez@uclm.es}
}
\concept{Direct, Indirect and Total Effects.}
\keyword{@export}
\keyword{Direct}
\keyword{Indirect}
\keyword{SAR,}
\keyword{covariates.}
\keyword{effects,}
\keyword{non-parametric}
