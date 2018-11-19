#' sptpsar
#'
#' Estimation and inference of spatial and spatio-temporal
#' semiparametric models including spatial or spatio-temporal
#' non-parametric trends, parametric and non-parametric covariates 
#' and, possibly, a spatial lag for the dependent variable and 
#' temporal correlation in the noise. 
#' The spatio-temporal trend can be decomposed
#' in ANOVA way including main and interaction functional terms.
#' Use of SAP algorithm to estimate the spatial or spatio-temporal
#' trend and non-parametric covariates.
#'
#' @docType package
#' @name sptpsar
#' @author Román Mínguez

#' @importFrom graphics abline  lines matplot par plot
#' @importFrom MASS ginv mvrnorm
#' @importFrom Matrix bdiag crossprod Diagonal Matrix solve t
#' @importFrom methods as
#' @importFrom nloptr nloptr
#' @importFrom numDeriv hessian
#' @importFrom plyr alply
#' @importFrom rootSolve multiroot
#' @importFrom spdep trW lagsarlm mat2listw
#' @importFrom splines spline.des
#' @importFrom stats AIC BIC logLik loess lm model.response predict terms
#' @importFrom stats pnorm pt qnorm rnorm sd var update 
#' @importFrom stats .getXlevels printCoefmat weights
#' @useDynLib sptpsar
#' @exportPattern "^[[:alpha:]]+"
NULL
