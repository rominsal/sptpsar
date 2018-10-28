#' sptpsar
#' 
#' Estimation and inference of spatial and spatio-temporal 
#' semiparametric models including spatial or spatio-temporal 
#' non-parametric trends, parametric and non-parametric covariates and, 
#' possibly, a spatial lag for the dependent variable and temporal correlation 
#' in the noise. The spatio-temporal trend can be decomposed  
#' in ANOVA way including main and interaction functional terms.Use 
#' of SAP algorithm to estimate the spatial or spatio-temporal 
#' trend and non-parametric covariates. 
#' 
#' @docType package
#' @name sptpsar
#' @author Román Mínguez
#' @import Rcpp RcppEigen
#' @importFrom Rcpp evalCpp sourceCpp 
#' @useDynLib sptpsar
#' @exportPattern "^[[:alpha:]]+"
NULL  
