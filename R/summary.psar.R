#' @name summary.psar
#' @rdname summary.psar
#'
#' @title Summary method for objects of class psar.
#'
#' @description This method summarizes both spatial (2-dimension) and
#'   spatio-temporal (3-dimension) \emph{psar} objects.
#'   The tables include information of:
#'   \itemize{
#'      \item The spatial (or spatio-temporal) trends. When the model is ANOVA
#'        the trend is decomposed in main and interaction effects.
#'      \item The parametric and non-parametric covariates.
#'      \item The \eqn{\rho} parameter when the model is SAR.
#'      \item The \eqn{\phi} parameter when the model is spatio-temporal
#'        with a first-order autorregressive in the noise.
#'  }
#'
#' @param object \emph{psar} object fitted using \code{\link{psar}} function.
#' @param ... further arguments passed to or from other methods.
#'
#' @return An object of class \emph{summary.psar}
#'
#' @author Roman Minguez \email{roman.minguez@@uclm.es}
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{psar}} estimate spatial or spatio-temporal semiparametric PS-SAR
#'   regression models.
#'   \item \code{\link{print.summary.psar}} print objects of class \emph{summary.psar}
#' }
#'
#' @examples
#'  See examples for \code{\link{psar}} function.
#'
#' @export
summary.psar <- function(object,...)
{
 z <- object

 names_var <- labels(z$terms)
 names_varspt <- names_var[grepl("pspt",names_var)]
 nvarspt <- length(names_varspt)
 names_varnopar <- names_var[grepl("pspl",names_var)]
 nvarnopar <- length(names_varnopar)
 names_varpar <- names(z$bfixed) # To include intercept (if there is...)
 names_varpar <- names_varpar[!grepl("pspl",names_varpar) & !grepl("pspt",names_varpar)]
 nvarpar <- length(names_varpar)

 rdf <- z$df.residual
 r <- z$residuals
 f <- z$fitted.values
 rss <- sum(r^2)
 resvar <- rss/rdf
 names_varpar <- gsub("fixed_","",names_varpar)
 match_names_varpar <- unique(grep(paste(c("spt","_main","_int"),
                                   collapse="|"),
                             names_varpar, value=TRUE))
 names_varpar <- names_varpar[ !(names_varpar %in% match_names_varpar)]
 if ("Intercept" %in% names_varpar){
   names_varpar <- c("Intercept",names_varpar[!(names_varpar %in% "Intercept")])
 }
 names(z$bfixed) <- gsub("fixed_","",names(z$bfixed))
 names(z$se_bfixed) <- gsub("fixed_","",names(z$se_bfixed))
 est_par <- z$bfixed[names_varpar]
 se_par <- z$se_bfixed[names_varpar]
 tval_par <- est_par/se_par
 z$coef_par_table <- cbind(est_par, se_par, tval_par,
                      2 * pt(abs(tval_par),rdf, lower.tail = FALSE))
 colnames(z$coef_par_table) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
 rownames(z$coef_par_table) <- names_varpar
 if(z$sar)
 {
   t_rho <- z$rho / z$se_rho
   pval_rho <- 2 * pt(abs(t_rho),rdf, lower.tail = FALSE)
   z$coef_par_table <- rbind(z$coef_par_table,
                         c(z$rho,z$se_rho,t_rho,pval_rho))
   n <- length(rownames(z$coef_par_table))
   rownames(z$coef_par_table)[n] <- c("rho")
 }
 if (!is.null(z$time)){
   if(z$ar1)
   {
     t_phi <- z$phi / z$se_phi
     pval_phi <- 2 * pt(abs(t_phi),rdf, lower.tail = FALSE)
     z$coef_par_table <- rbind(z$coef_par_table,
                               c(z$phi,z$se_phi,t_phi,pval_phi))
     n <- length(rownames(z$coef_par_table))
     rownames(z$coef_par_table)[n] <- c("phi")
   }
 }

 if (!is.null(z$edfspt))
 {
   if(z$psanova) { # is.null(time)
     if(is.null(z$time)){
       z$coef_spttrend_table <- matrix(0,nrow=3,ncol=2)
       z$coef_spttrend_table[1,1] <- z$edfspt[c("f1_main")]
       z$coef_spttrend_table[2,2] <- z$edfspt[c("f2_main")]
       z$coef_spttrend_table[3,c(1,2)] <- z$edfspt[c("f12.1","f12.2")]
       colnames(z$coef_spttrend_table) <- c("EDF sp1","EDF sp2")
       rownames(z$coef_spttrend_table) <- c("f1","f2","f12")
     } else {# !is.null(time)
       z$coef_spttrend_table <- matrix(0,nrow=7,ncol=3)
       z$coef_spttrend_table[1,1] <- z$edfspt[c("f1_main")]
       z$coef_spttrend_table[2,2] <- z$edfspt[c("f2_main")]
       z$coef_spttrend_table[3,3] <- z$edfspt[c("ft_main")]
       z$coef_spttrend_table[4,c(1,2)] <- z$edfspt[c("f12.1","f12.2")]
       z$coef_spttrend_table[5,c(1,3)] <- z$edfspt[c("f1t.1","f1t.2")]
       z$coef_spttrend_table[6,c(2,3)] <- z$edfspt[c("f2t.1","f2t.2")]
       z$coef_spttrend_table[7,c(1,2,3)] <- z$edfspt[c("f12t.1","f12t.2","f12t.3")]
       colnames(z$coef_spttrend_table) <- c("EDF sp1","EDF sp2","EDF time")
       rownames(z$coef_spttrend_table) <- c("f1","f2","ft","f12","f1t","f2t","f12t")
     }

   } else { # PSANOVA = FALSE
     if (is.null(z$time)){
       z$coef_spttrend_table <- matrix(0,nrow=2,ncol=1)
       z$coef_spttrend_table[1,1] <- z$edfspt[c("sp1")]
       z$coef_spttrend_table[2,1] <- z$edfspt[c("sp2")]
       colnames(z$coef_spttrend_table) <- c("EDF")
       rownames(z$coef_spttrend_table) <- c("sp1","sp2")
     } else {
       z$coef_spttrend_table <- matrix(0,nrow=3,ncol=1)
       z$coef_spttrend_table[1,1] <- z$edfspt[c("sp1")]
       z$coef_spttrend_table[2,1] <- z$edfspt[c("sp2")]
       z$coef_spttrend_table[3,1] <- z$edfspt[c("time")]
       colnames(z$coef_spttrend_table) <- c("EDF")
       rownames(z$coef_spttrend_table) <- c("sp1","sp2","time")
     }
   }
 } else { z$coef_spttrend_table <- NULL}

 if( !is.null(z$coef_spttrend_table)){
   if (any(is.na(z$coef_spttrend_table))) {
     index_na_table <- which(is.na(z$coef_spttrend_table))
     z$coef_spttrend_table[c(index_na_table)] <- 0
   }
 }

 if(!is.null(z$edfnopar))
 {
   z$coef_nopar_table <- matrix(z$edfnopar,ncol=1)
   colnames(z$coef_nopar_table) <- c("EDF")
   rownames(z$coef_nopar_table) <- names_varnopar
 }
 z$sigma <- sqrt(resvar)
 z$edftot <- z$edftot
 class(z) <- c("summary.psar",class(z))
 z
}
