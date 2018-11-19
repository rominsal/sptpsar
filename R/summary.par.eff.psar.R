#' @name summary.par.eff.psar
#' @rdname summary.par.eff.psar
#'
#' @title Summary method for object of class par.eff.psar.
#'
#' @description This method summarizes direct, indirect and total effects (or impacts)
#'   for continous parametric covariates in PS-SAR regression models.
#'
#' @param object \emph{par.eff.psar} object fitted using \code{\link{psar}} function.
#' @param ... further arguments passed to or from other methods.
#'
#' @return An object of class \emph{summary.par.eff.psar}
#'
#' @author Roman Minguez \email{roman.minguez@@uclm.es}
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{eff_par}} Compute direct, indirect and total
#'     effects (or impacts)
#'     for continous parametric covariates in PS-SAR regression models.
#'   \item \code{\link{print.summary.par.eff.psar}} print objects of
#'     class \emph{summary.par.eff.psar}
#' }
#'
#' @examples
#'   See examples for \code{\link{effects_par}} function.
#' @export
summary.par.eff.psar <- function(object,...)
{
 z <- object
 tot <- z$tot_eff
 dir <- z$dir_eff
 ind <- z$ind_eff
 varpar <- rownames(tot)
 nrep <- ncol(tot)

 mean_dir <- apply(dir,1,mean)
 mean_tot <- apply(tot,1,mean)
 mean_ind <- apply(ind,1,mean)
 sd_dir <- apply(dir,1,sd)
 sd_tot <- apply(tot,1,sd)
 sd_ind <- apply(ind,1,sd)
 t_dir <- mean_dir / sd_dir
 t_tot <- mean_tot / sd_tot
 t_ind <- mean_ind / sd_ind


 z$tot_table <- cbind(mean_tot, sd_tot, t_tot,
                      2 * pnorm(abs(t_tot),mean=0,sd=1,lower.tail = FALSE))
 colnames(z$tot_table) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
 rownames(z$tot_table) <- varpar

 z$dir_table <- cbind(mean_dir, sd_dir, t_dir,
                      2 * pnorm(abs(t_dir),mean=0,sd=1,lower.tail = FALSE))
 colnames(z$dir_table) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
 rownames(z$dir_table) <- varpar

 z$ind_table <- cbind(mean_ind, sd_ind, t_ind,
                      2 * pnorm(abs(t_ind),mean=0,sd=1,lower.tail = FALSE))
 colnames(z$ind_table) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
 rownames(z$ind_table) <- varpar
 class(z) <- c("summary.par.eff.psar",class(z))
 z
}
