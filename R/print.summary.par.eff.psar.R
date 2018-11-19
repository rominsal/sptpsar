#' @name print.summary.psar
#' @rdname print.summary.psar
#'
#' @title Print method for objects of class summary.par.eff.psar.
#'
#' @param x object of class \emph{summary.par.eff.psar}.
#' @param digits number of digits to show in printed tables.
#'   Default: max(3L, getOption("digits") - 3L).
#' @param ... further arguments passed to or from other methods.
#'
#' @author Roman Minguez \email{roman.minguez@@uclm.es}
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{eff_par}} Compute direct, indirect and
#'     total effects (or impacts).
#'     for continous parametric covariates in PS-SAR regression models.
#'   \item \code{\link{summary.par.eff.psar}} Summary method
#'     for \emph{par.eff.psar} objects.
#' }
#'
#' @examples
#'   See examples for \code{\link{effects_par}} function.
#' @export
print.summary.par.eff.psar <- function(x, 
      digits = max(3L, getOption("digits") - 3L), ...)
{
  if(!is.null(x$tot_table)) {
    cat("\n Total Parametric Effects \n")
    printCoefmat( x$tot_table, P.values=FALSE, has.Pvalue=FALSE)
  }
  if(!is.null(x$dir_table)) {
    cat("\n Direct Parametric Effects \n")
    printCoefmat( x$dir_table, P.values=FALSE, has.Pvalue=FALSE)
  }
  if(!is.null(x$ind_table)) {
    cat("\n Indirect Parametric Effects \n")
    printCoefmat( x$ind_table, P.values=FALSE, has.Pvalue=FALSE)
  }
  invisible(x)
}
