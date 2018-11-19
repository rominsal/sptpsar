#' @name coef.psar
#' @rdname coef.psar
#'
#' @title Extract method for model coefficients from objects of class psar.
#'
#' @param z object of class \emph{psar}.
#' @param ... further arguments passed to or from other methods.
#'
#' @return A vector of the coefficients of the model including fixed
#'   and random effects.
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{logLik.psar}} Extract log-likelihood or restricted
#'     log-likelihood for fitted \emph{psar} objects.
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
#' gamsar <- psar(form1,data=unemp_it,sar=TRUE,Wsp=Wsp_it)
#' summary(gamsar)
#' coef(gamsar)
#' @export

coef.psar <- function(object, ...)
{
  coeff <- c(object$bfixed,object$brandom)
  names_coeff <- names(coeff)
  if(object$sar) {
    coeff <- c(coeff,object$rho)
    names_coeff <- c(names_coeff,"rho")
  }
  if(!is.null(object$ar1) && object$ar1) {
    coeff <- c(coeff,object$phi)
    names_coeff <- c(names_coeff,"phi")
  }
  names(coeff) <- names_coeff
  coeff
}
