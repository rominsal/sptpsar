#' @name logLik.psar
#' @rdname logLik.psar
#'
#' @title Extract method log-likelihood or restricted log-likelihood
#'   from objects of class psar.
#'
#' @param object object of class \emph{psar}.
#' @param REML logical value indicating if the log-likelihood is
#'   restricted. Default: FALSE.
#' @param ... further arguments passed to or from other methods.
#'
#' @return An object of class \emph{logLik}. Apart from the
#'   log-likelihood (or restricted log-likelihood) value, this object
#'   has some additional attributes:
#' \tabular{ll}{
#'   \code{nall} \tab Whole number of observations. \cr
#'   \code{nobs} \tab Effective number of observations. \cr
#'   \code{df} \tab Effective degrees of freedom. \cr
#'  }
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{coef.psar}} Extract vector of coefficients for
#'     fitted \emph{psar} objects including fixed and random effects.
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
#' llik_gamsar <- logLik(gamsar)
#' llikreml_gamsar <- logLik(gamsar, REML = TRUE)
#' cat("Log-Likelihood: ",llik_gamsar,"\n")
#' cat("Restricted Log-Likelihood: ",llikreml_gamsar,"\n")
#' @export
logLik.psar <- function(object, REML = FALSE, ...)
{
  res <- object$residuals
  edftot <- object$edftot
  N <- length(res)
  N0 <- N
  if (REML) {
    N <- N - edftot
    val <- object$llik_reml
  } else val <- object$llik
  attr(val, "nall") <- N0
  attr(val, "nobs") <- N
  attr(val, "df") <- edftot
  class(val) <- "logLik"
  val
}
