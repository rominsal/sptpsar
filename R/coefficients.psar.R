coef.psar <- function(object, ...)
{
  coeff <- c(object$bfixed,object$brandom)
  names_coeff <- names(coeff)
  if(object$sar) {
    coeff <- c(coeff,object$rho)
    names_coeff <- c(names_coeff,"rho")
  }
  if(object$ar1) {
    coeff <- c(coeff,object$phi)
    names_coeff <- c(names_coeff,"phi")
  }
  names(coeff) <- names_coeff
  coeff
}
