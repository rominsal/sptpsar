logLik.psar <- function(object, REML = FALSE, ...)
{
  res <- object$residuals
  edftot <- object$edftot
  N <- length(res)
  # Uncomment when weights are allowing...
  # if (is.null(w <- object$weights)) {
  #   w <- rep.int(1, N)
  # }
  # else {
  #   excl <- w == 0
  #   if (any(excl)) {
  #     res <- res[!excl]
  #     N <- length(res)
  #     w <- w[!excl]
  #   }
  # }
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
