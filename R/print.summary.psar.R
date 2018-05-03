print.summary.psar <- function(z, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("Call:\n")
  print(z$call)
  if(!is.null(z$coef_par_table)) {
    cat("\n Parametric Terms \n")
    printCoefmat(z$coef_par_table, P.value=TRUE, has.Pvalue=TRUE)
  }
  if(!is.null(z$coef_nopar_table)) {
    cat("\n Non-Parametric Terms \n")
    printCoefmat( z$coef_nopar_table, P.value=FALSE, has.Pvalue=FALSE)
  }
  if(!is.null(z$coef_spttrend_table)) {
    cat("\n Non-Parametric Spatio-Temporal Trend \n")
    printCoefmat( round(z$coef_spttrend_table,3),
                  P.value=FALSE, has.Pvalue=FALSE)
  }
  cat("\n Goodness-of-Fit \n")
  cat("\nEDF Total:",formatC(z$edftot,digits=6,width=6),
      " Sigma:",formatC(z$sigma,digits=6,width=6))
  cat("\nAIC:      ",formatC(z$aic,digits=6,width=6),
      "BIC: ",formatC(z$bic,digits=6,width=6))

  invisible(z)
}
