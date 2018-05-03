print.summary.par.eff.psar <- function(z, digits = max(3L, getOption("digits") - 3L), ...)
{

  if(!is.null(z$tot_table)) {
    cat("\n Total Parametric Effects \n")
    printCoefmat( z$tot_table, P.value=FALSE, has.Pvalue=FALSE)
  }
  if(!is.null(z$dir_table)) {
    cat("\n Direct Parametric Effects \n")
    printCoefmat( z$dir_table, P.value=FALSE, has.Pvalue=FALSE)
  }
  if(!is.null(z$ind_table)) {
    cat("\n Indirect Parametric Effects \n")
    printCoefmat( z$ind_table, P.value=FALSE, has.Pvalue=FALSE)
  }
  invisible(z)
}
