anova.psar <- function(object,..., dispersion = NULL, test = NULL, freq = FALSE)
{
  objects <- list(object, ...)
  ns <- sapply(objects, function(x) length(x$residuals))
  if (any(ns != ns[1L]))
    stop("models were not all fitted to the same size of dataset")
  nmodels <- length(objects)
  resdf <- resdev <- resrss <- ressigmasq <- numeric(nmodels)
  rlogLiks <- logLiks <- aics <- bics <- edfs <- numeric(nmodels)
  for (i in 1:nmodels){
    resdf[i] <- objects[[i]]$df.residual
    resrss[i] <- sum(objects[[i]]$residuals^2)
    ressigmasq[i] <- resrss[i] / resdf[i]
    llik_i <- logLik(objects[[i]], REML = FALSE)
    logLiks[i] <- as.numeric(llik_i)
    rllik_i <- logLik(objects[[i]], REML = TRUE)
    rlogLiks[i] <- as.numeric(rllik_i)
    resdev[i] <- -2*logLiks[i]
    aics[i] <- as.numeric(AIC(llik_i))
    bics[i] <- as.numeric(BIC(llik_i))
    edfs[i] <- attr(llik_i,"df")
    rm(llik_i,rllik_i)
  }
  table <- data.frame(edfs,logLiks,rlogLiks,aics,bics,
                      resdf,ressigmasq)
  dimnames(table) <- list(1L:nmodels, c("EDF","logLik",
                                        "RlogLik","AIC","BIC",
                                        "Res.Df", "Sigma.Sq"))
  print(table) # Change the format...
  class(table) <- "anova.psar"
  invisible(table)
}
