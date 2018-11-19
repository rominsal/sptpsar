#' @name anova.psar
#' @rdname anova.psar
#'
#' @title Build an ANOVA table of fitted semiparametric models.
#' @description Method to compare goodness-of-fit measures of
#'   fitted spatio-temporal semiparametric PS-SAR regression models.
#'
#' @param object a fitted semiparametric model of class \emph{psar}.
#' @param ... additional fitted models to be compared.
#' @return A class \emph{anova.psar} table including
#'   goodness-of-fit measures for each fitted model.
#'   This table is also printed.
#' @author Roman Minguez \email{roman.minguez@@uclm.es}
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
#' gam <- psar(form1,data=unemp_it)
#' gamsar <- psar(form1,data=unemp_it,sar=TRUE,Wsp=Wsp_it)
#' form2 <- unrate ~ partrate + agri + cons +
#'   pspl(serv,nknots=15) + pspl(empgrowth,nknots=20) +
#'   pspt(long,lat,nknots=c(20,20),psanova=TRUE)
#' gamsp2d <- psar(form2,data=unemp_it)
#' gamsp2dsar <- psar(form2,data=unemp_it,sar=TRUE,Wsp=Wsp_it)
#' anova(gam,gamsar,gamsp2d,gamsp2dsar)
#' @export
anova.psar <- function(object, ...)
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
  print(table)
  class(table) <- "anova.psar"
  invisible(table)
}
