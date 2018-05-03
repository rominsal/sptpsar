plot_terms <- function(fit_terms,data,conflevel=0.95){
  # Function to plot fitted terms X*\hat{\beta} + Z*\hat{\alpha}

  fit <- fit_terms$fitted_terms
  se_fit <- fit_terms$se_fitted_terms
  fit_fixed <- fit_terms$fitted_terms_fixed
  se_fit_fixed <- fit_terms$se_fitted_terms_fixed
  fit_random <- fit_terms$fitted_terms_random
  se_fit_random <- fit_terms$se_fitted_terms_random
  variables <- colnames(fit)
  crval <- qnorm((1-conflevel)/2,mean=0,sd=1,lower.tail=FALSE)

  for (i in 1:length(variables)) {
    name_var <- variables[i]

    fit_var <- matrix(fit[,c(name_var)],ncol=1)
    colnames(fit_var) <- name_var
    se_fit_var <- matrix(se_fit[,c(name_var)],ncol=1)
    colnames(se_fit_var) <- name_var
    up_fit_var <- fit_var + crval*se_fit_var
    colnames(up_fit_var) <- name_var
    low_fit_var <- fit_var - crval*se_fit_var
    colnames(low_fit_var) <- name_var

    fit_var_fixed <- matrix(fit_fixed[,c(name_var)],ncol=1)
    colnames(fit_var_fixed) <- name_var
    se_fit_var_fixed <- matrix(se_fit_fixed[,c(name_var)],ncol=1)
    colnames(se_fit_var_fixed) <- name_var
    up_fit_var_fixed <- fit_var_fixed + crval*se_fit_var_fixed
    colnames(up_fit_var_fixed) <- name_var
    low_fit_var_fixed <- fit_var_fixed - crval*se_fit_var_fixed
    colnames(low_fit_var_fixed) <- name_var

    fit_var_random <- matrix(fit_random[,c(name_var)],ncol=1)
    colnames(fit_var_random) <- name_var
    se_fit_var_random <- matrix(se_fit_random[,c(name_var)],ncol=1)
    colnames(se_fit_var_random) <- name_var
    up_fit_var_random <- fit_var_random + crval*se_fit_var_random
    colnames(up_fit_var_random) <- name_var
    low_fit_var_random <- fit_var_random - crval*se_fit_var_random
    colnames(low_fit_var_random) <- name_var

    var <- matrix(data[,c(name_var)],ncol=1)
    colnames(var) <- name_var
    ord <- order(var)
    par(mfrow=c(2,1))
    plot(var[ord],fit_var[ord],type="l",
         ylab=paste("f(",name_var,")"),xlab=name_var,
         ylim=c(min(low_fit_var),max(up_fit_var)),cex.lab=1.0,col=2,lty=1,lwd=2,
         cex.main=1.0,
         main=paste("Term: ",name_var),
         sub="Confidence Intervals in dashed lines")
    lines(var[ord],up_fit_var[ord],xlab="",ylab="",type="l",col=2,lty=2,lwd=1.5)
    lines(var[ord],low_fit_var[ord],xlab="",ylab="",type="l",col=2,lty=2,lwd=1.5)
    abline(a=0,b=0)


    plot(var[ord],fit_var[ord],type="l",
         ylab=paste("f(",name_var,")"),xlab=name_var,
         ylim=c(min(low_fit_var),max(up_fit_var)),cex.lab=1.0,col=2,lty=1,lwd=2,
         cex.main=1.0,
         main=paste("Global (red), Fixed (green) and Random (blue) terms"))
    #lines(var[ord],up_fit_var[ord],xlab="",ylab="",type="l",col=2,lty=2)
    #lines(var[ord],low_fit_var[ord],xlab="",ylab="",type="l",col=2,lty=2)
    lines(var[ord],fit_var_fixed[ord],xlab="",ylab="",type="l",col=3,lty=2,lwd=2)
    #lines(var[ord],up_fit_var_fixed[ord],xlab="",ylab="",type="l",col=3,lty=2)
    #lines(var[ord],low_fit_var_fixed[ord],xlab="",ylab="",type="l",col=3,lty=2)
    lines(var[ord],fit_var_random[ord],xlab="",ylab="",type="l",col=4,lty=3,lwd=2)
    #lines(var[ord],up_fit_var_random[ord],xlab="",ylab="",type="l",col=4,lty=2)
    #lines(var[ord],low_fit_var_random[ord],xlab="",ylab="",type="l",col=4,lty=2)
    abline(a=0,b=0)
    readline(prompt="Press [enter] to continue")
  }
}




