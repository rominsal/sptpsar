#' Compute Total, Direct and Indirect effects functions for  geoadditive
#' spatial or spatio-temporal semiparametric PS-SAR regression models.
#'
#' Compute effects function for non-parametric covariates in semiparametric
#' models.
#'
#' @param fitted_nopar A matrix including the non-parametric fitted
#'    functions (GAM) for each covariate.
#' @param sd_fitted_nopar A matrix including the standard deviations of
#'    non-parametric functions(GAM) for each convariate.
#' @param Wsp A neighbour spatial matrix.
#'    The dimension of the matrix is always \emph{nxn}, where \emph{n} is the
#'    dimension of each spatial coordinate. Default NULL.
#' @param rho Spatial parameter of the spatial lag of the dependent
#'    variable for SAR models. Default 0.
#' @param durbin A logical value indicating if the model include a spatial lag
#'    of each covariate (either parametric or non-parametric). Default FALSE.
#' @param conf_level Numerical value for the confidence interval of the
#'    effect functions.
#'
eff_nopar <- function(fitted_nopar,sd_fitted_nopar,rho=0,
                      Wsp,conf_level=0.95,durbin=FALSE)
{
    #rho <- spt_model$rho
    #Xnopar <- spt_model$Xnopar
    #fitted_nopar <- spt_model$fitted_nopar
    #sd_fitted_nopar <- spt_model$sd_fitted_nopar
    nspt <- nrow(fitted_nopar)
    I_n <- Matrix::Diagonal(nspt)
    #I_rhoWinv <- solve(I-rho*W)
    if(!durbin) { p <- ncol(fitted_nopar) } else { p <- ncol(fitted_nopar)/2 }
    ADE <- matrix(NA,nrow=n,ncol=p)
    if(!is.null(colnames(fitted_nopar))) colnames(ADE) <- colnames(fitted_nopar)[1:p]
    ATE <- ADE; AIE <- ADE
    upper_ADE <- ADE; upper_ATE <- ATE; upper_AIE <- AIE
    lower_ADE <- ADE; lower_ATE <- ATE; lower_AIE <- AIE
    crit.val.norm01 <- qnorm((1-conf_level)/2,mean=0,sd=1,lower.tail=FALSE)
    #for(i in 1:nrow(fitted_nopar))
    #{
        for(j in 1:p)
        {
            fitted.vector <- fitted_nopar[,j]
            if(durbin) fitted.vector <- fitted.vector + fitted_nopar[,(p+j)]
            ADE[,j] <- diag(solve(I_n-rho*Wsp,diag(fitted.vector)))
            ATE[,j] <- solve(I_n-rho*Wsp,fitted.vector)
            AIE[,j] <- ATE[,j] - ADE[,j]

            upper_fitted.vector <- fitted_nopar[,j] + crit.val.norm01*sd_fitted_nopar[,j]
            if(durbin) upper_fitted.vector <- upper_fitted.vector +
                    (fitted_nopar[,(p+j)] + crit.val.norm01*sd_fitted_nopar[,(p+j)])
            upper_ADE[,j] <- diag(solve(I_n-rho*Wsp,diag(upper_fitted.vector)))
            upper_ATE[,j] <- solve(I_n-rho*Wsp,upper_fitted.vector)
            upper_AIE[,j] <- upper_ATE[,j] - upper_ADE[,j]

            lower_fitted.vector <- fitted_nopar[,j] - crit.val.norm01*sd_fitted_nopar[,j]
            if(durbin) lower_fitted.vector <- lower_fitted.vector +
                (fitted_nopar[,(p+j)] - crit.val.norm01*sd_fitted_nopar[,(p+j)])
            lower_ADE[,j] <- diag(solve(I_n-rho*Wsp,diag(lower_fitted.vector)))
            lower_ATE[,j] <- solve(I_n-rho*Wsp,lower_fitted.vector)
            lower_AIE[,j] <- lower_ATE[,j] - lower_ADE[,j]

            #S_W <- I_rhoWinv %*% (I*fitted_nopar[i,j])
            #S_W <- solve(I_n-rho*W, I_n*fitted_nopar[i,j])
            #ADE[,j] <- sum(diag(S_W))/n
            #ATE[i,j] <- sum(S_W)/n
            #AIE[i,j] <- ATE[i,j]-ADE[i,j]
            #upper_S_W <- solve(I_n-rho*W,
            #             I_n*(fitted_nopar[i,j]+crit.val.norm01*sd_fitted_nopar[i,j]))
            #upper_S_W <- I_rhoWinv %*%
            #    (I*(fitted_nopar[i,j]+crit.val.norm01*sd_fitted_nopar[i,j]))
            #upper_ADE[i,j] <- sum(diag(upper_S_W))/n
            #upper_ATE[i,j] <- sum(upper_S_W)/n
            #upper_AIE[i,j] <- upper_ATE[i,j]-upper_ADE[i,j]
            #lower_S_W <- I_rhoWinv %*%
            #    (I*(fitted_nopar[i,j]-crit.val.norm01*sd_fitted_nopar[i,j]))
            #lower_S_W <- solve(I_n-rho*W,
            #                   I_n*(fitted_nopar[i,j]+crit.val.norm01*sd_fitted_nopar[i,j]))
            #lower_ADE[i,j] <- sum(diag(lower_S_W))/n
            #lower_ATE[i,j] <- sum(lower_S_W)/n
            #lower_AIE[i,j] <- lower_ATE[i,j]-lower_ADE[i,j]
        }
    #}
    res <- list(ATE = ATE, ADE = ADE, AIE = AIE,
                upper_ATE = upper_ATE, upper_ADE = upper_ADE, upper_AIE = upper_AIE,
                lower_ATE = lower_ATE, lower_ADE = lower_ADE, lower_AIE = lower_AIE)
    res
} # end of function eff.npar



# for(i in 1:ncol(mat.gam.nlin))
# {
#   x <- mat.gam.nlin[,i]
#   name <- colnames(mat.gam.nlin)[i]
#   ord <- order(x)
#   maximoy <- max(eff$upper_ADE[,i]-mean(eff$ADE[,i]))
#   minimoy <- min(eff$lower_ADE[,i]-mean(eff$ADE[,i]))
#   #name_file <- paste("spt_anova_sar_nlin_ADE_AIE_",name,".png",sep="")
#   #png(name_file)
#   ADE <- (eff$ADE[,i]-mean(eff$ADE[,i]))
#   ADE.smooth <- predict(loess(ADE[ord]~x[ord], span = 0.1),method = "loess()")
#   upper.ADE <- eff$upper_ADE[,i][ord]-mean(eff$ADE[,i])
#   upper.ADE.smooth <- predict(loess(upper.ADE~x[ord], span = 0.1),method = "loess()")
#   lower.ADE <- eff$lower_ADE[,i][ord]-mean(eff$ADE[,i])
#   lower.ADE.smooth <- predict(loess(lower.ADE~x[ord], span = 0.1),method = "loess()")
#
#   AIE <- (eff$AIE[,i]-mean(eff$AIE[,i]))
#   AIE.smooth <- predict(loess(AIE[ord]~x[ord], span = 0.1),method = "loess()")
#   upper.AIE <- eff$upper_AIE[,i][ord]-mean(eff$AIE[,i])
#   upper.AIE.smooth <- predict(loess(upper.AIE~x[ord], span = 0.1),method = "loess()")
#   lower.AIE <- eff$lower_AIE[,i][ord]-mean(eff$AIE[,i])
#   lower.AIE.smooth <- predict(loess(lower.AIE~x[ord], span = 0.1),method = "loess()")
#
#   ATE <- (eff$ATE[,i]-mean(eff$ATE[,i]))
#   ATE.smooth <- predict(loess(ATE[ord]~x[ord], span = 0.1),method = "loess()")
#   upper.ATE <- eff$upper_ATE[,i][ord]-mean(eff$ATE[,i])
#   upper.ATE.smooth <- predict(loess(upper.ATE~x[ord], span = 0.1),method = "loess()")
#   lower.ATE <- eff$lower_ATE[,i][ord]-mean(eff$ATE[,i])
#   lower.ATE.smooth <- predict(loess(lower.ATE~x[ord], span = 0.1),method = "loess()")
#
#   plot(x[ord],ADE.smooth,type="l",
#        #y[ord],xlab=name,
#        ylab="",xlab="Pollution Level",
#        #ylab=paste("ADE(",name,")",sep=""),
#        ylim=c(minimoy,maximoy),cex.lab=1.5,,col=2,lty=1)
#        #cex.main=1.5,main="Direct Effects")
#   lines(x[ord],upper.ADE.smooth,xlab="",ylab="",type="l",col=2,lty=1)
#   lines(x[ord],lower.ADE.smooth,xlab="",ylab="",type="l",col=2,lty=1)
#   lines(x[ord],AIE.smooth,xlab="",ylab="",type="l",col=4,lty=3)
#   lines(x[ord],upper.AIE.smooth,xlab="",ylab="",type="l",col=4,lty=3)
#   lines(x[ord],lower.AIE.smooth,xlab="",ylab="",type="l",col=4,lty=3)
#   abline(a=0,b=0)
#
#   # PLOT OF TOTAL EFFECTS
#   maximoy <- max(eff$upper_ATE[,i]-mean(eff$ATE[,i]))
#   minimoy <- min(eff$lower_ATE[,i]-mean(eff$ATE[,i]))
#
#   plot(x[ord],ATE.smooth,type="l",
#        #y[ord],xlab=name,
#        ylab="",xlab="Pollution Level",
#        #ylab=paste("ADE(",name,")",sep=""),
#        ylim=c(minimoy,maximoy),cex.lab=1.5,,col=2,lty=1)
#   #cex.main=1.5,main="Direct Effects")
#   lines(x[ord],upper.ATE.smooth,xlab="",ylab="",type="l",col=2,lty=1)
#   lines(x[ord],lower.ATE.smooth,xlab="",ylab="",type="l",col=2,lty=1)
#
#   dev.off()
# }

#
# for(i in 1:ncol(Xnopar))
# {
#   x <- Xnopar[,i]
#   y <- ATE[,i]
#   name <- colnames(Xnopar)[i]
#   ord <- order(x)
#   maximoy <- max(upper_ATE[,i])
#   minimoy <- min(lower_ATE[,i])
#   name_file <- paste("spt_anova_sar_nlin_ATE_",name,".pdf",sep="")
#   pdf(name_file)
#   plot(x[ord],y[ord],xlab=name,
#        ylab=paste("ATE(",name,")",sep=""),type="l",
#        ylim=c(minimoy,maximoy),cex.lab=1.5,
#        cex.main=1.5,main="Total Effects")
#   lines(x[ord],upper_ATE[,i][ord],xlab="",ylab="",col=1,lty=4)
#   lines(x[ord],lower_ATE[,i][ord],xlab="",ylab="",col=1,lty=4)
#   dev.copy(which=4)
# }
#
# for(i in 1:ncol(Xnopar))
# {
#   x <- Xnopar[,i]
#   y <- AIE[,i]
#   name <- colnames(Xnopar)[i]
#   ord <- order(x)
#   maximoy <- max(upper_AIE[,i])
#   minimoy <- min(lower_AIE[,i])
#   name_file <- paste("spt_anova_sar_nlin_AIE_",name,".pdf",sep="")
#   pdf(name_file)
#   plot(x[ord],y[ord],xlab=name,
#        ylab=paste("AIE(",name,")",sep=""),type="l",
#        ylim=c(minimoy,maximoy),cex.lab=1.5,
#        cex.main=1.5,main="Indirect Effects")
#   lines(x[ord],upper_AIE[,i][ord],xlab="",ylab="",col=1,lty=4)
#   lines(x[ord],lower_AIE[,i][ord],xlab="",ylab="",col=1,lty=4)
#   dev.copy(which=4)
# }

impacts_par <- function(n,beta_par,cov.beta_par,rho,sd_rho,durbin=FALSE,
                            traceW=NULL,W=NULL,nrep=1000,seed=1111,tol=0.1)
    # Function to compute marginal effects in Spatio-Temporal Models
    # for parametric linear covariates.

# REPASAR: SALE SIEMPRE EL MISMO ESTADÍSTICO T PARA TODOS LOS EFECTOS
# ADE,ATE Y AIE DE CADA VARIABLE (AUNQUE LA MEDIA Y LA DESV. TÍPICA SON DISTINTAS)
{
    if(!is.null(seed)) set.seed(seed)
    if(is.null(traceW) && !is.null(W)) spdep::trW(W,m=100,p=50, type="MC")
    rho_sim <- rnorm(nrep,rho,sd_rho)
    beta_par.sim <- MASS::mvrnorm(nrep,beta_par,cov.beta_par,tol=tol)
    beta_par.sim <- t(beta_par.sim)
    rownames(beta_par.sim) <- names(beta_par)
    if(!durbin) { p <- length(beta_par) } else { p <- length(beta_par)/2 }
    # Code based on the definition (I-rho*W)^(-1)*(I_n*\beta_i+W*\theta_i)
    ADE <- matrix(NA,nrow=p,ncol=nrep)
    ATE <- matrix(NA,nrow=p,ncol=nrep)
    # for(i in 1:p)
    # {
    #     for(j in 1:nrep)
    #     {
    #         if(!durbin)
    #         {
    #             imp.matrix.ij <- Matrix::solve(Diagonal(nspt)-rho*Wdist.mat.tot,
    #                                     Diagonal(nspt)*beta_par.sim[i,j])
    #         } else {
    #             imp.matrix.ij <- Matrix::solve(Diagonal(nspt)-rho*Wdist.mat.tot,
    #                                     Diagonal(nspt)*beta_par.sim[i,j] +
    #                                     W*beta_par.sim[p+i,j])
    #         }
    #         ADE[i,j] <- sum(diag(as.matrix(imp.matrix.ij))) / n
    #         ATE[i,j] <- sum(as.matrix(imp.matrix.ij)) / n
    #     }
    # }
    #
    # # Code based on  LeSage and Page (2009) pp. 114-115
    if(!durbin)
    {
        q <- length(traceW)
    #    p <- length(beta_par)
    } else {
        q <- length(traceW) - 1
    #    p <- length(beta_par)/2
    }
    #I_n <- diag(n)
    #ADE[i,j] <- sum(diag(solve(I_n-rho*W,I_n*beta_par.sim[i,j])))
    #ATE[i,j] <- solve(I_n-rho*W,fitted.vector)

    mT <- matrix(c(1,traceW[1:q]/n),nrow=1)
    if(durbin) mT <- rbind(mT,traceW[1:(q+1)]/n)
#    g <- 1
#    for(i in 1:q) { g <- c(g,rho^i) }
    a <- matrix(1,nrow=q+1,ncol=1)
#   G <- diag(g)
    mP <- list()
    for (i in 1:nrep)
    {
        mP[[i]] <- matrix(beta_par.sim[1:p,i],ncol=1)
        if(durbin) mP[[i]] <- cbind(mP[[i]],beta_par.sim[(p+1):(2*p),i])
        g <- 1
        for(j in 1:q) { g <- c(g,rho_sim[i]^j) }
        G <- diag(g)
        ADE[,i] <- mP[[i]] %*% mT %*% G %*% a
        ATE[,i] <- matrix(rowSums(mP[[i]]),ncol=1) %*% (matrix(g,nrow=1) %*% a)
    }
    AIE <- ATE - ADE
    rownames(ADE) <- names(beta_par)[1:p]
    rownames(ATE) <- names(beta_par)[1:p]
    rownames(AIE) <- names(beta_par)[1:p]
    ADE.mean <- apply(ADE,1,mean)
    ATE.mean <- apply(ATE,1,mean)
    AIE.mean <- apply(AIE,1,mean)
    ADE.sd <- apply(ADE,1,sd)
    ATE.sd <- apply(ATE,1,sd)
    AIE.sd <- apply(AIE,1,sd)

    ADE.tstat <- ADE.mean / ADE.sd
    ATE.tstat <- ATE.mean / ATE.sd
    AIE.tstat <- AIE.mean / AIE.sd

    ADE.pvalue <- 2*pnorm(abs(ADE.tstat),mean=0,sd=1,lower.tail=FALSE)
    ATE.pvalue <- 2*pnorm(abs(ATE.tstat),mean=0,sd=1,lower.tail=FALSE)
    AIE.pvalue <- 2*pnorm(abs(AIE.tstat),mean=0,sd=1,lower.tail=FALSE)

    imp.ADE <- cbind(ADE.mean,ADE.sd,ADE.tstat,ADE.pvalue)
    imp.ATE <- cbind(ATE.mean,ATE.sd,ATE.tstat,ATE.pvalue)
    imp.AIE <- cbind(AIE.mean,AIE.sd,AIE.tstat,AIE.pvalue)
    colnames(imp.ADE) <- c("Dir. Imp.","Std. Dev.","t-stat","p-value")
    colnames(imp.ATE) <- c("Tot. Imp.","Std. Dev.","t-stat","p-value")
    colnames(imp.AIE) <- c("Ind. Imp.","Std. Dev.","t-stat","p-value")

    res <- list(ADE.sim=ADE,ATE.sim=ATE,AIE.sim=AIE,
                table.imp.ADE=imp.ADE,table.imp.ATE=imp.ATE,table.imp.AIE=imp.AIE)
    res
}
