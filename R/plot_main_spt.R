plot_main_spt <- function(spttrend,sp1,sp2,nT,time=NULL,conflevel=0.95){
# Function to plot main effects in spatio-temporal trends

  fit_spt <- spttrend$fitted_terms
  se_fit_spt <- spttrend$se_fitted_terms
  n <- nrow(fit_spt)
  if (is.null(time)){ # 2d
    main_effects <- c("f1_main","f2_main")
    par(mfrow=c(2,1))
  } else {
    main_effects <- c("f1_main","f2_main","ft_main")
    par(mfrow=c(3,1))
  }
  crval <- qnorm((1-conflevel)/2,mean=0,sd=1,lower.tail=FALSE)
  for (i in 1:length(main_effects)){
    main <- main_effects[i]
    if (main!="ft_main"){
      seq_sel <- seq(from=1,to=n,by=nT)
    } else {
      seq_sel <- seq(from=1,to=nT,by=1)
    }
    fmain <- fit_spt[seq_sel,main]
    se_fmain <- se_fit_spt[seq_sel,main]
    up_fmain <- fmain + crval*se_fmain
    low_fmain <- fmain - crval*se_fmain
    if (main=="f1_main") {
      coord <- sp1[seq_sel]
      names(coord) <- "sp1"
    } else if (main=="f2_main") {
      coord <- sp2[seq_sel]
      names(coord) <- "sp2"
    } else {
      coord <- time[seq_sel]
      names(coord) <- "time"
    }
    ord <- order(coord)
    plot(coord[ord],fmain[ord],type="l",
         ylab=main,xlab=names(coord),
         ylim=c(min(low_fmain),max(up_fmain)),
         cex.lab=1.0,col=2,lty=1,lwd=2,cex.main=1.0,
         main=paste("Spatio-Temporal ANOVA Effect: ",main))
    lines(coord[ord],up_fmain[ord],xlab="",ylab="",type="l",col=2,lty=2,lwd=1.5)
    lines(coord[ord],low_fmain[ord],xlab="",ylab="",type="l",col=2,lty=2,lwd=1.5)
    abline(a=0,b=0)
  }
  par(mfrow=c(1,1))
}
