plot_timetrend <- function(sptsarfit,data,
                           time="year",spat="region",
                           xlab="Time",ylab="Time Trend",
                           title="Non-Parametric Time Trend by Spatial Unit"){
  # Function to plot time trend for each spatial unit.

  terms_spt_trend <- fit_terms(sptsarfit,"spttrend")
  spt_trend <- terms_spt_trend$fitted_terms
  se_spt_trend <- terms_spt_trend$se_fitted_terms
  data$spttrend <- spt_trend[,c("spttrend")]
  data$se_spttrend <- se_spt_trend[,c("spttrend")]
  data$spat_var <- as.factor(data[,colnames(data)==spat])
  data$time_var <- as.factor(data[,colnames(data)==time])
  spttrend_by_year <- NULL
  for (i in 1:length(levels(data$spat_var)))
  {
    name_reg <- levels(data$spat_var)[i]
    spttrend_reg <- data$spttrend[data$spat_var==name_reg]
    spttrend_by_year <- cbind(spttrend_by_year,spttrend_reg)
  }
  colnames(spttrend_by_year) <- levels(data$spat_var)
  matplot(levels(data$time_var),spttrend_by_year,type="l",
          xlab=xlab,ylab=ylab,
          main=title)
  res <- list(data=data,
              terms_spt_trend=terms_spt_trend)
  res

}
