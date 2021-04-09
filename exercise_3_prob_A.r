library(forecast)

source("ex3-additional-files/probAdata.R")
source("ex3-additional-files/probAhelp.R")


beta.LS.bootstrap = function(ts, B=2000, p=2) {
  n = length(ts)
  
  # Calculate LS coefficients:
  beta.ls = ARp.beta.est(ts, p=p)$LS

  # Calculate LS residuals:
  # Has length (n - p), assumed i.i.d:
  ls.resid = ARp.resid(ts, beta.ls)
  
  # Generate n i.i.d residuals down B columns:
  bootstrap.resid = matrix(sample(1:(n-p), B*n, replace=T), 
                           nrow=n, ncol=B)

  start.indices = sample(1:(n - p + 1), B, replace=T)
  
  # Store the bootstrap coefficents (beta0, beta1)
  # down the rows:
  bootstrap.betas = matrix(NA, nrow=B, ncol=2)  
  
  # Could think of making a time series of length ceil(1.2*n)
  # and then only estimate the beta-coefficients on the last 
  # n elements, to hopefully avoid problems of non-stationarity.
  for (i in 1:B) {
    # Generate an AR(p) time series of length n:
    # Use the originally estimated "beta.ls" coefficients
    # to "drive" the new time series:
    x0 = ts[start.indices[i]:(start.indices[i] + p - 1)]
    
    # Filter 'n' steps forward, and then remove the initial p steps:
    bootstrap.ts = ARp.filter(x0, beta.ls, bootstrap.resid[ , i])[-(1:p)]
    
    # Calculate LS beta coeffiecients on that series:
    betas = ARp.beta.est(bootstrap.ts, p)$LS
    
    # Store:
    bootstrap.betas[i, ] = betas
  }
  
  return (bootstrap.betas)
}

beta.LA.bootstrap = function(ts, B=2000, p=2) {
  n = length(ts)
  
  # Calculate LA coefficients:
  beta.la = ARp.beta.est(ts, p=p)$LA
  
  # Calculate LS residuals:
  # Has length (n - p), assumed i.i.d:
  la.resid = ARp.resid(ts, beta.la)
  
  # Generate n i.i.d residuals down B columns:
  bootstrap.resid = matrix(sample(1:(n-p), B*n, replace=T), 
                           nrow=n, ncol=B)
  
  start.indices = sample(1:(n - p + 1), B, replace=T)
  
  # Store the bootstrap coefficents (beta0, beta1)
  # down the rows:
  bootstrap.betas = matrix(NA, nrow=B, ncol=2)  
  
  # Could think of making a time series of length ceil(1.2*n)
  # and then only estimate the beta-coefficients on the last 
  # n elements, to hopefully avoid problems of non-stationarity.
  for (i in 1:B) {
    # Generate an AR(p) time series of length n:
    # Use the originally estimated "beta.la" coefficients
    # to "drive" the new time series:
    x0 = ts[start.indices[i]:(start.indices[i] + p - 1)]
    
    # Filter 'n' steps forward, and then remove the initial p steps:
    bootstrap.ts = ARp.filter(x0, beta.la, bootstrap.resid[ , i])[-(1:p)]
    
    # Calculate LS beta coeffiecients on that series:
    betas = ARp.beta.est(bootstrap.ts, p)$LA
    
    # Store:
    bootstrap.betas[i, ] = betas
  }
  
  return (bootstrap.betas)
}
# 
A1.main = function() {
  time.series = ts(data3A$x)
  plot(time.series)
  
  betas.full = ARp.beta.est(time.series, p=2)
  LS.betas = betas.full$LS
  LS.betas.bootstrap = beta.LS.bootstrap(time.series)
  LS.bias.est = colMeans(LS.betas.bootstrap) - LS.betas
  
  cat("\nOriginal LS.betas:\n")
  print(LS.betas)
  
  cat("\nEstimated bias of LS betas-statistic:\n")
  print(LS.bias.est)
  
  cat(paste("\n", "Bootstrapped LS.betas covariance:\n", sep=""))
  print(cov(LS.betas.bootstrap))
  par(mfrow = c(1, 2))
  hist(LS.betas.bootstrap[ , 1], probability = T, main="LS.beta0")
  hist(LS.betas.bootstrap[ , 2], probability = T, main="LS.beta1")
  par(mfrow= c(1, 1))
  
  LA.betas.bootstrap = beta.LA.bootstrap(time.series)
  LA.betas = betas.full$LA
  LA.bias.est = colMeans(LA.betas.bootstrap) - LA.betas
  
  cat("\n\nOriginal LA.betas:\n")
  print(LA.betas)
  
  cat("\nEstimated bias of LA betas-statistic:\n")
  print(LA.bias.est)
  
  cat(paste("\n", "Bootstrapped LA.betas covariance:\n", sep=""))
  print(cov(LA.betas.bootstrap))
  par(mfrow = c(1, 2))
  hist(LA.betas.bootstrap[ , 1], probability = T, main="LA.beta0")
  hist(LA.betas.bootstrap[ , 2], probability = T, main="LA.beta1")
  par(mfrow= c(1, 1))
}

A1.main()