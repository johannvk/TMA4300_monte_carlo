library(forecast)
library(ggplot2)
library(latex2exp)

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
    actual_residuals = ls.resid[bootstrap.resid[ , i]]
    bootstrap.ts = ARp.filter(x0, beta.ls, actual_residuals)[-(1:p)]
    
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
    actual_residuals = la.resid[bootstrap.resid[ , i]]
    bootstrap.ts = ARp.filter(x0, beta.la, actual_residuals)[-(1:p)]
    
    # Calculate LS beta coeffiecients on that series:
    betas = ARp.beta.est(bootstrap.ts, p)$LA
    
    # Store:
    bootstrap.betas[i, ] = betas
  }
  
  return (bootstrap.betas)
}

one.step.bootstrap = function(time.series, B=2000, p=2, 
                              include.noise=T, norm=c("LA", "LS")) {
  
  n = length(time.series)
  last.p.values = time.series[(n - p + 1):n]
  
  # betas = c(beta1, beta2), e ~ WN(sigma^2)
  if (include.noise) {
    next.value = function(betas, e) return(sum(last.p.values*rev(betas)) + e)
  } else {
    next.value = function(betas, e) return(sum(last.p.values*rev(betas)))
  }
  
  betas.full = ARp.beta.est(time.series, p)[[norm]] 
  
  # Calculate LS residuals:
  # Has length (n - p), assumed i.i.d:
  full.resid = ARp.resid(time.series, betas.full)
  
  # Generate n i.i.d residuals indices down B columns:
  bootstrap.resid = matrix(sample(1:(n-p), B*n, replace=T), 
                           nrow=n, ncol=B)
  
  # Where to start bootstrap time-series from:
  start.indices = sample(1:(n - p + 1), B, replace=T)
  
  # Which i.i.d. residual to use when time-stepping. 
  iid.indices = sample(1:(n - p + 1), B, replace=T)
  
  # Store the estimated step (beta0, beta1)
  # down the rows:
  bootstrap.next.step = double(B)  
  
  for (i in 1:B) {
    # Generate an AR(p) time series of length n:
    # Use the originally estimated "betas.full" coefficients
    # to "drive" the new time series:
    x0 = time.series[start.indices[i]:(start.indices[i] + p - 1)]
    
    # Filter 'n' steps forward, and then remove the initial p steps:
    actual_residuals = full.resid[bootstrap.resid[ , i]]
    bootstrap.ts = ARp.filter(x0, betas.full, actual_residuals)[-(1:p)]
    
    # Calculate beta coeffiecients on that series:
    bootstrap.betas = ARp.beta.est(bootstrap.ts, p)[[norm]]
    
    # Find new bootstrap residulas, length = (n - p):
    bootstrap.ts.resid = ARp.resid(bootstrap.ts, bootstrap.betas) 
    
    # Calculate the next step with a random noise term,
    # and the bootstrapped beta-coefficients:
    noise.i = bootstrap.ts.resid[iid.indices[i]]
    bootstrap.next.step[i] = next.value(bootstrap.betas, noise.i)
  }
  #Returns bootstrapped one step evolution
  return (bootstrap.next.step[!is.na(bootstrap.next.step)])
}


emp.quantile.interval = function(x, alpha) {
  # Generate a '(1-alpha)' confidence interval.
  lower.quant = alpha/2.0
  upper.quant = 1 - alpha/2.0
  
  return (quantile(x, probs=c(lower.quant, upper.quant)))
}
  

A1.main = function() {
  time.series = ts(data3A$x)
  plot(time.series)
  
  pTS <- ggplot(data.frame(ts=data3A$x, id=(1:100)), aes(y=ts,x=id)) +
    geom_line() + 
    xlab(TeX("$t$")) + 
    ylab(TeX("$x_t$")) +
    theme_bw()
  print(pTS)
  
  betas.full = ARp.beta.est(time.series, p=2)
  
  LS.betas = betas.full$LS
  LS.betas.bootstrap = beta.LS.bootstrap(time.series)
  LS.bias.est = colMeans(LS.betas.bootstrap) - LS.betas
  
  ls.resid.full = ARp.resid(time.series, LS.betas)
  #hist(ls.resid.full, breaks=50, probability = T, main="LS Residuals")
  
  pLSres <- ggplot(data.frame(res=ls.resid.full), aes(x=res)) + 
    geom_histogram(bins=40, aes(y=..count../sum(..count..))) + 
    xlab(TeX("LS $\\hat{\\epsilon}_t$")) + 
    ylab("Density") +
    labs(title = "LS Residuals") + 
    theme_bw()
  print(pLSres)
  
  cat("\nOriginal LS.betas:\n")
  print(LS.betas)
  
  cat("\nEstimated bias of LS betas-statistic:\n")
  print(LS.bias.est)
  
  cat(paste("\n", "Bootstrapped LS.betas covariance:\n", sep=""))
  print(cov(LS.betas.bootstrap))
  
  
  pLSbeta0 <- ggplot(data.frame(betas=LS.betas.bootstrap[ , 1]), aes(x=betas)) + 
    geom_histogram(bins=40, aes(y=..count../sum(..count..))) + 
    xlab(TeX("$\\beta_1$")) + 
    ylab("Density") +
    labs(title = TeX("LS $\\beta_1^*$")) + 
    theme_bw()
  print(pLSbeta0)
  pLSbeta1 <- ggplot(data.frame(betas=LS.betas.bootstrap[ , 2]), aes(x=betas)) + 
    geom_histogram(bins=40, aes(y=..count../sum(..count..))) + 
    xlab(TeX("$\\beta_2$")) + 
    ylab("Density") +
    labs(title = TeX("LS $\\beta_2^*$")) + 
    theme_bw()
  print(pLSbeta1)
  
  LA.betas.bootstrap = beta.LA.bootstrap(time.series)
  LA.betas = betas.full$LA
  
  LA.bias.est = colMeans(LA.betas.bootstrap) - LA.betas

  la.resid.full = ARp.resid(time.series, LA.betas)
  #hist(ls.resid.full, breaks=50, probability = T, main="LA Residuals")
  
  pLAres <- ggplot(data.frame(res=ls.resid.full), aes(x=res)) + 
    geom_histogram(bins=40, aes(y=..count../sum(..count..))) + 
    xlab(TeX("LA $\\hat{\\epsilon}_t$")) + 
    ylab("Density") +
    labs(title = "LA Residuals") + 
    theme_bw()
  print(pLAres)
  
  cat("\n\nOriginal LA.betas:\n")
  print(LA.betas)
  
  cat("\nEstimated bias of LA betas-statistic:\n")
  print(LA.bias.est)
  
  cat(paste("\n", "Bootstrapped LA.betas covariance:\n", sep=""))
  print(cov(LA.betas.bootstrap))
  
  pLAbeta0 <- ggplot(data.frame(betas=LA.betas.bootstrap[ , 1]), aes(x=betas)) + 
    geom_histogram(bins=40, aes(y=..count../sum(..count..))) + 
    xlab(TeX("$\\beta_1$")) + 
    ylab("Density") +
    labs(title = TeX("LA $\\beta_1^*$")) + 
    theme_bw()
  print(pLAbeta0)
  pLAbeta1 <- ggplot(data.frame(betas=LA.betas.bootstrap[ , 2]), aes(x=betas)) + 
    geom_histogram(bins=40, aes(y=..count../sum(..count..))) + 
    xlab(TeX("$\\beta_2$")) + 
    ylab("Density") +
    labs(title = TeX("LA $\\beta_2^*$")) + 
    theme_bw()
  print(pLAbeta1)
  
  
}


A2.main = function() {
  # Think we are meant to use 'newly computed'
  # residuals from each bootstrap time series
  # to estimate x_{101} for each bootstrap-run.
  
  # Get very (!) wide prediction intervals for 
  # x_{101} when doing that. Something seems wrong?

  set.seed(8)
  
  time.series = ts(data3A$x)
  
  print("LS:")
  next.steps.boot = one.step.bootstrap(time.series, B=3000, 
                                       include.noise=T, norm="LS")
  
  alpha = 0.05
  emp.conf.interval = emp.quantile.interval(next.steps.boot, alpha)
  
  cat("\n95% confidence interval for x_{101}:\n")
  print(emp.conf.interval)
  
  hist(next.steps.boot, breaks=100, probability = T, xlim=c(14, 18))
  
  #The coolest plot:
  x_101 = 1.5528106*data3A$x[100]-0.5680178*data3A$x[99]
  tsPred = c(data3A$x,x_101)
  print(paste("x_101",x_101))
  print(paste("mean",mean(next.steps.boot)))
  lowerConfint = c(data3A$x,emp.conf.interval[1])
  upperConfint = c(data3A$x,emp.conf.interval[2])
    
    pTS <- ggplot(data.frame(ts=tsPred, id=(1:101)), aes(y=ts,x=id)) +
    geom_line() + 
      geom_ribbon(data=data.frame(t=(1:101),lowerConfint=lowerConfint,upperConfint=upperConfint),
                  aes(x=t, ymin=lowerConfint,ymax=upperConfint),
                  inherit.aes = FALSE, alpha=0.3) +
    xlab(TeX("$t$")) + 
    xlim(80,101) +
    ylim(0,60) + 
    ylab(TeX("$x_t$")) +
    theme_bw()
  print(pTS)
  
  # geom_ribbon(data=mc.eta, aes(x=t, ymin=conf_int_lower,ymax=conf_int_upper),
  #             inherit.aes = FALSE, alpha=0.3)
  
  
  
  
  print("LA:")
  next.steps.boot = one.step.bootstrap(time.series, B=3000, 
                                       include.noise=T, norm="LA")
  
  alpha = 0.05
  emp.conf.interval = emp.quantile.interval(next.steps.boot, alpha)
  
  cat("\n95% confidence interval for x_{101}:\n")
  print(emp.conf.interval)

}

#A1.main()
A2.main()
