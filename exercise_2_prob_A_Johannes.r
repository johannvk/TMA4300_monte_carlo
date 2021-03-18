library(boot)
library(ggplot2)
library(coda)
library(forecast)

# Using our own functions from previous exercises:
# Sampling from Normal distribution:
source("exercise_1_prob_A.r")

# Sampling from Gamma distribution:
source("exercise_1_prob_B.r")


cumulative.coal.disasters = function(t, coal.df){
  # Find the cumulative number of disasters before time t.

  # First and last time in the dataset are not disasters,
  # but the 'start'- and 'end'-time for the dataset.
  tot_num_disasters = length(coal.df$date) - 2

  if (t < coal.df$date[[2]]) {
    return(0)
  } 
  index_of_first_time_after_t = Position(function(x) x > t, coal.df$date, 
                                         nomatch=NA)
  # Did not find a time in the dataset greater than 't':
  if (is.na(index_of_first_time_after_t)){
    return(tot_num_disasters)
  }
  # Return number of disasters before this time:  
  num_disasters = index_of_first_time_after_t - 2
  return(num_disasters)  
}

plot.coal.disasters = function(coal.df) {
  cumulative_disasters = mapply(function(t) 
                                cumulative.coal.disasters(t, coal.df), 
                                coal.df$date, SIMPLIFY=T)
  par(mfrow = c(1, 1))
  plot(coal.df$date, cumulative_disasters, type='l', col="blue",
       xlab="Time (Years)", 
       ylab="Cumulative Coal Disasters (10 or more deaths)", 
       main="Overview of Coal Mining Disasters in the UK (Jarrett, 1979)",
       ylim=c(0L, length(coal.df$date)*1.05), 
       # xaxs="i", yaxs="i",
       # panel.first=  grid(nx=16, ny=20, lty=1)
       )
  x_ticks <- axis(1, labels = FALSE)
  y_ticks <- axis(2, labels = FALSE)
  abline(v = x_ticks, h = y_ticks, lwd = 2.5, 
         lty = 3, col = "lightgray")
}

beta.sample = function(lambda0, lambda1){
  # Sample from Inverse-Gamma distribution
  # with alpha = 4, beta = lambda_0 + lambda_1 + 1.
  alpha.IG = 4
  beta.IG = lambda0 + lambda1 + 1.0
  beta.value = 1.0/gamma.sample(1, alpha.IG, beta.IG)
  return(beta.value)
}

lambda0.sample = function(t1, beta, coal.df) {
  # Sample from a Gamma distribution,
  # alpha.G = x(t1) - x(t0) + 2,
  # beta.G  = t1 - t0 + 1/beta
  t0 = coal.df$date[[1]]
  alpha.G = cumulative.coal.disasters(t1, coal.df) - 
            cumulative.coal.disasters(t0, coal.df) + 2
  beta.G = t1 - t0 + 1.0/beta
  lambda0.value = gamma.sample(1, alpha.G, beta.G)
  return(lambda0.value)
}

lambda1.sample = function(t1, beta, coal.df) {
  # Sample from full conditonal for lambda1:
  # A Gamma distribution with parameters
  # alpha.G = x(t2) - x(t1) + 2,
  # beta.G  = t2 - t1 + 1/beta
  t2 = coal.df$date[[length(coal.df$date)]]
  alpha.G = cumulative.coal.disasters(t2, coal.df) - 
            cumulative.coal.disasters(t1, coal.df) + 2
  beta.G = t2 - t1 + 1.0/beta
  lambda1.value = gamma.sample(1, alpha.G, beta.G)
  return(lambda1.value)
}

t1.kernel = function(t1, lambda0, lambda1, coal.df, log.scale=T) {
  x_t1 = cumulative.coal.disasters(t1, coal.df)
  y0 = x_t1; y1 = 189 - x_t1
  
  log.f = (y0 + 1)*log(lambda0) + (y1 + 1)*log(lambda1)
  log.f = log.f + (lambda1 - lambda0)*t1
  if (log.scale) return(log.f)
  else return(exp(log.f))
}

t1.array.kernel = function(t1.array, lambda0, lambda1, coal.df, log.scale=T) {
  res = double(length(t1.array))
  for (i in 1:length(t1.array)) {
      res[[i]] = t1.kernel(t1.array[[i]], lambda0, lambda1, coal.df,
                           log.scale = log.scale)
  }
  return(res)
}

plot.t1.kernel = function(lam0, lam1, t1.old=1906, alpha=F) {
  coal.df = boot::coal
  alpha.t1.post = function(t) {
       t1.array.kernel(t, lam0, lam1, coal.df = coal.df) - 
       t1.array.kernel(t1.old, lam0, lam1, coal.df = coal.df)
  }
  
  t1.post.dens = function(t) t1.array.kernel(t, lam0, lam1, coal.df = coal.df)
  
  if (alpha) {
    curve(alpha.t1.post, from=coal.df$date[[1]], 
          to=coal.df$date[[length(coal.df$date)]], 
          ylab="log(alpha.t1)", xlab="Years")
  } else {
    curve(t1.post.dens, from=coal.df$date[[1]], 
          to=coal.df$date[[length(coal.df$date)]], 
          ylab="log(t1.kernel)")
  }
}

t1.RW.Metropolis.sample = function(t1, sigma2.t, lambda0, lambda1, coal.df) {
  t.start = coal.df$date[[1]]
  t.end = coal.df$date[[length(coal.df$date)]]

  t1.star = normal.sample(1, t1, sigma2.t)

  log.f.prop = t1.kernel(t1.star,  lambda0, lambda1, coal.df=coal.df)
  log.f.old  = t1.kernel(t1,       lambda0, lambda1, coal.df=coal.df)

  alpha = exp(log.f.prop - log.f.old)
  accepted = runif(1) < alpha
  viable = (t.start <= t1.star) && (t1.star <= t.end)
  if (accepted && viable) { t1.value = t1.star }
  else { t1.value = t1; accepted = F }
  
  ret.list = list(t1.value, accepted)
  names(ret.list) = c("value", "accepted")
  return(ret.list)
}

t1.Unif.Metropolis.sample = function(t1, sigma2.t, lambda0, lambda1, coal.df) {
  t.start = coal.df$date[[1]]
  t.end = coal.df$date[[length(coal.df$date)]]
  
  
  t1.star = runif(1, t1 - sqrt(sigma2.t), t1 + sqrt(sigma2.t))
  
  log.f.prop = t1.kernel(t1.star,  lambda0, lambda1, coal.df=coal.df)
  log.f.old  = t1.kernel(t1,       lambda0, lambda1, coal.df=coal.df)
  
  alpha = exp(log.f.prop - log.f.old)
  accepted = runif(1) < alpha
  if (accepted) { t1.value = t1.star }
  else { t1.value = t1 }
  
  ret.list = list(t1.value, accepted)
  names(ret.list) = c("value", "accepted")
  return(ret.list)
}

hybrid.Gibbs.MCMC.coal.sampling = function(N, t1, lambda0, lambda1, beta,
                                           sigma2.t, coal.df, t1.prop="norm") {
  # Generate a Markov chain of length 'N' for the four parameters 
  # 't1', 'lambda0', 'lambda1', and 'beta'.
  
  t.start = coal.df$date[[1]]
  t.end = coal.df$date[[length(coal.df$date)]]

  # Initialize storage vectors:
  t1.vec = double(N)
  lambda0.vec = double(N)
  lambda1.vec = double(N)
  beta.vec = double(N)
  
  num.accepted.t1 = 0 

  for (i in 1:N) {
    if (t1.prop == "norm") {
      t1.list = t1.RW.Metropolis.sample(t1, sigma2.t, lambda0, 
                                        lambda1, coal.df)
    } else if (t1.prop == "unif") {
      t1.list = t1.Unif.Metropolis.sample(t1, sigma2.t, lambda0, 
                                          lambda1, coal.df)
    }
    t1 = t1.list$value
    num.accepted.t1 = num.accepted.t1 + t1.list$accepted 

    lambda0 = lambda0.sample(t1, beta, coal.df)
    lambda1 = lambda1.sample(t1, beta, coal.df)
    beta = beta.sample(lambda0, lambda1)
    
    # Store results:
    t1.vec[[i]]      = t1
    lambda0.vec[[i]] = lambda0
    lambda1.vec[[i]] = lambda1
    beta.vec[[i]]    = beta
  }

  t1.acceptance.rate = num.accepted.t1/N

  ret.list = list(t1.vec, lambda0.vec, lambda1.vec, beta.vec, 
                  t1.acceptance.rate, sigma2.t)
  names(ret.list) = c("t1", "lambda0", "lambda1", "beta", 
                      "t1.accept.rate", "sigma2.t")
  return(ret.list)
}

first.block.update = function(params, sigma2.t, coal.df) {
  t.start = coal.df$date[[1]]
  t.end = coal.df$date[[length(coal.df$date)]]

  tot.disasters = length(coal.df$date) - 2
  t1 = params[["t1"]]; beta = params[["beta"]]
  lam0 = params[["lambda0"]]; lam1 = params[["lambda1"]]

  y0 = cumulative.coal.disasters(t1, coal.df)
  y1 = tot.disasters - y0
  
  # Calculate Gamma parameters for the full conditional 
  # distributions of lambda0 and lambda1.
  lam0.alpha = y0 + 2; lam0.beta = t1 - t.start + 1.0/beta
  lam1.alpha = y1 + 2; lam1.beta = t.end - t1 + 1.0/beta
  
  # Sample new values:
  t1.star = normal.sample(1, t1, sigma2.t)
  
  lam0.star = lambda0.sample(t1.star, beta, coal.df)
  lam1.star = lambda1.sample(t1.star, beta, coal.df)
  
  # Calculate Gamma parameters for the full conditional 
  # distributions of lambda0.star and lambda1.star
  y0.star = cumulative.coal.disasters(t1.star, coal.df)
  y1.star = tot.disasters - y0.star
  
  # Uncomment: 
  lam0.star.alpha = y0.star + 2; lam0.star.beta = t1.star - t.start + 1.0/beta
  lam1.star.alpha = y1.star + 2; lam1.star.beta = t.end - t1.star + 1.0/beta
  
  # Calculate acceptance probability 'alpha': 
  log.numerator  = lam0.alpha*log(lam0.beta) +
                   lgamma(lam0.star.alpha) +
                   lam1.alpha*log(lam1.beta) + 
                   lgamma(lam1.star.alpha)

  log.denominator = lam0.star.alpha*log(lam0.star.beta) +
                    lgamma(lam0.alpha) +
                    lam1.star.alpha*log(lam1.star.beta) + 
                    lgamma(lam1.alpha)
  
  log.alpha = log.numerator - log.denominator

  # Restrict the movement of t1-parameter.
  viable.time = (t.start <= t1.star) && (t1.star <= t.end)
  accepted = runif(1) < exp(log.alpha)

  if (!is.na(accepted) && viable.time && accepted) {
    params[["t1"]] = t1.star
    params[["lambda0"]] = lam0.star
    params[["lambda1"]] = lam1.star
  } else { accepted = F }

  res.list = list(params, accepted)
  names(res.list) = c("params", "accepted")
  return(res.list)
}

second.block.update = function(params, sigma2.beta, coal.df) {
  t.start = coal.df$date[[1]]
  t.end = coal.df$date[[length(coal.df$date)]]
  
  t1 = params[["t1"]]; beta = params[["beta"]]
  lam0 = params[["lambda0"]]; lam1 = params[["lambda1"]]

  tot.disasters = length(coal.df$date) - 2
  y0 = cumulative.coal.disasters(t1, coal.df)
  y1 = tot.disasters - y0
  
  # Sample new beta value, capped at the lowest value of 0.0.
  beta.star = normal.sample(1, beta, sigma2.beta)
  if (beta.star < 0.0) beta.star = 0.0
  
  lam0.star = lambda0.sample(t1, beta.star, coal.df)
  lam1.star = lambda1.sample(t1, beta.star, coal.df)
  
  # Calculate Gamma parameters for the full conditional 
  # distributions of lambda0 and lambda1. Before and after sampling.
  lam0.alpha = y0 + 2; lam0.beta = t1 - t.start + 1.0/beta
  lam1.alpha = y1 + 2; lam1.beta = t.end - t1 + 1.0/beta
  
  # Calculate acceptance probability 'alpha': 
  lam0.star.alpha = y0 + 2; lam0.star.beta = t1 - t.start + 1.0/beta.star
  lam1.star.alpha = y1 + 2; lam1.star.beta = t.end - t1 + 1.0/beta.star
  
  log.prop.numerator  = lam0.alpha*log(lam0.beta) +
                        lgamma(lam0.star.alpha) +
                        lam1.alpha*log(lam1.beta) + 
                        lgamma(lam1.star.alpha)
  
  log.prop.denominator = lam0.star.alpha*log(lam0.star.beta) +
                         lgamma(lam0.alpha) +
                         lam1.star.alpha*log(lam1.star.beta) + 
                         lgamma(lam1.alpha)
  # Correct!
  log.alpha = 5*(log(beta) - log(beta.star)) + 1/beta - 1/beta.star +
              log.prop.numerator - log.prop.denominator
  
  viable = 0 <= beta.star
  accepted = (runif(1) < exp(log.alpha)) 
  if (!is.na(accepted) && accepted && viable) {
    params["lambda0"] = lam0.star
    params["lambda1"] = lam1.star
    params["beta"] = beta.star
  } else { accepted = F }
  
  
  res.list = list(params, accepted)
  names(res.list) = c("params", "accepted")
  return(res.list)
}

block1.MCMC.coal.sampling = function(N, params, sigma2.t, sigma2.beta, coal.df) {
  t1.vec = double(N)
  lambda0.vec = double(N)
  lambda1.vec = double(N)
  beta.vec = double(N)
  
  num.first.block.accepted = 0
  
  for (i in 1:N) {
    # First block update:
    block.1 = first.block.update(params, sigma2.t, coal.df) 
    params = block.1$params
    
    # Gibbs-Update for beta:
    params[["beta"]] = beta.sample(params[["lambda0"]], params[["lambda1"]])

    t1.vec[[i]]      = params[["t1"]]
    lambda0.vec[[i]] = params[["lambda0"]]
    lambda1.vec[[i]] = params[["lambda1"]]
    beta.vec[[i]]    = params[["beta"]]
    
    num.first.block.accepted = num.first.block.accepted +
                               block.1$accepted
  }
  first.accept.rate = num.first.block.accepted/N

  ret.list = list(t1.vec, lambda0.vec, lambda1.vec, beta.vec, 
                  first.accept.rate)
  names(ret.list) = c("t1", "lambda0", "lambda1", "beta", 
                      "first.accept.rate")
  return(ret.list)  
}

block2.MCMC.coal.sampling = function(N, params, sigma2.t, sigma2.beta, coal.df) 
{
  t1.vec = double(N)
  lambda0.vec = double(N)
  lambda1.vec = double(N)
  beta.vec = double(N)
  
  num.t1.accepted = 0
  num.second.block.accepted = 0
  
  for (i in 1:N) {
    
    # block.1 = first.block.update(params, sigma2.t, coal.df) 
    # params = block.1$params
    
    # Block update:
    block.2 = second.block.update(params, sigma2.beta, coal.df)
    params = block.2$params
    num.second.block.accepted = num.second.block.accepted +
                            block.2$accepted
    # Gibbs-Update for t1:
    t1.star = t1.RW.Metropolis.sample(params[["t1"]], sigma2.t,
                  params[["lambda0"]], params[["lambda1"]], coal.df)
    num.t1.accepted = num.t1.accepted + t1.star$accepted
    params[["t1"]] = t1.star$value
    
    t1.vec[[i]]      = params[["t1"]]
    lambda0.vec[[i]] = params[["lambda0"]]
    lambda1.vec[[i]] = params[["lambda1"]]
    beta.vec[[i]]    = params[["beta"]]
    
  }
  t1.accept.rate = num.t1.accepted/N
  second.accept.rate = num.second.block.accepted/N
  
  ret.list = list(t1.vec, lambda0.vec, lambda1.vec, beta.vec, 
                  t1.accept.rate, second.accept.rate)
  names(ret.list) = c("t1", "lambda0", "lambda1", "beta", 
                      "t1.accept.rate", "second.accept.rate")
  return(ret.list)  
}

two.block.coal.MCMC = function(N, params, sigma2.t, sigma2.beta, coal.df) {
  t1.vec = double(N)
  lambda0.vec = double(2*N)
  lambda1.vec = double(2*N)
  beta.vec = double(N)
  
  num.first.block.accepted = 0
  num.second.block.accepted = 0
  
  for (i in 1:N) {
    # First block update:
    block.1 = first.block.update(params, sigma2.t, coal.df) 
    params = block.1$params
    num.first.block.accepted = num.first.block.accepted +
                               block.1$accepted
    
    # Store results from block 1:
    t1.vec[[i]]      = params[["t1"]]
    lambda0.vec[[2*i - 1]] = params[["lambda0"]]
    lambda1.vec[[2*i - 1]] = params[["lambda1"]]
    
    # Second block update:
    block.2 = second.block.update(params, sigma2.beta, coal.df)
    params = block.2$params
    num.second.block.accepted = num.second.block.accepted +
                                block.2$accepted
    
    # Store results from block 2:
    lambda0.vec[[2*i]] = params[["lambda0"]]
    lambda1.vec[[2*i]] = params[["lambda1"]]
    beta.vec[[i]]    = params[["beta"]]
  }

  first.accept.rate = num.first.block.accepted/N
  second.accept.rate = num.second.block.accepted/N
  
  ret.list = list(t1.vec, lambda0.vec, lambda1.vec, beta.vec, 
                  first.accept.rate, second.accept.rate,
                  sigma2.t, sigma2.beta)
  names(ret.list) = c("t1", "lambda0", "lambda1", "beta", 
                      "first.accept.rate", "second.accept.rate",
                      "sigma2.t", "sigma2.beta")
  return(ret.list)  
}

run.hybrid.Gibbs.MCMC = function(N=1.0e4L, normalize.time=F, sigma2.t=100.0,
                                 t1.prop="norm") {
  normalize_time = F
  coal.df = boot::coal
  
  if (normalize.time){
    # Makes numerics nicer, but will affect 
    # the Poisson process rates.
    translate = min(coal.df)
    scale = max(coal.df) - min(coal.df)
    coal.df = (coal.df - translate)/scale
  }
  
  t.start = coal.df$date[[1]]
  t.end = coal.df$date[[length(coal.df$date)]]   
  
  ######### Run MCMC chain with Hybrid Gibbs sampler #########
  # Initial values:
  t1 = (t.end + t.start)/2.0
  lambda0 = 2.0
  lambda1 = 2.0
  beta = 1.0
  
  # Hyper-parameter:
  # sigma2.t = 100.0
  
  gibbs.mcmc.res = hybrid.Gibbs.MCMC.coal.sampling(N, 
                                                   t1, lambda0, lambda1, beta, 
                                                   sigma2.t, coal.df, t1.prop)
  print(paste("t1 acceptance rate:", gibbs.mcmc.res$t1.accept.rate))
  return(gibbs.mcmc.res)
}

run.block.MCMC = function(sigma2.t = 9.0, sigma2.beta = 4.0,
                          params = c(1870, 1.0, 3.0, 1.0), N=1.0e4L,
                          normalize_time = F) {
  names(params) = c("t1", "lambda0", "lambda1", "beta")
  
  coal.df = boot::coal
  if (normalize_time){
    # Makes numerics nicer, but will affect 
    # the Poisson process rates.
    translate = min(coal.df)
    scale = max(coal.df) - min(coal.df)
    coal.df = (coal.df - translate)/scale
  }
  
  block.res = two.block.coal.MCMC(N, params, sigma2.t, 
                                  sigma2.beta, coal.df)
  print(paste("First block acceptance rate:", block.res$first.accept.rate))
  print(paste("Second block acceptance rate:", block.res$second.accept.rate))
  return(block.res)
}

plot.parameter.results = function(mcmc.result) {
  
  init = par(no.readonly=TRUE)
  par(mfrow = c(1, 3))
  
  t1.res = mcmc.result$t1; beta.res = mcmc.result$beta
  lambda0.res = mcmc.result$lambda0; lambda1.res = mcmc.result$lambda1
    
  # Time results:
  print("Result Mosaic for the 'splitting time' t1.")
  plot(seq_along(t1.res), t1.res, type="l", col="red",
       main="t1: Trace plot", ylab="MCMC Samples of t1", xlab="Step")
  hist(t1.res, breaks=100, probability = T, main="Histogram Density", 
       ylab="Density", xlab="Years", 
       # xlim=c(1, 1900), ylim=c(0.0, 0.45),
       col="blue")
  forecast::Acf(t1.res, lag.max = 10, main="ACF for t1")
  
  # Lambda0 results:
  print("Result Mosaic for the first poisson process intensity, lambda0:")
  plot(seq_along(lambda0.res), lambda0.res, type="l", col="red",
       main="lambda0: Trace plot", ylab="MCMC Samples of lambda0", xlab="Step")
  hist(lambda0.res, breaks=100, probability = T, main="Histogram Density", 
       ylab="Density", xlab="Number of accidents per year between [t0, t1]", 
       col="blue")
  forecast::Acf(lambda0.res, lag.max = 10, main="ACF for lambda0")
  
  # Lambda1 results:
  print("Result Mosaic for the second poisson process intensity, lambda1:")
  plot(seq_along(lambda1.res), lambda1.res, type="l", col="red",
       main="lambda1: Trace plot", ylab="MCMC Samples of lambda1", xlab="Step")
  hist(lambda1.res, breaks=100, probability = T, main="Histogram Density", 
       ylab="Density", xlab="Number of accidents per year between [t1, t2]", 
       col="blue")
  forecast::Acf(lambda1.res, lag.max = 10, main="ACF for lambda0")
    
  # Beta results:
  print("Result Mosaic for the hyper-parameter, beta:")
  plot(seq_along(beta.res), beta.res, type="l", col="red",
       main="beta: Trace plot", ylab="MCMC Samples of lambda1", xlab="Step")
  hist(beta.res, breaks=200, probability = T, main="Histogram Density", 
       ylab="Density", xlab="scale parameter of Poisson process' intensity", 
       xlim=c(0, 10),
       col="blue")
  forecast::Acf(beta.res, lag.max = 10, main="ACF for beta")

  # Possible args to change font sizes:
  # cex.lab = 2,
  # cex.axis = 2,
  # cex.main = 2,
  
  par(mfrow = c(1, 1))
}

samples.hpd = function(samples, prob=0.95) {
  attr(samples, "class") = "mcmc"
  samples.hpd = coda::HPDinterval(samples, prob=prob)
  dimnames(samples.hpd)[[1]][1] = c("HPD")
  return(samples.hpd)
}

mcmc.diagnostics = function(MCMC.res) {
  t1.mcmc = coda::mcmc(MCMC.res$t1)
  densplot(t1.mcmc, main="Estimated posterior density for t1", 
           ylab="Density", xlab="Years (+ observations)",
           xlim=c(1884.5, 1900))

  lam0.mcmc = coda::mcmc(MCMC.res$lambda0)
  densplot(lam0.mcmc, main="Estimated posterior density for lambda0", 
           ylab="Density", xlab="Disasters/year (+ observations)")

  lam1.mcmc = coda::mcmc(MCMC.res$lambda1)
  densplot(lam1.mcmc, main="Estimated posterior density for lambda1", 
           ylab="Density", xlab="Disasters/year (+ observations)")
  
  beta.mcmc = coda::mcmc(MCMC.res$beta)
  densplot(beta.mcmc, main="Estimated posterior density for beta", 
           ylab="Density", xlab="Disasters/year (+ observations)",
           xlim=c(0, 6.2))
}

running.mean = function(data, main) {
  run.mean = cumsum(data)/1:length(data)
  plot(1:length(data), run.mean, ylab="Running Mean", 
       xlab="Iterations", type="l", main=main)
}

trace.plot = function(data, main) {
  plot(seq_along(data), data, type="l", col="red",
       main=main, ylab="MCMC Samples", xlab="Step")
}

facet.trace = function(MCMC.res, gibbs=T) {
  layout_matr = matrix(c(1, 2, 4, 3), nrow=2, ncol=2)
  layout(layout_matr)
  if (gibbs) {
    prefix = "Gibbs:"
    sigma2.t = MCMC.res$sigma2.t
    trace.plot(MCMC.res$t1, main=paste(prefix, "t1,", "sigma2.t =", sigma2.t))
    trace.plot(MCMC.res$lambda0, main=paste(prefix, "lam0"))
    trace.plot(MCMC.res$lambda1, main=paste(prefix, "lam1"))
    trace.plot(MCMC.res$beta, main=paste(prefix, "beta"))
  } else {
    prefix = "Block:"
    sigma2.t = MCMC.res$sigma2.t; sigma2.beta = MCMC.res$sigma2.beta
    trace.plot(MCMC.res$t1, main=paste(prefix, "t1,", "sigma2.t =", sigma2.t))
    trace.plot(MCMC.res$lambda0, main=paste(prefix, "lam0"))
    trace.plot(MCMC.res$lambda1, main=paste(prefix, "lam1"))
    trace.plot(MCMC.res$beta, main=paste(prefix, "beta,", 
                                           "sigma2.beta =", sigma2.beta))
  }
}

facet.running.means = function(MCMC.res, gibbs=T) {
  
  params = data.frame(t1 = MCMC.res$t1, lambda0 = MCMC.res$lambda0,
                      lambda1 = MCMC.res$lambda1, beta=MCMC.res$beta)
  layout_matr = matrix(c(1, 2, 4, 3), nrow=2, ncol=2)
  layout(layout_matr) 
  if (gibbs) {
    prefix = "Gibbs:"
    sigma2.t = MCMC.res$sigma2.t
    running.mean(MCMC.res$t1, main=paste(prefix, "t1,", "sigma2.t =", sigma2.t))
    running.mean(MCMC.res$lambda0, main=paste(prefix, "lam0"))
    running.mean(MCMC.res$lambda1, main=paste(prefix, "lam1"))
    running.mean(MCMC.res$beta, main=paste(prefix, "beta"))
  } else {
    prefix = "Block:"
    sigma2.t = MCMC.res$sigma2.t; sigma2.beta = MCMC.res$sigma2.beta
    running.mean(MCMC.res$t1, main=paste(prefix, "t1,", "sigma2.t =", sigma2.t))
    running.mean(MCMC.res$lambda0, main=paste(prefix, "lam0"))
    running.mean(MCMC.res$lambda1, main=paste(prefix, "lam1"))
    running.mean(MCMC.res$beta, main=paste(prefix, "beta,", 
                                           "sigma2.beta =", sigma2.beta))
  }
}



main_A = function() {
  set.seed(2)
  
  normalize_time = F
  coal.df = boot::coal
  
  if (normalize_time){
    # Makes numerics nicer, but will affect
    # the Poisson process rates.
    translate = min(coal.df)
    scale = max(coal.df) - min(coal.df)
    coal.df = (coal.df - translate)/scale
  }

  # plot.coal.disasters(coal.df)
  
  # Best acceptance rate for 'unif': sigma2.t = 100
  # Best acceptance rate for 'norm': sigma2.t = 25
  N_gibbs = 1.0e4L
  gibbs.hybrid.res = run.hybrid.Gibbs.MCMC(N=N_gibbs, sigma2.t = 25.0, 
                                           t1.prop="norm")

  # plot.parameter.results(gibbs.hybrid.res)
  facet.running.means(gibbs.hybrid.res, gibbs=T)
  facet.trace(gibbs.hybrid.res, gibbs=T)
  
  N_block = 1.0e5L
  block.res = run.block.MCMC(sigma2.beta=16.0, N=N_block)

  facet.running.means(block.res, gibbs=F)
  facet.trace(block.res, gibbs=F)
  
  print(paste("Total samples:", N_block))
  print(paste("Effective sample size, t1, blocking:",
               coda::effectiveSize(block.res$t1)))
  print(paste("Effective sample size, beta, blocking:",
              coda::effectiveSize(block.res$beta)))
  print(paste("Effective sample size, lambda0, blocking:",
              coda::effectiveSize(block.res$lambda0)))
  print(paste("Effective sample size, lambda1, blocking:",
              coda::effectiveSize(block.res$lambda1)))
}

main_A()
# MCMC.res = main_A()
# t1.res = MCMC.res$t1
# mcmc.list.t1.res = coda::as.mcmc.list(t1.res, start=1, end=length(t1.res))
# start(mcmc.t1.res)
# mc1 = coda::geweke.plot(coda::as.mcmc(t1.res))
