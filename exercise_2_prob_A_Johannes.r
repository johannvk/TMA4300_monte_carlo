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
  plot(coal.df$date, cumulative_disasters, type='l', col="blue",
       xlab="Time (Years)", 
       ylab="Cumulative Coal Disasters (10 or more deaths)", 
       main="Overview of Coal Mining Disasters in the UK (Jarrett, 1979)",
       ylim=c(0L, length(coal.df$date)*1.1)
       )
}

beta.sample = function(lambda0, lambda1){
  # Sample from Inverse-Gamma distribution
  # with alpha = 4, beta = lambda_0 + lambda_1 + 1.
  alpha.IG = 4
  beta.IG = lambda0 + lambda1 + 1.0
  beta.value = 1.0/gamma.sample(1, alpha.IG, beta.IG)
  return(beta.value)
}

lambda0.sample = function(t1, t0, beta, coal.df) {
  # Sample from a Gamma distribution,
  # alpha.G = x(t1) - x(t0) + 2,
  # beta.G  = t1 - t0 + 1/beta
  alpha.G = cumulative.coal.disasters(t1, coal.df) - 
            cumulative.coal.disasters(t0, coal.df) + 2
  beta.G = t1 - t0 + 1.0/beta
  lambda0.value = gamma.sample(1, alpha.G, beta.G)
  return(lambda0.value)
}

lambda1.sample = function(t2, t1, beta, coal.df) {
  # Sample from full conditonal for lambda1:
  # A Gamma distribution with parameters
  # alpha.G = x(t2) - x(t1) + 2,
  # beta.G  = t2 - t1 + 1/beta
  alpha.G = cumulative.coal.disasters(t2, coal.df) - 
            cumulative.coal.disasters(t1, coal.df) + 2
  beta.G = t2 - t1 + 1.0/beta
  lambda1.value = gamma.sample(1, alpha.G, beta.G)
  return(lambda1.value)
}

t1.kernel = function(t1, lambda0, lambda1, coal.df, log.scale=T) {
  log.f = cumulative.coal.disasters(t1, coal.df)*log(lambda0/lambda1)
  log.f = log.f - (lambda0 - lambda1)*t1
  if (log.scale) return(log.f)
  else return(exp(log.f))
}

t1.RW.Metropolis.sample = function(t1, sigma2.t, lambda0, lambda1, coal.df) {
  t1.star = normal.sample(1, t1, sigma2.t)

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

    lambda0 = lambda0.sample(t1, t.start, beta, coal.df)
    lambda1 = lambda1.sample(t.end, t1, beta, coal.df)
    beta = beta.sample(lambda0, lambda1)
    
    # Store results:
    t1.vec[[i]]      = t1
    lambda0.vec[[i]] = lambda0
    lambda1.vec[[i]] = lambda1
    beta.vec[[i]]    = beta
  }

  t1.acceptance.rate = num.accepted.t1/N
  print(paste("t1 acceptance rate:", t1.acceptance.rate))

  ret.list = list(t1.vec, lambda0.vec, lambda1.vec, beta.vec)
  names(ret.list) = c("t1", "lambda0", "lambda1", "beta")
  return(ret.list)
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
       xlim=c(1885, 1900), ylim=c(0.0, 0.45),
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
  
  par(init)
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
  return(gibbs.mcmc.res)
}

main_A = function() {
  
  normalize_time = F
  coal.df = boot::coal
  
  if (normalize_time){
    # Makes numerics nicer, but will affect 
    # the Poisson process rates.
    translate = min(coal.df)
    scale = max(coal.df) - min(coal.df)
    coal.df = (coal.df - translate)/scale
  }

  test.t = 1901.6325
  num.accidents = cumulative.coal.disasters(test.t, coal.df)
  print(paste("Num of accidents before time", test.t, ":", num.accidents))
  
  plot.coal.disasters(coal.df)
  
  # Best acceptance rate for 'unif': sigma2.t = 100
  # Best acceptance rate for 'norm': sigma2.t = 25
  gibbs.hybrid.res = run.hybrid.Gibbs.MCMC(N=1.0e4L, sigma2.t = 25.0, 
                                           t1.prop="norm")
  plot.parameter.results(gibbs.hybrid.res)
  return(gibbs.hybrid.res)
}

samples.hpd = function(samples, prob=0.95) {
  attr(samples, "class") = "mcmc"
  samples.hpd = coda::HPDinterval(samples, prob=prob)
  dimnames(samples.hpd)[[1]][1] = c("HPD")
  return(samples.hpd)
}

MCMC.res = main_A()
