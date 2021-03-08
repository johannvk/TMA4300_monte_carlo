library(boot)
library(ggplot2)

# Using our own functions from previous exercises:

# Sampling from Normal distribution:
source("exercise_1_prob_A.r")

# Sampling from Gamma distribution:
source("exercise_1_prob_B.r")



cumulative.coal.disasters = function(t, coal.df){
  # Find the cumulative number of disasters before time t.
  # User decides whether or not to normalize the time of 
  # disasters to the interval [0, 1]

  # First and last time in the dataset are not disasters,
  # but the 'start'- and 'end'-time for the dataset.
  tot_num_disasters = length(coal.df$date) - 2

  if (t < coal.df$date[[2]]) {
    return(0)
  } 
  index_of_first_time_after_t = Position(function(x) x >= t, coal.df$date, 
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
  log.f.prop = t1.kernel(t1.star, lambda0, lambda1, coal.df)
  log.f.old = t1.kernel(t1, lambda0, lambda1, coal.df)
  alpha = exp(log.f.prop - log.f.old)
  
  if (runif(1) < alpha){
    return(t1.star)
  }
  else {
    return(t1)
  }
}

hybrid.Gibbs.MCMC.coal.sampling = function(N, t1, lambda0, lambda1, beta,
                                           sigma2.t, coal.df) {
  # Generate a Markov chain of length 'N' for the four parameters 
  # 't1', 'lambda0', 'lambda1', and 'beta'.
  
  t.start = coal.df$date[[1]]
  t.end = coal.df$date[[length(coal.df$date)]]

  # Initialize storage vectors:
  t1.vec = double(N)
  lambda0.vec = double(N)
  lambda1.vec = double(N)
  beta.vec = double(N)
  
  for (i in 1:N) {
    t1 = t1.RW.Metropolis.sample(t1, sigma2.t, lambda0, lambda1, coal.df)
    lambda0 = lambda0.sample(t1, t.start, beta, coal.df)
    lambda1 = lambda1.sample(t.end, t1, beta, coal.df)
    beta = beta.sample(lambda0, lambda1)
    
    # Store results:
    t1.vec[[i]]      = t1
    lambda0.vec[[i]] = lambda0
    lambda1.vec[[i]] = lambda1
    beta.vec[[i]]    = beta
  }
  
  ret.list = list(t1.vec, lambda0.vec, lambda1.vec, beta.vec)
  names(ret.list) = c("t1", "lambda0", "lambda1", "beta")
  return(ret.list)
}


main_A = function() {
  
  t1.test = t1.RW.sample(1900, 4.0)
  print(t1.test)

  normalize_time = F
  coal.df = boot::coal
  
  if (normalize_time){
    # Makes numerics nicer, but will affect 
    # the Poisson process rates.
    translate = min(coal.df)
    scale = max(coal.df) - min(coal.df)
    coal.df = (coal.df - translate)/scale
  }

  t.start = coal.df$date[[1]]
  t.end = coal.df$date[[length(coal.df$date)]]   
  
  num.accidents = cumulative.coal.disasters(1851.6325, coal.df)
  print(num.accidents)
  
  plot.coal.disasters(coal.df)
  
  ######### Run MCMC chain with Hybrid Gibbs sampler #########
  # Initial values:
  t1 = (t.end + t.start)/2.0
  lambda0 = 2.0
  lambda1 = 2.0
  beta = 1.0
  
  # Hyper-parameter:
  sigma2.t = 9.0

  N = 1.0e5L
  gibbs.mcmc.res = hybrid.Gibbs.MCMC.coal.sampling(N, t1, lambda0, lambda1, 
                                                   beta, sigma2.t, coal.df)
  t1.res = gibbs.mcmc.res$t1
  plot(seq_along(t1.res), t1.res, type="l")  
}

main_A()
