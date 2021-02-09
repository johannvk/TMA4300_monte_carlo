# Optimization function:
library(stats)

################ Problem D 1: ################
#Function for calculating
posterior.kernel = function(theta) {
  # Assume theta in (0, 1).
  ys = c(125, 18 + 20, 34)
  probs = (2+theta)^ys[1]*(1-theta)^(ys[2])*theta^(ys[3])
  return(probs)
}

log.posterior.kernel = function(theta) {
  # Assume theta in (0, 1).
  ys = c(125, 18 + 20, 34)
  log_probs = log(c(2 + theta, 1 - theta, theta))
  return(sum(ys*log_probs))
}

rejection.sampling.D1 = function(n) {
  
  # Find the log of the enveloping constant c, being the 
  # maximum value of the log posterior 
  # ys = c(y1, y2, y3, y4) = c(125, 18, 20, 34).
  opt.result = stats::optimize(log.posterior.kernel, interval=c(0, 1), 
                                   maximum=T) 

   # log(c):
  env_const_log = opt.result$objective
  # print(env_const_log)
  
  thetas = vector(length = n)
  for (i in (1:n)){
    # print(i)
    us = c(1,0)
    while(log(us[1])>log.posterior.kernel(us[2])-env_const_log){
      us = runif(2)
      #print(us)
    }
    thetas[i]=us[2]
  }
  return(thetas)
}

################ Problem D 2: ################

mc.mean.D= function(n){ #Posterior mean estimated by Monte Carlo integration
  thetas = rejection.sampling.D1(n)
  mc_mean = 1/n*sum(thetas)
  return(mc_mean)
}

#Theoretical distribution function
posterior = function(theta){
  normalizer = integrate (posterior.kernel, 0,1)
  return(posterior.kernel(theta)/normalizer$value)
}

analytical.mean = function(){
  #numerical estimation of the theoretical mean
  
  #Integrand for use in numerical estimation of the theoretical mean
  integrand = function(theta){return(posterior(theta)*theta)}  
  
  num_mean = integrate(integrand,0,1)
  return(num_mean)
}

################ Problem D 3: ################

rejection.sampling.counter = function(n) {
  
  # Find the log of the enveloping constant c, being the 
  # maximum value of the log posterior 
  # ys = c(y1, y2, y3, y4) = c(125, 18, 20, 34).
  opt.result = stats::optimize(log.posterior.kernel, interval=c(0, 1), 
                               maximum=T) 
  # log(c):
  env_const_log = opt.result$objective
  #print(env_const_log)
  
  thetas = vector(length = n)
  attempts_needed = vector(length = n)
  for (i in (1:n)){
    # print(i)
    attempt=0
    us = c(1,0)
    while(log(us[1])>log.posterior.kernel(us[2])-env_const_log){
      us = runif(2)
      attempt=attempt+1
    }
    attempts_needed[i]=attempt
    thetas[i]=us[2]
  }
  print(paste("Average number of trials needed for acceptance:", 
              sum(attempts_needed)/n))
  return(thetas)
}

theoretical.acceptance = function(){
  
  #Calculates the log of the enveloping constant
  opt.result = stats::optimize(log.posterior.kernel, interval=c(0, 1), maximum=T) 
  
  # log(c):
  env_const_log = opt.result$objective
  
  integral = integrate(posterior.kernel,0,1)$value
  # print(integral)
  accept_prob = integral/exp(env_const_log)
  print(paste("Expected number of draws per accepted sample:", 
              format(1.0/accept_prob, digits=4)))
  return(accept_prob)
}
