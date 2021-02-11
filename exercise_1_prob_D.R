set.seed(5)
# Optimization function:
library(stats)
library(ggplot2)

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
  
  # Finding the log of the enveloping constant c, being the 
  # maximum value of the log posterior kernel: 
  opt.result = stats::optimize(log.posterior.kernel, interval=c(0, 1), 
                                   maximum=T) 

   # log(c):
  envelope_const_log = opt.result$objective
  
  #Prepares a vector to store the realizations
  thetas = vector(length = n)
  
  for (i in (1:n)){
    
    us = c(1,0) #Initialization of the while loop
    
    #The while loop tests the rejection condition
    while(log(us[1])>log.posterior.kernel(us[2])-envelope_const_log){
      #Generating the testing value u and the proposal value theta respectively:
      us = runif(2)
    }
    #Adding accepted theta to the storage vector:
    thetas[i]=us[2]
  }
  return(thetas)
}

################ Problem D 2: ################

#Posterior mean estimated by Monte Carlo integration
mc.mean.D= function(n){ 
  #Simulates n realisations of the posterior using the rejection algorithm    
  thetas = rejection.sampling.D1(n)
  
  #Computing the MC mean
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

verification.plot = function(n){
  fsample = rejection.sampling.D1(n)
  hist(fsample, probability=T,
       breaks=50, main=paste("Histogram of", n, 
                              "samples from", expression(f(θ|y))),
       xlab="θ", ylab="f(θ|y)",
       ylim=c(0, 9), xlim=c(0, 1))
  curve(posterior, add=T)
  #abline(v=0.6226, add=T)
  points(t(c(0.6226,-0.25)),pch="|",col="red")
}


################ Problem D 3: ################

rejection.sampling.counter = function(n) {
  
  # Finding the log of the enveloping constant c, being the 
  # maximum value of the log posterior 
  opt.result = stats::optimize(log.posterior.kernel, interval=c(0, 1), maximum=T) 
  
  # log(c):
  env_const_log = opt.result$objective

  #Preparing vectors for storage the realizations, and attempts needed for generation
  thetas = vector(length = n)
  attempts_needed = vector(length = n)
  
  for (i in (1:n)){
    
    attempt=0 #Variable for counting attempts
    
    us = c(1,0)#Initialization of the while loop
    
    #The while loop tests the rejection condition
    while(log(us[1])>log.posterior.kernel(us[2])-env_const_log){
      #Generating the testing value u and the proposal value theta respectively:
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
  
  #Calculating the log of the enveloping constant
  opt.result = stats::optimize(log.posterior.kernel, interval=c(0, 1), maximum=T) 
  envelope_const_log = opt.result$objective
  
  #Computes the integral numerically
  integral = integrate(posterior.kernel,0,1)$value
  
  #Calculates the acceptance probability
  accept_prob = integral/exp(envelope_const_log)
  
  #Prints the number of attempts until success, the reciprocal of the probability
  print(paste("Expected number of random pairs per accepted sample:", 
              format(1.0/accept_prob, digits=4)))
  return(accept_prob)
}

################ Problem D 4: ################

mc.importance.beta = function(n){
  #Sampling the proposal distribution
  thetas = rejection.sampling.D1(n)
  
  #Declaring and normalizing the new posterior distribution
  new.posterior.kernel = function(theta){return(posterior.kernel(theta)*(1-theta)^4)}
  normalising.constant = integrate(new.posterior.kernel,0,1)$value
  new.posterior = function(theta){return(new.posterior.kernel(theta)/normalising.constant)}
  
  #Computing the weights using the ratio of the probability densities:
  gs = posterior(thetas) 
  fs = new.posterior(thetas)
  ws = fs/gs 
  
  #Computing the estimate:
  est = (1/n)*sum(thetas*ws) 
  
  return(est)
}

beta5.analytical = function(){
  new.posterior.kernel = function(theta){return(posterior.kernel(theta)*(1-theta)^4)}
  normalising.constant = integrate(new.posterior.kernel,0,1)$value
  integrand = function(theta){return(theta*new.posterior.kernel(theta)/normalising.constant)}
  b5analyticmean = integrate(integrand,0,1)$value
  return(b5analyticmean)
}

