set.seed(10)
#########Functions from earlier##########
box.muller.stdnorm.sample = function(n){
  num_samples = ceiling(n/2)
  u1s = runif(num_samples)
  u2s = runif(num_samples)
  magnitudes = sqrt(-2*log(u1s))
  xs = c(magnitudes*cos(2*pi*u2s), magnitudes*sin(2*pi*u2s))
  return(xs[1:n])
}

################ Problem C 1: ################
mc.normal.integral = function(n,h){
  stdnorm_sample = box.muller.stdnorm.sample(n) #draws a std normal sample of size n
  vals = h(stdnorm_sample) #evaluates the function h at the sampled points 
  #plot(stdnorm_sample,vals)
  est = (1/n)*sum(vals) #Computes the MC estimate to be returned
  #plot((vals-est)^2)
  var_est = 1/(n-1)*sum((vals-est)^2)#Computes the sample variance
  return(c(est, var_est)) #returns a vector containing the MC estimate and the sample variance
}

#Indicator function for the set of real numbers x>4
indicator.g4 = function(x){return(x>4)}

##Estimation of Theta = Prob(X>4)
theta_vec = mc.normal.integral(100000,indicator.g4)
theta = theta_vec[1]
var_theta = theta_vec[2]
theta_interval = c(theta+10^(-5/2)*qnorm(0.025)*sqrt(var_theta),theta+10^(-5/2)*qnorm(0.975)*sqrt(var_theta))

#Test function for montecarlo
square = function(x){
  return(x^2)
}


################ Problem C 2: ################
importance.sampling.expnor = function(n,h){
  xs = sqrt(16-2*log(runif(n))) #Samples the proposal distribution
  gs = exp(8)*xs*exp(-1/2*xs^2) 
  fs = (2*pi)^(-1/2)*exp(-1/2*xs^2)
  ws = fs/gs #computes the weights using the ratio of the probability densities
  vals = h(xs)
  est = (1/n)*sum(vals*ws) #computes the estimate
  var_est = 1/(n-1)*sum((vals*ws-est)^2) #Computes an estimate for the variance of the MC estimate
  return(c(est,var_est))
}

theta2_vec = importance.sampling.expnor(100000,indicator.g4)
theta2 = theta2_vec[1]
var_theta2 = theta2_vec[2]

theta2_interval = c(theta2+10^(-5/2)*qnorm(0.025)*sqrt(var_theta2),theta2+10^(-5/2)*qnorm(0.975)*sqrt(var_theta2))
width2 = theta2_interval[2]-theta2_interval[1]


################ Problem C 3: ################
antithetical.importance = function(n,h){
  us = runif(n)
  xs = c(sqrt(16-2*log(us)),sqrt(16-2*log(1-us))) #Samples antithetically from the distribution g
  gs = exp(8)*xs*exp(-1/2*xs^2) 
  fs = (2*pi)^(-1/2)*exp(-1/2*xs^2)
  ws = fs/gs #computes the weights using the ratio of the probability densities
  vals = h(xs)
  est = (1/(2*n))*sum(vals*ws) #computes the estimate
  var_est = 1/(2*n-1)*sum((vals*ws-est)^2) #Computes an estimate for the variance of the MC estimate
  return(c(est,var_est))
}

theta3_vec = antithetical.importance(50000, indicator.g4)
theta3 = theta3_vec[1]
var_theta3 = theta3_vec[2]

theta3_interval = c(theta3+10^(-5/2)*qnorm(0.025)*sqrt(var_theta3),theta3+10^(-5/2)*qnorm(0.975)*sqrt(var_theta3))
