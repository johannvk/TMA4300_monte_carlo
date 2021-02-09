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
  est = (1/n)*sum(vals) #Computes the MC estimate to be returned
  
  #Computing the sample variance, divided by n to account for CLT:
  var_est = 1/(n-1)*sum((vals-est)^2)/n 
  return(c(est, var_est)) #returns a vector containing the MC estimate and the sample variance
}

#Indicator function for the set of real numbers x>4
indicator.g4 = function(x){return(x>4)}

##Estimation of Theta = Prob(X>4)
theta_vec = mc.normal.integral(100000,indicator.g4)
theta = theta_vec[1]
var_theta = theta_vec[2]

mc.normal.confint.95 = function(meanest, varest){
  return(c(meanest+qnorm(0.025)*sqrt(varest),meanest+qnorm(0.975)*sqrt(varest))) 
}
theta_interval = mc.normal.confint.95(theta,var_theta)

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
  
  #Computing a variance estimate, divided by n to account for CLT:
  var_est = 1/(n-1)*sum((vals*ws-est)^2)/n 
  return(c(est,var_est))
}

theta2_vec = importance.sampling.expnor(100000,indicator.g4)
theta2 = theta2_vec[1]

theta2_interval = mc.normal.confint.95(theta2_vec[1],theta2_vec[2])
width2 = theta2_interval[2]-theta2_interval[1]


################ Problem C 3: ################

antithetical.importance = function(n,h){
  us = runif(n)
  
  #Sampling antithetically from the distribution g:
  xs1 = sqrt(16-2*log(us))
  xs2 = sqrt(16-2*log(1-us))
  gs1 = exp(8)*xs1*exp(-1/2*xs1^2)
  gs2 = exp(8)*xs2*exp(-1/2*xs2^2)
  fs1 = (2*pi)^(-1/2)*exp(-1/2*xs1^2)
  fs2 = (2*pi)^(-1/2)*exp(-1/2*xs2^2)
  
  #computing the weights using the ratio of the probability densities
  ws1 = fs1/gs1
  ws2 = fs2/gs2
  vals1 = h(xs1)
  vals2 = h(xs2)
  
  #computes the estimates for 1 and 2
  est1 = (1/(n))*sum(vals1*ws1) 
  est2 = (1/(n))*sum(vals2*ws2) 
  
  #Computes an estimate for the variance, and covariance of the estimators
  var_est1 = 1/(n-1)*sum((vals1*ws1-est1)^2) 
  var_est2 = 1/(n-1)*sum((vals2*ws2-est2)^2) 
  cov_est = 1/(n-1)*sum((vals1*ws1-est1)*(vals2*ws2-est2)) 
  
  #Calculating the mean and variance of the total averaged estimator
  tot_est = 1/2*(est1+est2)
  tot_var = 1/4*var_est1+1/4*var_est2+1/2*cov_est
  
  return(c(tot_est,tot_var))
}

theta3_vec = antithetical.importance(50000, indicator.g4)
theta3 = theta3_vec[1]

theta3_interval = mc.normal.confint.95(theta3_vec[1],theta3_vec[2])

