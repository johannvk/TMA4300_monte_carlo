set.seed(5)
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
