# Optimization function:
library(stats)

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
  print(env_const_log)
  
  thetas = vector(length = n)
  for (i in (1:n)){
    print(i)
    us = c(1,0)
    while(log(us[1])>log.posterior.kernel(us[2])-env_const_log){
      us = runif(2)
      #print(us)
    }
    thetas[i]=us[2]
  }
  return(thetas)
}
mc.mean.D= function(n){ #Posterior mean estimated by Monte Carlo integration
  thetas = rejection.sampling.D1(n)
  mc_mean = 1/n*sum(thetas)
  return(mc_mean)
}
