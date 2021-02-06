# Optimization function:
library(stats)


rejection.sampling.posterior = function(n,h){
  y1 = 125
  y2=18
  y3=20
  y4=34
  u = runif() 
  
}

log.likelihood = function(theta) {
  # Assume theta in (0, 1).
  ys = c(125, 18 + 20, 34)
  log_probs = log(c(2 + theta, 1 - theta, theta))
  return(sum(ys*log_probs))
}

bayesian.rejection.sampling = function() {
  
  # Find the 'envelope' constant c, being the 
  # maximum likelihood for the given data: 
  # ys = c(y1, y2, y3, y4) = c(125, 18, 20, 34).
  opt.result = stats::optimize(log.likelihood, interval=c(0, 1), 
                                   maximum=T)

    # log(c):
  max.log.li = opt.result$objective
  
  
}
