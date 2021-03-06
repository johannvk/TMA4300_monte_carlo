library(boot)
library(ggplot2)

# Using our own functions from previous exercises:
source("exercise_1_prob_B.r")

cumulative.coal.disasters = function(t, coal_df){
  # Find the cumulative number of disasters before time t.
  # User decides whether or not to normalize the time of 
  # disasters to the interval [0, 1]

  # First and last time in the dataset are not disasters,
  # but the 'start'- and 'end'-time for the dataset.
  tot_num_disasters = length(coal_df$date) - 2

  if (t < coal_df$date[[2]]) {
    return(0)
  } 
  index_of_first_time_after_t = Position(function(x) x >= t, coal_df$date, 
                                         nomatch=NA)
  # Did not find a time in the dataset greater than 't':
  if (is.na(index_of_first_time_after_t)){
    return(tot_num_disasters)
  }
  # Return number of disasters before this time:  
  num_disasters = index_of_first_time_after_t - 2
  return(num_disasters)  
}

plot.coal.disasters = function(coal_df) {
  cumulative_disasters = mapply(function(t) 
                                cumulative.coal.disasters(t, coal_df), 
                                coal_df$date, SIMPLIFY=T)
  plot(coal_df$date, cumulative_disasters, type='l', col="blue",
       xlab="Time (Years)", 
       ylab="Cumulative Coal Disasters (10 or more deaths)", 
       main="Overview of Coal Mining Disasters in the UK (Jarrett, 1979)"
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

lambda0.sample = function(t1, t0, beta, coal_df) {
  # Sample from a Gamma distribution,
  # alpha.G = x(t1) - x(t0) + 2,
  # beta.G  = t1 - t0 + 1/beta
  alpha.G = cumulative.coal.disasters(t1) - cumulative.coal.disasters(t0) + 2
  beta.G = t1 - t0 + 1.0/beta
  lambda0.value = gamma.sample(1, alpha.G, beta.G)
  return(lambda0.value)
}

lambda1.sample = function(t2, t1, beta, coal_df) {
  # Sample from full conditonal for lambda1:
  # A Gamma distribution with parameters
  # alpha.G = x(t2) - x(t1) + 2,
  # beta.G  = t2 - t1 + 1/beta
  alpha.G = cumulative.coal.disasters(t2) - cumulative.coal.disasters(t1) + 2
  beta.G = t2 - t1 + 1.0/beta
  lambda1.value = gamma.sample(1, alpha.G, beta.G)
  return(lambda1.value)
}

main_A = function() {
  normalize_time = F
  coal_df = boot::coal
  
  if (normalize_time){
    # Makes numerics nicer, but will affect 
    # the Poisson process rates.
    translate = min(coal_df)
    scale = max(coal_df) - min(coal)
    coal_df = (coal_df - translate)/scale
  }

  t_start = coal_df$date[[1]]
  t_end = coal_df$date[[length(coal_df$date)]]   
  
  num_accidents = cumulative.coal.disasters(1851.6325, coal_df)
  print(num_accidents)
  
  plot.coal.disasters(coal_df)
}

main_A()
