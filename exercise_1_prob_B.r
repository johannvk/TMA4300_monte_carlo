# Load functions used for problem A:
source("exercise_1_prob_A.r")


gamma.density = function(xs, alpha, beta) {
  # Return the gamma(alpha, beta) density
  # evaluated at points in the vector 'xs'.
  # Working on the log-scale.
  log_f = alpha*log(beta) - log(gamma(alpha))
  log_f = log_f + (alpha - 1)*log(xs) - beta*xs
  return(exp(log_f))
}

gamma.B.1.sample = function(n, alpha) {
  # Generate 'n' samples from a gamma(x | alpha, beta=1) 
  # distribution, with 0 < alpha < 1, using rejection 
  # sampling from the distribtion g(x|alpha) in Prob. A.2.
  
  # Define the 'envelope' constant:
  kappa = (1.0/gamma(alpha))*(1.0/alpha + 1.0/exp(1))
  draws = double(n)
  tot_draws = 0

  while (tot_draws < n) {
    # Draw the expected number of required draws each time:
    N = ceiling((n - tot_draws)*kappa)
    g.A_draws = g.A.sample(N, alpha)
    envelope_proportions = gamma.density(g.A_draws, alpha, beta=1.0) / 
                            (kappa*g.A.density(g.A_draws, alpha)) 
    accepted_draws = g.A_draws[runif(N) < envelope_proportions]

    # Find out how many draws to add to the end of 'draws':
    new_draws = min(length(accepted_draws), n - tot_draws)
    draws[(tot_draws+1):(tot_draws+new_draws)] = accepted_draws[1:new_draws]
    tot_draws = tot_draws + new_draws
  }
  return(draws)
}

# xs = gamma.B.1.sample(100, 0.96)
