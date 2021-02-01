library(ggplot2)

# Load functions used for problem A:
source("exercise_1_prob_A.r")

################ Problem B 1: ################
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
    if (new_draws > 0){
    draws[(tot_draws+1):(tot_draws+new_draws)] = accepted_draws[1:new_draws]
    tot_draws = tot_draws + new_draws }
  }
  return(draws)
}


################ Problem B 2: ################

gamma.B.2.sample = function(n, alpha) {
  # Generate 'n' samples from a gamma(x | alpha, beta=1) 
  # distribution, with 1 < alpha, using the 'ratio of 
  # uniforms' method.
  # Think these are correct: Had forgotten the ^(1/2)...
  log_a = 0.5*(alpha-1)*(log(alpha-1) - 1)
  log_b_plus = 0.5*(alpha+1)*(log(alpha+1) - 1)
  
  log.f = function(x) (alpha-1)*log(x) - x
  
  C_f.criterion = function(log_x1_x2) { 
    log_x1 = log_x1_x2[1]; log_x2 = log_x1_x2[2]
    return(log_x1 <= 0.5*log.f(exp(log_x2 - log_x1))) }
  
  sample_transform = function(log_x1_x2) {
    log_x1 = log_x1_x2[1]; log_x2 = log_x1_x2[2]
    return(exp(log_x2 - log_x1))
  }
  # Sample 'N' pairs of (x1, x2) from [0, a]X[0, b_plus]
  # each pass. Only accept the pairs if they are in C_f.
  draws = double(n); tot_draws = 0; num_attempts = 0
  while (tot_draws < n) {
    # Draw the remaining number of samples needed each time:
    N = n - tot_draws; num_attempts = num_attempts + N
    log_x1s = log_a + log(runif(N))
    log_x2s = log_b_plus + log(runif(N))
    
    tup_list = Map(c, log_x1s, log_x2s)
    accepted_draws = Filter(C_f.criterion, tup_list)
    accepted_draws = unlist(Map(sample_transform, accepted_draws))
    
    # Find out how many draws to add to the end of 'draws':
    new_draws = min(length(accepted_draws), n - tot_draws)
    if (new_draws > 0) {
      draws[(tot_draws+1):(tot_draws+new_draws)] = accepted_draws[1:new_draws]
      tot_draws = tot_draws + new_draws
    }
  }
  result = list(draws, num_attempts)
  names(result) = c("samples", "num.attempts")
  return(result)
}


display.attempts.per.sample = function(N = 1000) {
  num_alphas = 200
  alphas = seq(1.1, 2000, length.out=num_alphas)
  attempts_vec = double(num_alphas)
  # Find how many attempts are required to get
  # N samples for alpha in (1, 2000].
  for (i in 1:length(alphas)) {
    sample_result = gamma.B.2.sample(N, alpha=alphas[i])
    # print(sample_result$num.attempts)
    attempts_vec[[i]] = sample_result$num.attempts/N
  }
  p = ggplot() + aes(alphas, attempts_vec) + geom_point() + 
      labs(y=paste("Number of proposed values per point sampled."), 
           x="Value of Alpha") + 
      ggtitle("Depence on alpha when sampling from Gamma(alpha, beta=1)")
  p
}


test.gamma.B.2.sample = function() {
  alpha_t = 2000.0
  N = 5000
  result = gamma.B.2.sample(N, alpha=alpha_t)
  print(paste("It took", result$num.attempts, "attempts to draw", 
              length(result$samples), "samples."))

  hist(result$samples, breaks=100, probability = T, xlab="x",
       main=paste("Histogram of", N, "samples from Gamma(", 
                  alpha_t, ", 1 ) with density overlay."))
  test_gamma = function(x) dgamma(x, shape=alpha_t, scale=1.0)
  curve(test_gamma, add=T)
}

# display.attempts.per.sample()
# test.gamma.B.2.sample()


################ Problem B 3: ################

gamma.B.3.sample = function(n, alpha) {
  # Sample from a Gamma(alpha, beta=1) distribution,
  # Fractional part taken care of by gamma.B.1.sample.
  int_alpha = floor(alpha)
  frac_alpha = alpha - int_alpha
  
  # Add up samples from distribution gamma(alpha=1, beta=1) = Exp(1).
  x = unlist(Map(function(x) sum(exp.A.1.sample(int_alpha, lambda=1)), 
                 1:n))
  if (0.005 < frac_alpha) { 
    # Practical lower limit for frac_alpha:
    # lim (alpha -> 0) Gamma(alpha, beta=1) = 0.
    x = x + gamma.B.1.sample(n, alpha=frac_alpha)
  }
  return(x)
}

# alpha_t = 8.5
# gamma_xs = gamma.B.3.sample(10000, alpha_t)
# mean(gamma_xs)
# hist(gamma_xs, breaks=100, probability=T)


################ Problem B 4: ################

gamma.B.4.sample = function(n, alpha, beta) {
  xs_1 = gamma.B.3.sample(n, alpha)
  return(xs_1/beta)
}

# N = 10000
# gamma.4.result = gamma.B.4.sample(N, alpha=5.02, beta=0.5)
# hist(gamma.4.result, breaks=100, probability = T)
# mean(gamma.4.result)


################ Problem B 5: ################

beta.B.5.sample = function(n, alpha, beta) {
  # Sample 'n' values from Gamma(alpha, 1)
  # and Gamma(beta, 1).
  xs = gamma.B.4.sample(n, alpha=alpha, beta=1)
  ys = gamma.B.4.sample(n, alpha=beta,  beta=1)
  return(xs / (xs + ys))
}


test.beta.B.5 = function() {
  alpha = 0.8; beta = 0.8
  beta_test.result = beta.B.5.sample(10000, alpha, beta)
  hist(beta_test.result, breaks=100, probability = T)
  est_mean = mean(beta_test.result)
  exact_mean = alpha/(alpha + beta)
  print(paste("Est.mean - exact.mean:", format(est_mean - exact_mean, 
                                               scientific = T)))
}
