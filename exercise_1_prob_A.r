
################ Problem A 1: ################
exp.A.1.sample = function(n, lambda) {
  # Generates 'n' samples from the Exponential 
  # distribution f(x) = lambda*exp(-lambda*x), x>0.
  xs = -(1.0/lambda)*log(runif(n))
  return(xs)
}

test.Exp.A.1.sample = function(lambda=2.5, breaks=100) {
  N = 1.0e5L
  test.Exp.density = function(x) dexp(x, rate=lambda)
  hist(exp.A.1.sample(N, lambda), breaks=breaks, probability=T, 
       xlim=c(0.0, 2.5),
       main=paste("Comparison of Exp-density and sample histogram, lambda =",
                  lambda), 
       xlab=paste("x ~ Exp(x | lambda),", "N =", N, "samples."), ylab="f(x)")
  curve(test.Exp.density, add=T)
  exp_sample = exp.A.1.sample(N, lambda=0.5)
  est_mu = mean(exp_sample)
  est_var = var(exp_sample)
  delta_mu = 1.96*sqrt(est_var/N)
  print(paste("A 95% confidence interval for the mean of an X ~ Exp(1/2)",
              "distribution, with analytical answer E[X] = 2: [", 
              format(est_mu-delta_mu, digits=4), ",",
              format(est_mu + delta_mu, digits=4), "]."))
}


################ Problem A 2: ################
g.A.density = function(x, alpha) {
  c0 = (1.0/alpha + 1.0/exp(1))^(-1)
  if (is.vector(x)){
    ys = numeric(length(x))
    ys[x < 1.0] = c0*(x[x < 1.0])^(alpha - 1)
    ys[1.0 <= x] = c0*exp(-x[1.0 <= x])
    return(ys)
  } else {
    y = NaN
    if (x < 1){
      y = c0*(x)^(alpha - 1)
    } else if (1 <= x) {
      y = c0*exp(-x)
    }
    return(y)  
  }
}

g.A.inv = function(u, alpha) {
  # Assume 0 <= u < 1, 0 < alpha < 1.
  # Adapted to vectors u.
  threshold = 1.0/(1.0 + alpha/exp(1))
  
  if (is.vector(u)) { 
    xs = numeric(length(u))
    xs[u < threshold] = (u[u < threshold]/threshold)^(1.0/alpha)
    xs[threshold <= u] = 1 - log((1 + exp(1)/alpha)*(1 - u[threshold<=u]))
    return(xs)
  }
  else {
    # Scalar version:
    if (u < threshold){
    x = (u/threshold)^(1.0/alpha)
    } 
    else {
      x = 1 - log((1 + exp(1)/alpha)*(1 - u))
    }
    return(x)
  }
}

g.A.sample = function(n, alpha) {
  return(g.A.inv(runif(n), alpha))
}

test.g.A.sample = function(alpha=0.7, breaks=200) {
  N = 100000
  test.g.density = function(x) g.A.density(x, alpha)
  hist(g.A.sample(N, alpha), breaks=breaks, probability=T, xlim=c(0.0, 3.0),
       main=paste("Comparison of g-density and sample histogram, alpha =",
                  alpha), 
       xlab=paste("x ~ g(x),", "N =", N, "samples"), ylab="g(x)")
  curve(test.g.density, from=0.0, to=3.0, n=2000, add=T, 
        ylab="Problem A: g-density")  
}
# test.g.A.sample()


################ Problem A 3: ################
# f(x) = c0*exp(alpha*x)/(1 + exp(alpha*x))^2
# c0 = alpha
# F(x) = (1/alpha)*log_e(y/(1 - y))

f.A.density = function(xs, alpha) {
  # Assume 'xs' are of vector type.
  log_ys = log(alpha) + alpha*xs - 2*log(1 + exp(alpha*xs))
  return(exp(log_ys))
}

f.A.sample = function(n, alpha) {
  us = runif(n)
  xs = (1.0/alpha)*(log(us) - log(1 - us))
  return(xs)
}

test.f.A.sample = function(N = 100000,
                           alpha=2.0) {
  test.f.density = function(xs) f.A.density(xs, alpha)
  hist(f.A.sample(N, alpha), probability=T,
       breaks=100, main=paste("Histogram of", N, 
                              "samples from f( x | alpha=", alpha,")"),
       xlab="x", ylab="f (x | alpha )",
       ylim=c(0, test.f.density(0.0)), xlim=c(-5.0, 5.0))
  curve(test.f.density, add=T)
}
# Samle funksjonene vi vil kjøre til slutt i en "main"-fil?
# test.f.A.sample(alpha=1.5)


################ Problem A 4: ################
# Use Box-Muller method to generate samples
# from the standard normal N(0, 1) distribution.

box.muller.stdnorm.sample = function(n){
  num_samples = ceiling(n/2)
  u1s = runif(num_samples)
  u2s = runif(num_samples)
  magnitudes = sqrt(-2*log(u1s))
  xs = c(magnitudes*cos(2*pi*u2s), magnitudes*sin(2*pi*u2s))
  return(xs[1:n])
}


################ Problem A 5: ################
# Transform a standard-normal vector of samples
# to a multivariate normal N(mu, Sigma).

# Ganske streng: Lett for å lage matriser A * A^T som
# i prinsippet burde være pos.def. men som i praksis
# mister den egenskapen grunnet numeriske feil.
one.multivar.norm.sample = function(mu, Sigma, A=NULL) {
  # Perform a single sampling from the multivariate
  # normal distribution N(mu, Sigma). 
  # Expect: mu in R^{N}, Sigma in R^{NxN} > 0.
  n = length(mu)
  if (is.null(A)) { A = t(chol(Sigma)) }
  std_norm_xs = box.muller.stdnorm.sample(n)
  multivar_xs = A %*% std_norm_xs + mu
  return(multivar_xs)
}

multivar.norm.sample = function(n, mu, Sigma){
  # Sample 'n' times from a multivariate
  # normal distribution with mean vector
  # 'mu' and covariance matrix 'Sigma'.
  
  A = t(chol(Sigma))
  dim = length(mu)
  samples = matrix(rep(0, n*dim), nrow=n, ncol=dim)
  for(i in 1:n){
    samples[i, ] = one.multivar.norm.sample(mu, 0, A)
  }
  return(samples)
}


test.std.normal.sampling = function() {
  # Generate n=1.e5L samples from the standard normal
  # distribution N(0, 1), and compare the probability 
  # histogram with the probability density.
  N = 100000
  test.f.density = function(xs) dnorm(xs)

  norm_sample = box.muller.stdnorm.sample(N)
  hist(norm_sample, probability=T,
       breaks=100, main=paste("Histogram of", N, 
                              "samples from N( mu=0 | sigma=1 )"),
       xlab="x", ylab="f_N (x | mu=0, sigma=1 )",
       ylim=c(0, test.f.density(0.0)), xlim=c(-4.0, 4.0))
  curve(test.f.density, add=T)
  
  # Should be 1:
  est_mean = mean(norm_sample)
  est_var = var(norm_sample)
  
  est_skewness = (1/N)*sum((norm_sample - est_mean)^3)/est_var^(3/2)
  
  print(paste("Sample characteristics:", "Mean:", 
              format(est_mean, scientific=T), "Var:", 
              format(est_var, scientific=T),
              "Skewness:", format(est_skewness, scientific=T)))
  delta_mean = 1.96*sqrt(est_var/N)
  mean.conf.int = c(est_mean - delta_mean, est_mean + delta_mean)
  print(paste("Mean 95% conf.int.: (", mean.conf.int[1], ",", 
              mean.conf.int[2], ")."))
}


test.multivariate.normal.samling = function(){
  N = 1.0e5L
  mu = c(0, 2, -2)
  Sigma = matrix(c(7, -3, 2, -3, 6, 1, 2, 1, 4), nrow=3, ncol=3, byrow=T)
  
  samples = multivar.norm.sample(N, mu, Sigma)  
  
  est_mean = unlist(Map(function(i) mean(samples[, i]), c(1, 2, 3)))
  est_Sigma = cov(samples, samples)
  print("Estimated mean:")
  print(est_mean)
  print("Estimated Sigma:")
  print(est_Sigma)
  
  diff_Sigma_Inf_norm = norm(Sigma - est_Sigma, "I")
  print(paste("Infinity norm of (Sigma - Sigma_hat):", 
              format(diff_Sigma_Inf_norm, scientific=T)))
  rel_Sigma_error = norm(Sigma - est_Sigma, "2")/norm(Sigma, "2")
  print(paste("Relative covariance matrix error:", rel_Sigma_error))
}
