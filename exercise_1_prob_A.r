
################ Problem A 1: ################
exp.A.1.sample = function(n, lambda) {
  # Generates 'n' samples from the Exponential 
  # distribution f(x) = lambda*exp(-lambda*x), x>0.
  xs = -(1.0/lambda)*log(runif(n))
  return(xs)
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
       breaks=100,
       xlim=c(-5, 5), ylim=c(0, test.f.density(0.0)))
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
multivar.norm.sample = function(mu, Sigma) {
  # Perform a single sampling from the multivariate
  # normal distribution N(mu, Sigma). 
  # Expect: mu in R^{N}, Sigma in R^{NxN} > 0.
  n = length(mu)
  A = chol(Sigma)
  std_norm_xs = box.muller.stdnorm.sample(n)
  multivar_xs = A %*% std_norm_xs + mu
  return(multivar_xs)
}
