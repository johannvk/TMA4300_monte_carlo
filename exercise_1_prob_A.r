

g_density = function(x, alpha) {
  c = (1.0/alpha + 1.0/exp(1))^(-1)
  if (is.vector(x)){
    ys = numeric(length(x))
    ys[x < 1.0] = c*(x[x < 1.0])^(alpha - 1)
    ys[1.0 <= x] = c*exp(-x[1.0 <= x])
    return(ys)
  } else {
    y = NaN
    if (x < 1){
      y = c*(x)^(alpha - 1)
    } else if (1 <= x) {
      y = c*exp(-x)
    }
    return(y)  
  }
}

g_inv = function(u, alpha) {
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

g_sample = function(n, alpha) {
  return(g_inv(runif(n), alpha))
}

test_g_sample = function(alpha=0.7, breaks=200) {
  N = 100000
  test_g_density = function(x) g_density(x, alpha)
  hist(g_sample(N, alpha), breaks=breaks, probability=T, xlim=c(0.0, 3.0),
       main=paste("Comparison of g-density and sample histogram, alpha =",
                  alpha), 
       xlab=paste("x ~ g(x),", "N =", N, "samples"), ylab="g(x)")
  curve(test_g_density, from=0.0, to=3.0, n=2000, add=T, 
        ylab="Problem A: g-density")  
}

test_g_sample()
