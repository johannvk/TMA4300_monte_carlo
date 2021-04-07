
read.data = function() {
  z = read.table("ex3-additional-files/z.txt", 
               header=F, col.names=c('z'))
  u = read.table("ex3-additional-files/u.txt", 
               header=F, col.names=c('u'))

  observed = list(unlist(z), unlist(u))
  names(observed) = c("z", "u")
  return (observed)
}

EM.step = function(lambdas, z, u, c0, c1) {
  # Update the estimates for the parameters
  # lam.0.t and lam.1.t:
  n = length(z)
  
  lam0 = lambdas[1]
  kappa0 = sum((1 - u)*(1/lam0 - z/(exp(lam0*z) - 1)))
  
  lam1 = lambdas[2]
  kappa1 = sum(u*(1/lam1 - z/(exp(lam1*z) - 1)))
  
  lambdas.next = n/c(c0 + kappa0, c1 + kappa1)
  return(lambdas.next)
}

l2.norm = function(x) sqrt(sum(x^2))

EM.optimization = function(lambdas, z, u, rtol=1.0e-5, atol=1.0e-3) {
  # lambdas = c(lam.0, lam.1)
  if (length(lambdas) != 2) {
    stop("length(lambdas) must equal 2.")
  }
  # Use a relative tolerance, forward error stopping criterion.
  # When:
  #  ||lambdas_(t+1) - lambdas_(t)||/||lambdas_(t) - lambdas_(t-1)|| > (1-eps),
  # we say we have converged.
  # Or, converged when:
  #   ||lambdas_(t+1) - lambdas_(t)|| < eps.
  
  c0 = sum(u*z)
  c1 = sum(z) - c0
  
  iterate = function(lambdas) EM.step(lambdas, z, u, c0, c1)
  
  lambdas.prev = lambdas
  lambdas.next = iterate(lambdas)
  
  # Add storage of the lambda values:
  storage = cbind(lambdas.prev, lambdas.next)
  
  while (l2.norm(lambdas.next - lambdas.prev) > atol) {
    lambdas.prev = lambdas.next
    lambdas.next = iterate(lambdas.next)
    
    storage = cbind(storage, lambdas.next)
  }
  
  return (list(lambdas=lambdas.next, iterates=storage))  
}

C1.main = function() {
  observed = read.data()
  z = observed$z; u = observed$u
  
  # Initial guess:
  init.lambdas = c(1.0, 1.0)
  ml.lambdas = EM.optimization(init.lambdas, z, u)
  
  lambdas = ml.lambdas$lambdas
  cat(paste("\nlam0: ", format(lambdas[1], digits=5), 
            " lam1: ", format(lambdas[2], digits=5), 
            "\n\n", sep="")
      )
  print(dim(ml.lambdas$iterates))
  
  par(mfrow = c(1, 2))
  plot(ml.lambdas$iterates[1, ], ylab="lambda 0", ylim=c(0.8, 4.0))
  abline(h=lambdas[1],  
         col = "red",            # Modify color
         lty = "dashed",         # Modify line type
         lwd = 1.5)   
  
  plot(ml.lambdas$iterates[2, ], ylab="lambda 1", ylim=c(0.8, 10.0))
  abline(h=lambdas[2],  
         col = "red",            # Modify color
         lty = "dashed",         # Modify line type
         lwd = 1.5)
  mtext("Convergence of Exp(lamda) parameters", outer=T,  cex=1.5, line=-2.5)
}

C1.main()