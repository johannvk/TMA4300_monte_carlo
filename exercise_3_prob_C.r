
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

EM.optimization = function(lambdas, z, u, rtol=1.0e-5, atol=1.0e-3, 
                           store.iterates=F) {
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
  
  l2.difference = l2.norm(lambdas.next - lambdas.prev)
  
  # Add storage of the lambda values:
  if (store.iterates){
    stored.lambdas = cbind(lambdas.prev, lambdas.next)
    stored.l2.diffs = c(l2.difference)
  } 
  
  while (l2.difference > atol) {
    lambdas.prev = lambdas.next
    lambdas.next = iterate(lambdas.next)
    l2.difference = l2.norm(lambdas.next - lambdas.prev)

    if (store.iterates) {
      stored.lambdas = cbind(stored.lambdas, lambdas.next)
      stored.l2.diffs = c(stored.l2.diffs, l2.difference)
    }
  }

  if (store.iterates){
    return (list(lambdas=lambdas.next, iterates=stored.lambdas, 
                 l2.diffs=stored.l2.diffs))  
  } else {
    return (lambdas.next)
  }
}


lambdas.bootstrap = function(B, z, u) {
  # Sample indices with resampling from the length
  # of z.
  n = length(z)
  indices = sample(1:n, B*n, replace=T)
  
  # Sample bootstrap samples in columns:
  z.boot = matrix(z[indices], nrow=n, ncol=B)
  u.boot = matrix(u[indices], nrow=n, ncol=B)
  
  # Estimate the lambdas (lam0, lam1) B times: 
  # Store them in one column each.
  lambdas.boot = matrix(NA, nrow=B, ncol=2)

  init_lambdas = c(1.0, 1.0)
  for (i in 1:B) {
    lambdas_i = EM.optimization(init_lambdas, z.boot[ , i], u.boot[ , i])
    lambdas.boot[i, ] = lambdas_i
  }
  
  return (lambdas.boot)
}


lambdas.log.likelihood = function(lambdas, z, c0) {
  # c0 = sum(u).
  n = length(z)
  lam0 = lambdas[1]; lam1 = lambdas[2]
  log.lam0 = log(lam0); log.lam1 = log(lam1)

  log.u = -n*(log(lam0 + lam1) - log.lam0) +
          (log.lam1 - log.lam0)*c0
  
  exp.lam0 = exp(-lam0*z); exp.lam1 = exp(-lam1*z)
  log.z = sum(
    log(
          lam0*exp.lam0*(1 - exp.lam1) +
          lam1*exp.lam1*(1 - exp.lam0)
        )
    )
  
  log.likelihood = log.u + log.z
  return (log.likelihood)
}


C1.main = function() {
  cat(paste("\nProblem C.1:\n"))
  observed = read.data()
  z = observed$z; u = observed$u
  
  # Initial guess:
  init.lambdas = c(1.0, 1.0)
  ml.lambdas = EM.optimization(init.lambdas, z, u, store.iterates = T)
  
  lambdas = ml.lambdas$lambdas
  cat(paste("The EM-algorithm converged on the values:",
    "\nlambda0: ", format(lambdas[1], digits=5), 
    ", lambda1: ", format(lambdas[2], digits=5), 
    "\n\n", sep="")
      )
  num.iterates = length(ml.lambdas$iterates[1, ])
  
  # Plot L2-norm difference per iteration, on log-scale:
  par(mfrow = c(1, 1))
  l2.diffs = ml.lambdas$l2.diffs
  l2.conv.df = data.frame(x=1:(num.iterates-1), y=l2.diffs)
  l2.conv.lm = lm(log(y, base=exp(1)) ~ x, data=l2.conv.df)
  conv.slope = coef(l2.conv.lm)[["x"]]

  plot(l2.conv.df$x, log(l2.diffs, base=exp(1)), 
       main="EM Convergence in L2-Norm", xlab="Iterations",
       ylab="ln(|| lambdas^(t+1) - lambdas^(t) ||_2)")
  abline(l2.conv.lm, col="steelblue")
  text(x = 8, y = -1,                # Text with different color & size
       paste("Slope:", format(conv.slope, digits=4)),
       col = "#1b98e0"
       )

  # Convergence of each parameter:
  par(mfrow = c(1, 2))
  plot(ml.lambdas$iterates[1, ], ylab="lambda 0",  # ylim=c(0.8, 4.0), 
       xlab="Iterations", xlim=c(1, num.iterates))
  abline(h=lambdas[1],  
         col = "red",            # Modify color
         lty = "dashed",         # Modify line type
         lwd = 1.5)   
  
  plot(ml.lambdas$iterates[2, ], ylab="lambda 1", # ylim=c(0.8, 10.0), 
       xlab="Iterations", xlim=c(1, num.iterates))
  abline(h=lambdas[2],  
         col = "red",            # Modify color
         lty = "dashed",         # Modify line type
         lwd = 1.5)
  mtext("Convergence of Exp(lamda) parameters", outer=T,  cex=1.5, line=-2.5)
  par(mfrow = c(1, 1))
}


C2.main = function() {
  cat(paste("\nProblem C.2:\n"))
  
  observed = read.data()
  z = observed$z; u = observed$u
  lambdas.full = EM.optimization(c(1.0, 1.0), z, u)
  
  B = 1.0e4L
  lambdas.boot = lambdas.bootstrap(B, z, u)
  lam0.boot = mean(lambdas.boot[ , 1])
  lam1.boot = mean(lambdas.boot[ , 2])
  
  cat(paste("Bootstrap Mean:\nlambda0: ", 
            format(lam0.boot, digits=5), "\tlambda1: ", 
            format(lam1.boot, digits=5), "\n", sep=""))
  
  bias0 = lam0.boot - lambdas.full[1]
  bias1 = lam1.boot - lambdas.full[2]

  par(mfrow = c(1, 2))
  hist(lambdas.boot[ , 1], probability = T, main="Hist. of lambda0",
       xlab="Bootstrapped lambda0 values.", breaks=50)
  abline(v=lambdas.full[1],  
         col = "red",            # Modify color
         # lty = "dashed",         # Modify line type
         lwd = 2.0)
  abline(v=lam0.boot,  
         col = "blue",          # Modify color
         # lty = "dashed",         # Modify line type
         lwd = 2.0)
  cat(paste("The estimated bias of the lambda0 estimate is: ",
            format(bias0, digits=6), "\n", sep=""))
  
  hist(lambdas.boot[ , 2], probability = T, main="Hist. of lambda1",
       xlab="Bootstrapped lambda1 values.", breaks=50)
  abline(v=lambdas.full[2],  
         col = "red",            # Modify color
         # lty = "dashed",         # Modify line type
         lwd = 2.0)
  abline(v=lam1.boot,  
         col = "blue",          # Modify color
         # lty = "dashed",         # Modify line type
         lwd = 2.0)
  cat(paste("The estimated bias of the lambda0 estimate is: ",
            format(bias1, digits=6), "\n", sep=""))
  par(mfrow = c(1, 1))
  
  lambda.estimates.correlation = cor(lambdas.boot[ , 1], lambdas.boot[ , 2])
  lambdas.est.covar = cov(lambdas.boot)

  cat(paste("\nThe correlation between lambda0 and lambda1 from ", B, 
            " bootstrap samples, is:\n\t", 
            format(lambda.estimates.correlation, digits=6), "\n", sep=""))
  cat(paste("The bootstrapp-estimated covariance matrix for",
            " (lambda0, lambda1):\n", sep=""))
  print(lambdas.est.covar)
}

C3.main = function() {
  cat(paste("\nProblem C.3:\n"))
  
  observed = read.data()
  z = observed$z; u = observed$u
  c0 = sum(u)
  
  objective = function(lambdas) {
      lambdas.log.likelihood(lambdas, z, c0)
  }
  init.lambdas = c(2.0, 6.0)
  control = list(fnscale=-1.0)
  optim.res = optim(init.lambdas, objective, method="L-BFGS-B",
                    lower=c(0.01, 0.01), upper=c(100.0, 100.0), 
                    control=control,
                    hessian=T)
  # str(optim.res)
  ml.lambdas = optim.res$par
  lambdas.est.covar = solve(-optim.res$hessian)
  
  cat(paste("The maximum likelihood estimates for (lambda0, lambda1)",
            " are:\n( ", format(ml.lambdas[1], digits=5), 
            ", ", format(ml.lambdas[2],digits=5),
            " )\n", sep=""))
  cat(paste("The 'maximum likelihood'-estimated covariance matrix for",
            " (lambda0, lambda1):\n", sep=""))
  print(lambdas.est.covar)
}

C1.main()
C2.main()
C3.main()