library(tidyverse)
library(latex2exp)

yt <- read.table("Gaussiandata.txt") %>% 
  setNames("y_t") %>% 
  mutate(t = row_number())

Q.matrix <- function(theta){
  #Constructs the sparse presision matrix for the gaussian field
  Q = diag(c(1,5,rep(6,16),5,1)*theta)
  diag(Q[-1,]) <- c(-2,rep(-4,17),-2)*theta
  diag(Q[,-1]) <- c(-2,rep(-4,17),-2)*theta
  diag(Q[3:20,]) <- rep(1,18)*theta
  diag(Q[,3:20]) <- rep(1,18)*theta
  return(Q)
}



posterior.marginal.hyperparameter.approx <- function(N){
  
  yt <- read.table("Gaussiandata.txt") %>% 
    as.matrix()
  T=20 #Number of observations
  
  #Defining a multinormal density function, using log scale
  log.multinorm.density <- function(x,mu,Sigma){
    return(-T/2*log(2*pi)-1/2*log(det(Sigma))-0.5*t(x-mu)%*%solve(Sigma)%*%(x-mu))
  }
  #Defining a function for the latent gaussian field, using log scale
  log.gf.density <- function(etas,Q,theta){
    return((T-2)/2*log(theta)-0.5*t(etas)%*%Q%*%(etas))
  }
  
  post_theta_approx_ker <- function(theta){
  #Internal functions to compute the unnormalized values
    
    #Gaussian random field
    Q <- Q.matrix(theta)
    I <- diag(T)
    log_gf <- log.gf.density(etas, Q, theta)
    
    #Full conditional for eta
    log_post_eta <- log.multinorm.density(x = etas, mu = (solve(Q+I)%*%yt),Sigma = solve(Q+I))
    
    #Gamma prior
    log_gamma_prior <- (-theta)
    
    #Computing the log version of the proportonality relations
    log_sum <- log_gf+log_gamma_prior-log_post_eta
    
    #Returning the exponentiated value
    return(exp(log_sum))
  }
  
  #Choosing etas to avoid singularities, Monte Carlo estimates from B2
  etas <- c(-8.96740602, -6.10471353, -4.28967123, -3.16718189, -2.17194503, -1.29870432,
            -0.55681615,-0.38825264, -0.25562693, -0.13717433, -0.08126779,
            0.16168831,  0.54126891,  1.18293313, 1.73114604,  2.53866305,
            4.00250518,  6.16271147,  8.90207036, 11.52111310)
  
  #Creating grid of theta values
  thetas <- seq(from=0, to = 6, length.out = N)

  #Computing unnormalized values in each gridpoint
  post_marg_theta_unnormalized <- apply(as.matrix(thetas),1,FUN = post_theta_approx_ker)
  
  #Computing normalizing constant, using numerical interpolation and integration
  normalizing_integral <- pracma::piecewise(thetas,post_marg_theta_unnormalized)$area
  
  #Normalizing
  post_marg_theta_normalized <- post_marg_theta_unnormalized/normalizing_integral
  
  #Returning normalized values
  return(post_marg_theta_normalized)
}

marginal_theta_posterior <- posterior.marginal.hyperparameter.approx(100)


# B3 Plot -----------------------------------------------------------------

b3.hyperparameter.plot =function(N){
  thetas <- seq(from=0, to = 6, length.out = N)
  marginal_theta_posterior <- posterior.marginal.hyperparameter.approx(N) %>% 
    data.frame(param = marginal_theta_posterior) %>% 
    mutate(theta = thetas)
  
  gibbthetas = gibbs(100000)$thetas
  gibbthetas = gibbthetas[(100000/2):100000]
  gibbthetas <- data.frame(gibbthetas, param2=gibbthetas)
  
  p <- ggplot(data=marginal_theta_posterior, aes(x=theta, y=param)) + 
    #geom_histogram(binwidth=0.05,aes(y=..density..)) +
    geom_line() +
    geom_density(data = gibbthetas, aes(param2) ,inherit.aes = F, colour ="#FF9999") +
    ggtitle(TeX("INLA estimate for $\\pi(\\theta | \\textbf{y})$")) +
    xlab(TeX("$\\theta$")) +
    ylab("") +
    theme_bw()
  return(p)
}


# Etai posteriors ---------------------------------------------------------



marg.post.smooth <- function(etai,N){
  
  #Index with which we are concerned
  i=10
  T=20
  
  #Constructing the Theta grid
  thetas <- seq(from=0, to = 6, length.out = N)
  
  #Generating the marginal theta posterior
  marginal_theta_posterior <- posterior.marginal.hyperparameter.approx(N)
  
  #Constructing matrix linear in theta that can be scaled later
  Q1 <- Q.matrix(1)
  I <- diag(T)
  
  #Loading y data
  yt <- read.table("Gaussiandata.txt") %>% 
    as.matrix()
  
  #Looping through and computing from the full conditional of eta
  full_conditional_i <- vector(length = N)
  for (j in 1:N) {
    
    Q <- Q1*thetas[j]
    mu <- (solve(Q+I)%*%yt)[i]
    sig <- solve(Q+I)[i,i]
    
    full_conditional_i[j] = (2*pi*sig)^(-0.5)*exp(-0.5/sig*(etai-mu)^2)
  }

  #Computing the INLA integral
  int_est <- sum(full_conditional_i*marginal_theta_posterior)*6/N
  
  return(int_est)
}


test.smooth <- function(N,M){
  
  eta_grid <- seq(from=-2, to = 2, length.out = M)
  
  vec <- sapply(eta_grid, marg.post.smooth,N=N)
  #print(pracma::piecewise(eta_grid,vec*eta_grid)$area)
  return(vec)
}

smooth_test <- test.smooth(1000,100)
plot(seq(from=-5, to = 5, length.out = 100),smooth_test)
abline(v = -.137174)


# Smooth plot -------------------------------------------------------------



smooth.plot <- function(N,M){
  
  eta_grid <- seq(from=-2, to = 2, length.out = M)
  
  smooth_test <- test.smooth(N,M) %>% 
    data.frame(param = smooth_test) %>% 
    mutate(eta = eta_grid)
  
  p <- ggplot(data=smooth_test, aes(x=eta, y=param)) + 
    geom_line() +
    geom_vline(xintercept = -.137174, linetype="dotted", colour="blue") +
    geom_vline(xintercept = -0.214, linetype="dotted",) +
    ggtitle(TeX("INLA estimate for $\\pi(\\eta_{10} | \\textbf{y})$")) +
    xlab(TeX("$\\eta_{10}$")) +
    ylab("") +
    theme_bw()
  return(p)
}

