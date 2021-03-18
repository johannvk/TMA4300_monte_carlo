library(forecast)
library(tidyverse)
library(ggplot2)
library(latex2exp)
library(coda)
source("exercise_1_prob_A.r")
source("exercise_1_prob_B.r")

yt <- read.table("Gaussiandata.txt") %>% 
  setNames("y_t") %>% 
  mutate(t = row_number())

# Plotting the time series ------------------------------------------------

ggplot(yt, aes(x=t,y=y_t)) + 
  geom_point() +
  ggtitle("Gaussian data") +
  ylab(TeX("$y_t$")) +
  theme_bw()


# Setup for the parameters of the posteriors ------------------------------

Q.matrix <- function(theta){
  #Constructs the sparse presision matrix for the gaussian field
  Q = diag(c(1,5,rep(6,16),5,1)*theta)
  diag(Q[-1,]) <- c(-2,rep(-4,17),-2)*theta
  diag(Q[,-1]) <- c(-2,rep(-4,17),-2)*theta
  diag(Q[3:20,]) <- rep(1,18)*theta
  diag(Q[,3:20]) <- rep(1,18)*theta
  return(Q)
}


# Gibbs sampling ----------------------------------------------------------

gibbs = function(N){
  #Generate a gibbs sampling Markov Chain of length N
  T=20
  thetas = vector(length = N)
  etas = matrix(nrow = N, ncol = T)
  
  #initialize theta and eta
  theta = 1
  etavec = c(1:(T/2),(T/2):1) #Careful about initial conditions
  
  #Load the gaussian data
  y <- read.table("Gaussiandata.txt") %>% 
    as.matrix()
  
  for (i in 1:N) {
    #Theta update:
    ratepar = 1+1/2*t(etavec) %*% Q.matrix(1) %*% etavec
    #sampling from gamma distribution
    theta <- gamma.B.4.sample(1,T/2,ratepar)
    thetas[i] = theta

    #Eta update:
    Q <- Q.matrix(theta)
    I = diag(20)
    #sampling from multivariate normal
    etavec = one.multivar.norm.sample(mu = (solve(Q+I)%*%y),Sigma = solve(Q+I))
    etas[i,] = etavec
  }
  #returning values as a list
  return(list(thetas=thetas,etas=etas))
}

mc.estimate.smooth <- function(N){
  #Function for computing the posterior mean of the marginal etas
  
  etas = gibbs(N)$etas #Sampling the etas from our gibbs sampler
  
  #Computing column wise MC estimates, keeping half of the simulated chain
  est <- colMeans(etas[N/2:N,]) 
  
  #computing sample variances, in order to find conf int
  var_est <- apply(etas[N/2:N,],2,var)
  
  #Computing effective sample sizes, with a function from the coda package
  eff_size <- coda::effectiveSize(etas[N/2:N,])
  
  #Computing the limits of the confidence interval
  conf_int_lower = est -  1.96*sqrt(var_est)/sqrt(eff_size)
  conf_int_upper = est +  1.96*sqrt(var_est)/sqrt(eff_size)
  
  #Returns posterior mean and confidence interval
  return(list(est=est, conf_int_lower=conf_int_lower, conf_int_upper = conf_int_upper))
}

#mc.postmean.etas <- mc.estimate.smooth(10000)

eta.plot <- function(N){
  
  mc.eta <- data.frame(mc.estimate.smooth(N)) %>% 
    mutate(t = row_number())
  #print(mc.eta)
  yt <- read.table("Gaussiandata.txt") %>% 
    setNames("y_t") %>% 
    mutate(t = row_number())
  #print(yt)
  p <- ggplot(yt, aes(x=t,y=y_t)) + 
    geom_point() +
    geom_line(data=mc.eta, aes(x=t, y=est)) +
    geom_ribbon(data=mc.eta, aes(x=t, ymin=conf_int_lower,ymax=conf_int_upper),
                inherit.aes = FALSE, alpha=0.3) +
    ggtitle(TeX("Estimated smooth effect E$(\\eta_t|\\textbf{y})$, with 95% confidence interval")) +
    theme_bw()
  return(p)
}


# Hyperparameter plots ----------------------------------------------------

hyperparameter.plot =function(N){
  thetas = gibbs(N)$thetas
  thetas = thetas[(N/2):N]
  thetas <- data.frame(thetas, param=thetas)
  p <- ggplot(data=thetas, aes(param)) + 
    geom_histogram(binwidth=0.05,aes(y=..density..)) +
    geom_density() +
    ggtitle(TeX("Estimate for $\\pi(\\theta | \\textbf{y})$")) +
    xlab(TeX("$\\theta$")) +
    theme_bw()
  return(p)
}
 


