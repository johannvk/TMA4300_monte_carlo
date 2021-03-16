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
  
  #computing sample variances
  var_est <- apply(etas[N/2:N,],2,var)
  
  #Computing effective sample sizes
  #eff_size <- apply(etas[N/2:N,],2,effectiveSize)
  print(effectiveSize(etas[N/2:N,]))
  
  return(est)
}

mc.postmean.etas <- mc.estimate.smooth(10000)


