library(forecast)
library(tidyverse)
library(ggplot2)
library(latex2exp)
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
  Q = diag(c(1,5,rep(6,16),5,1))
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
  etavec = c(1:(T/2),(T/2):1)*100 #Careful about initial conditions
  print(paste("etavec",etavec))
  
  #Load the gaussian data
  y <- read.table("Gaussiandata.txt") %>% 
    as.matrix()
  
  for (i in 1:N) {
    #Theta update:
    ratepar = 1+1/2*t(etavec) %*% Q.matrix(1) %*% etavec
    print(paste("ratepar",ratepar))
    theta <- gamma.B.4.sample(1,T/2,ratepar)
    thetas[i] = theta
    print(paste("Theta",theta))

    #Eta update:
    Q <- Q.matrix(theta)
    print(Q)
    I = diag(20)
    etavec = one.multivar.norm.sample(mu = (solve(Q+I)%*%y),Sigma = solve(Q+I))
    etas[i,] = etavec
  }
  return(thetas)
}

