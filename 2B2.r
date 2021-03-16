library(forecast)
library(tidyverse)
library(ggplot2)
library(latex2exp)
source("exercise_1_prob_A.r")

yt <- read.table("Gaussiandata.txt") %>% 
  setNames("y_t") %>% 
  mutate(t = row_number())

# Plotting the time series ------------------------------------------------

ggplot(yt, aes(x=t,y=y_t)) + 
  geom_point() +
  ggtitle("Gaussian data") +
  ylab(TeX("$y_t$")) +
  theme_bw()


# Gibbs sampling ----------------------------------------------------------

# gibbs = function(N){
#   #Generate a gibbs sampling Markov Chain of length N
#   T=20
#   thetas = vector(length = N)
#   etas = matrix(nrow = N, ncol = T)
#   
#   
# }

