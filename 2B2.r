library(forecast)
library(tidyverse)

yt <- read.table("Gaussiandata.txt") %>% 
  as.matrix()

plot(yt)
