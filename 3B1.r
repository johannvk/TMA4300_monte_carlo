library(tidyverse)
#Reading in the data
bilirubin <- read.table("ex3-additional-files/bilirubin.txt",header=T)

#Splitting up the set wrt person
p1_bil <- filter(bilirubin, pers=="p1")
p2_bil <- filter(bilirubin, pers=="p2")
p3_bil <- filter(bilirubin, pers=="p3")

#Making a boxplot
boxplot(log(meas)~pers, data = bilirubin )
title("Logarithm of concentration for each individual")
# boxplot(log(meas)~pers, data = bilirubin, ylab = "Log of concentration",
# title("Logarithm of concentration for each individual") )

p <- ggplot(bilirubin, aes(x=pers, y=log(meas))) + 
  geom_boxplot() +
  labs(x = "Individual", y = "Log of concentration (mg/dL)") + 
  ylim(-2,0)+
  theme_bw()
#title = "Logarithm of concentration for each individual"
  

# F-test ------------------------------------------------------------------


linear_model <- lm(log(meas)~pers,data=bilirubin)
summary(linear_model)

Fval <- summary(linear_model)$fstatistic[1]

# Permutation function ----------------------------------------------------

permtest <- function(){
  
  #Reading in the data
  bilirubin <- read.table("ex3-additional-files/bilirubin.txt",header=T)
  meas <- bilirubin$meas
  pers <- bilirubin$pers
  
  #Permuting the data by sampling without replacement
  perm_meas <- sample(meas)
  perm_bilirubin <- data.frame(perm_meas,pers)
  
  #Fitting a linear model to the permuted data
  linear_model <- lm(log(perm_meas)~pers,data=perm_bilirubin)

  #Extracting and returning the F-statistic
  Fval <- summary(linear_model)$fstatistic[1]
  return(Fval)
}

Fs = vector(length = 999)
for (i in (1:999)) {
  Fs[i] <- permtest()
}

pval <- sum(Fs>Fval)/999
pval
