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

# Fitting a linear model
linear_model <- lm(log(meas)~pers,data=bilirubin)
summary(linear_model)

#Extracting the F-statistic
Fval <- summary(linear_model)$fstatistic[1]

#Computing the p-value
pvalue <- 1 - pf(Fval,2,26)   

# Permutation function ----------------------------------------------------

permTest <- function(){
  
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

Fs = vector(length = 1000)
Fs[1] = Fval
for (i in (2:1000)) {
  Fs[i] <- permTest()
}

pval <- sum(Fs>=Fval)/1000
pval

p2 <- ggplot(data=data_frame(Fs), aes(Fs)) + 
  geom_histogram(data=subset(data_frame(Fs),Fs<= Fval), bins = 30)+
  geom_histogram(data=subset(data_frame(Fs),Fs>Fval),fill ="red", bins = 30) +
  xlim(0,7.5) +
  labs(x = "F statistic", y = "Count") + 
  theme_bw()

# p2 <- ggplot(data=subset(data_frame(Fs),, aes(Fs)) + 
#                geom_histogram(bins = 20)+
#                xlim(0,7.5) +
#                geom_vline(xintercept = Fval, color="blue") +
#                theme_bw()
#              
# p2 <- ggplot(data=data_frame(Fs), aes(Fs)) + 
#   geom_density()+
#   theme_bw()


