#Reading in the data
bilirubin <- read.table("ex3-additional-files/bilirubin.txt",header=T)

linear_model <- lm(log(meas)~pers,data=bilirubin)
