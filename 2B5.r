library("INLA")
library("tidyverse")

#reading in the data
yt <- read.table("Gaussiandata.txt") %>% 
  setNames("y") %>% 
  mutate(t = row_number()) %>% 
  data.frame()

result <- inla(
  #Specifying the linear predictor, with the correct latent gaussian field, and 
  #correct choice of hyperparameter
  formula = y ~ -1 + f(t, model = "rw2", hyper = list(prec = list(prior = "loggamma", param = c(1, 1)))
                        constr = FALSE),
  #NB the parameter constr can be chosen true or
  
  
  #Distibution family of the likelihood
  family = "gaussian",
  
  #Gaussian data
  data = yt,
)

#Extracting our wanted posterior estimates
eta_10 = result$marginals.random$t$index.10
thetas = result$marginals.hyperpar$`Precision for t`

# m.adjust = m
# m.adjust[,1] = m.adjust[,1]+0.65 


# Thetaplot ---------------------------------------------------------------


thetaplot.INLA <- function(){
  
  result <- inla(
    formula = y ~ -1 + f(t, model = "rw2",
          hyper = list(prec = list(prior = "loggamma", param = c(1, 1)))),
    family = "gaussian",
    data = yt
  )
  
  inlathetas <- result$marginals.hyperpar$`Precision for t` %>% 
    data.frame() %>% 
    filter(x<8)
  
  #lines(thetas)
  
  
  thetas <- seq(from=0, to = 6, length.out = 100)
  marginal_theta_posterior <- posterior.marginal.hyperparameter.approx(100) %>% 
    data.frame(param = marginal_theta_posterior) %>% 
    mutate(theta = thetas)
  
  gibbthetas = gibbs(100000)$thetas
  gibbthetas = gibbthetas[(100000/2):100000]
  gibbthetas <- data.frame(gibbthetas, param2=gibbthetas)
  
  p <- ggplot(data=marginal_theta_posterior, aes(x=theta, y=param)) + 

    geom_line() +
    geom_density(data = gibbthetas, aes(param2) ,inherit.aes = F, colour ="#FF9999") +
    geom_line(data = inlathetas, inherit.aes = F, aes(x=x, y=y), colour ="#0000FF") +
    ggtitle(TeX("INLA estimate for $\\pi(\\theta | \\textbf{y})$")) +
    xlab(TeX("$\\theta$")) +
    ylab("") +
    theme_bw()
  return(p)
  
}


# Etaplot -----------------------------------------------------------------



etaplot.INLA<- function(N,M){
  
  result <- inla(
    formula = y ~ -1 + f(t, model = "rw2",
                         hyper = list(prec = list(prior = "loggamma", param = c(1, 1))),
                         constr = FALSE),
    family = "gaussian",
    data = yt
  )
  inlaetas <- result$marginals.random$t$index.10 %>% 
    data.frame()
  
  
  eta_grid <- seq(from=-2, to = 2, length.out = M)
  
  smooth_test <- test.smooth(N,M) %>% 
    data.frame(param = smooth_test) %>% 
    mutate(eta = eta_grid)
  
  p <- ggplot(data=smooth_test, aes(x=eta, y=param)) + 
    geom_line(colour="blue") +
    geom_line(data = inlaetas, inherit.aes = F, aes(x=x, y=y), colour ="#000000")+
    geom_vline(xintercept = -.137174, linetype="dotted", colour="red") +
    #geom_vline(xintercept = -0.218, linetype="dotted",) +
    ggtitle(TeX("INLA estimate for $\\pi(\\eta_{10} | \\textbf{y})$")) +
    xlab(TeX("$\\eta_{10}$")) +
    ylab("") +
    xlim(c(-3,3))
    theme_bw()
  return(p)
  
}

et1 = result$marginals.random$t$index.10
print(pracma::piecewise(et1[,1],et1[,2])$area)
plot(et1[,1],et1[,2])

#hyper = list(prec = list(prior = "gamma", param = c(1, 1)))

# result <- inla(
#   formula = y ~ -1 + f(t, model = "rw2", hyper = list(prec = list(prior = "loggamma", param = c(1, 1)))),
#   
#   family = "gaussian",
#   
#   data = list(y = yt$y, t = yt$t),
#   
#   verbose = TRUE
# )

