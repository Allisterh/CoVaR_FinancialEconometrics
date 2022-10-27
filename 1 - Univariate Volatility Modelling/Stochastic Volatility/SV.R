rm(list=ls()) 
setwd("/Users/tobiasbrammer/Library/Mobile Documents/com~apple~CloudDocs/Documents/Aarhus Uni/7. semester/FinancialEconometrics/Problem Sets/Problem Set 1")
## install the packages
 install.packages('stochvol')
 library('stochvol')

 # Simulating a time series with n = 500 observations
  sim = svsim(500,mu=-10,phi=0.99,sigma=0.2)

# Running the Monte Carlo Markov Chain
  draws = svsample(sim$y,draws=200000,burnin=1000,
                   thinpara=100,thinlatent=100,
                   priormu=c(-10,1),
                   priorphi=c(20,1.2),
                   priorsigma=0.2)
  
  
  # Simulate a highly persistent SV process of length 500
  sim <- svsim(500, phi = 0.99, sigma = 0.1)
  
  print(sim)
  summary(sim)
  plot(sim)
  
  ## Simulate an SV process with leverage
  sim <- svsim(200, phi = 0.94, sigma = 0.15, rho = -0.6)
  
  print(sim)
  summary(sim)
  plot(sim)
  
  ## Simulate an SV process with conditionally heavy-tails
  sim <- svsim(250, phi = 0.91, sigma = 0.05, nu = 5)
  
  print(sim)
  summary(sim)
  plot(sim)