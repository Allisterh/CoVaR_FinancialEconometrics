
## Method of moments and Maximum Likelihood 
## estimators for the Independent Random Volatility 

#simulate from the model
SimIRM <- function(iT, dZeta, dSigma2_w) {
  
  vU = rnorm(iT)
  vW = rnorm(iT, mean = dZeta, sd = sqrt(dSigma2_w))
  
  vY = exp(vW/2)*vU
  
  return(vY)
  
}

set.seed(123)

iT = 500
dZeta = 1
dSigma2_w = 0.5

vY = SimIRM(iT, dZeta, dSigma2_w)

plot(vY, type = "l")

MM_Est <- function(vY) {
  
  dMu2 = mean(vY^2)
  dMu4 = mean(vY^4)
  
  dZeta_MM = log(dMu2^2*sqrt(3)/sqrt(dMu4))
  dSigma2_w_MM = log(dMu4/(3*dMu2^2))
  
  vTheta_MM = c("zeta" = dZeta_MM, 
                "sigma2_w" = dSigma2_w_MM)
  
  return(vTheta_MM)
}

MM_Est(vY)

ML_Est <- function(vY) {
  
  vStarting = MM_Est(vY)
  
  optimizer = optim(vStarting, fn = function(vTheta, vY) {
    
    dLLK = 0
    for (t in 1:length(vY)) {
      
      dLLK = dLLK + log(integrate(function(dW, vTheta, dY) {
        
        dZeta = vTheta[1]
        dSigma2_w = vTheta[2]
        
        
        1/(2*pi*sqrt(dSigma2_w)) * exp(
          -0.5*(dW + dY^2/exp(dW) + (dW - dZeta)^2/dSigma2_w)
        )
        
      }, lower = -10, upper = 10, vTheta = vTheta, dY = vY[t])$value)
      
    }
    
    # robustness check
    if (!is.finite(dLLK)) {
      dLLK = -1e5
    }
    
    return(-dLLK)
  }, vY = vY, method = "L-BFGS-B", lower = c(-10, 1e-4), upper = c(10, 10))
  
  vTheta_ML = optimizer$par
  
  return(vTheta_ML)
}

# ML estimates
ML_Est(vY)
# MM estimates
MM_Est(vY)

### Small Monte Carlo experiment

# number of replicates
iB = 100

# sample size
iT = 100

# true parameters
dZeta = 1
dSigma2_w = 0.5

# array where to store estimates
aEst = array(NA, dim = c(iB, 2, 2), dimnames = list(NULL, c("zeta", "sigma2_w"), c("ML", "MM")))

for (b in 1:iB) {
  
  # simulate from the model
  vY = SimIRM(iT, dZeta, dSigma2_w)
  
  # ML estimates
  aEst[b,,"ML"] = ML_Est(vY)
  # MM estimates
  aEst[b,,"MM"] = MM_Est(vY)
  
  cat(round(b/iB,2), "\n")
}

par(mfrow = c(1,2))

plot(density(aEst[,"zeta","ML"]), main = "Zeta", ylim = c(0, 2.5))
lines(density(aEst[,"zeta","MM"]), col = "red")
abline(v = dZeta)
legend("topleft", col = c("black", "red"), lty = 1, legend = c("ML", "MM"))

plot(density(aEst[,"sigma2_w","ML"]), ylim = c(0, 1.5), main = c("Sigma2_w"))
lines(density(aEst[,"sigma2_w","MM"]), col = "red")
abline(v = dSigma2_w)

# Bias 
colMeans(aEst[,,"ML"]) - c(dZeta, dSigma2_w)
colMeans(aEst[,,"MM"]) - c(dZeta, dSigma2_w)

# RMSE

colMeans(t(t(aEst[,,"ML"]) - c(dZeta, dSigma2_w))^2)
colMeans(t(t(aEst[,,"MM"]) - c(dZeta, dSigma2_w))^2)




