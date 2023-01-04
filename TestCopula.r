source('/Users/tobiasbrammer/Library/Mobile Documents/com~apple~CloudDocs/Documents/Aarhus Uni/7. semester/FinancialEconometrics/Functions/DCC.r')

# Load the required packages
library(rugarch) # For GARCH modeling

# Set the seed for reproducibility
set.seed(123)

# Generate some synthetic data
nObs <- 1000
rho <- 0.5
sigma1 <- 0.1
sigma2 <- 0.2
mY <- matrix(NA, nrow = nObs, ncol = 2)
mY[1, ] <- c(0, 0)
for (t in 2:nObs) {
  mY[t, 1] <- rho * mY[t - 1, 2] + sqrt(sigma1) * rnorm(1)
  mY[t, 2] <- rho * mY[t - 1, 1] + sqrt(sigma2) * rnorm(1)
}

# Estimate the DCC model
DCCModel <- EstimateDCC(mY[, 1], mY[, 2])

# Print the estimated parameters
print(DCCModel$dLLK) # Log-likelihood
print(DCCModel$aR[, , nObs]) # Correlation at the final time step
print(DCCModel$dA) # dA parameter
print(DCCModel$dB) # dB parameter
