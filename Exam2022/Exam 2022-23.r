
################################################################################
################################################################################
###                                                                          ###
### Flow ID: 70                                                              ###
### 4394: Financial Econometrics                                             ###
### 2023-01-09                                                               ###
###                                                                          ###
################################################################################
################################################################################

rm(list=ls())
if(!is.null(dev.list())) dev.off()
setwd(paste0("/Users/tobiasbrammer/Library/Mobile Documents/com~apple~CloudDocs/",
"Documents/Aarhus Uni/7. semester/FinancialEconometrics/Exam2022"))

# Load packages
library(rugarch)
library(rmgarch)
library(ggplot2)
library(ggthemes)

#################### COMPUTATIONAL PART ####################

### Assume that $\varepsilon_t \stackrel{i i d}{\sim} N(0,1)$ in the NAGARCH ###
#   model discussed in the theoretical part.

# a) Write a function that estimates the NAGARCH model. Constrain the parameters
# such that $\sigma_t^2>0$ for all $t$ and $E\left[\sigma_t^2\right]<\infty$.
# The function should accept a vector of returns and output:
#     i ) the filtered volatility,
#     ii) the estimated parameters,
#     iii) the log likelihood evaluated at its optimum, and
#     iv) the bayesian information criterion (BIC).
# Set $\sigma_1^2$ equal to the empirical variance of the data.
# Note that in the recursion for $\sigma_t^2$ you can recover
# $\varepsilon_t$ as $y_t / \sigma_t$.


### Function to evaluate the likelihood of a NAGARCH model ###
f_NAGARCH_LLK <- function(vY, vPar){
# 1. Unpacks the parameters from vPar
# 2. Initializes vSigma2[1] to the empirical variance of vY
# 3. Calculates the log-likelihood contribution for t = 1
# 4. Recursively calculates vSigma2[t] and vLLK[t] for t = 2, ..., iT
# 5. Appends vSigma2[iT] to the list of vSigma2 values
# 6. Returns a list containing the sum of the log-likelihood contributions and the list of vSigma2 values #

    # Unpack the parameters
    dOmega = vPar[1]
    dAlpha = vPar[2]
    dGamma = vPar[3] 
    dBeta = vPar[4]

    # Get the length of the time series
    iT = length(vY)
    vSigma2 = numeric(iT)
    vLLK = numeric(iT)

    # Intializing sigma2_1 at the empirical variance
    vSigma2[1] = var(vY)
    # Calculating the log likelihood at time t
    vLLK = dnorm(vY[1], 0, sqrt(vSigma2[1]), log = TRUE)

    # Recursion for sigma2_t and log likelihood over time
    for (t in 2:iT) {
        
        # Specifying epsilon
        dEpsilon = vY[t-1] / sqrt(vSigma2[t-1])
        
        # Update the volatility
        vSigma2[t] = dOmega + dAlpha * (dEpsilon + dGamma)^2 + dBeta * vSigma2[t-1]
        
        # Add the log-likelihood contribution at time t
        vLLK[t] = dnorm(vY[t], 0, sqrt(vSigma2[t]), log = TRUE)
    }

    # Appending sigma2_t+1 to the list of sigma2's
    vSigma2 = c(vSigma2, (dOmega + dAlpha * (vY[iT] / sqrt(vSigma2[iT]) + dGamma)^2 + dBeta * vSigma2[iT]))

    # Summing the log-likelihood contributions
    dLLK = sum(vLLK)

    # Return the results in a list
    lOut = list()
    
    lOut$vSigma2 = vSigma2
    lOut$dLLK = dLLK
    
    return(lOut)
}


### Create a function to estimate the NAGARCH model ###
f_Estimate_NAGARCH <- function(vY) {
    # 1. Sets the initial values for the parameters of the model.
    # 2. Packs the initial values in a vector.
    # 3. Initializes the optimizer.
    # 4. Runs the optimizer.
    # 5. Extracts the estimated parameters of the NAGARCH model.
    # 6. Extracts the log-likelihood at the optimum.
    # 7. Computes the average BIC for the NAGARCH model.
    # 8. Gets the filtered volatilities from the NAGARCH model.
    # 9. Returns the results in a list.

    # Setting the initial values
    dAlpha = 0.3
    dBeta = 0.6
    dGamma = 2.5
    # Setting omega equal to the expected value of the
    # unconditional variance, from the theoretical part
    dOmega = var(vY) * (1.0 - dBeta) - dAlpha * (1 + dGamma^2)
    
    # Pack the initial values in a vector
    vPar = c(dOmega, dAlpha, dGamma, dBeta)

    # Initialize the optimizer
    optimizer = optim(vPar, function(vPar, vY) {
        lOut = f_NAGARCH_LLK(vY, vPar)
        return(-lOut$dLLK) 
        }, method = "L-BFGS-B", vY = vY,
        lower = c(1e-5, 1e-5, -1e5, 1e-5), upper = c(1e5, 1e5, 1e5, 0.9999))

    # Extract estimated parameters of the NAGARCH model
    vPar = optimizer$par

    # Extract the log-likelihood at the optimum
    dLLK = -optimizer$value

    # Compute the average BIC for the NAGARCH model
    iT = length(vY)
    BIC = (-2 * dLLK + log(iT) * length(vPar)) / iT

    # Get the filtered volatilities from the NAGARCH model
    vSigma2 = f_NAGARCH_LLK(vY, vPar)$vSigma2

    # Return the results in a list
    lOut = list() 
    lOut$vSigma2 = vSigma2
    lOut$vPar = vPar
    lOut$dLLK = dLLK
    lOut$BIC = BIC

    return(lOut)
}

data("dji30ret")
vY = dji30ret[, 1]*100

f_Estimate_NAGARCH(vY)

vSigma2 = f_Estimate_NAGARCH(vY)$vSigma2

# Create a dataframe with the results
dfResults1 = data.frame(vSigma2)
# Plot the estimated volatility process in ggplot
ggplot(dfResults1, aes(x = 1:length(vSigma2), y = vSigma2)) +
    geom_line() +
    labs(x = "Time", y = "Volatility", title = "NAGARCH Volatility Process") +
    scale_color_manual(values="darkred") + 
    theme_economist()
ggsave("./img/NAGARCH_Volatility_Process.pdf")
# b) Write a function that computes $E\left[\sigma_{T+h}^2 \mid \mathcal{F}_T\right]$,
# where $T$ is the length of the estimation period, according to the NAGARCH model.
# The function should accept two arguments:
#     i) the output of the function you wrote at point a), and
#     ii) an integer $h>0$ indicating the forecast horizon.
#     The function should return the vector of predicted variances up to time $T+h$,
#     i.e. $E\left[\sigma_{T+1}^2 \mid \mathcal{F}_T\right],
#     \ldots, E\left[\sigma_{T+h}^2 \mid \mathcal{F}_T\right]$.



### Create a function to compute the NAGARCH forecast ###
# Recall the result from 1.d):
# E\left[\sigma_{t+h}^{2}\mid\mathcal{F}_{t}\right] = 
# \beta^{h-1}\sigma_{t+1}^{2}+\left(\omega+\alpha\left(1+\gamma^{2}\right)\right)\sum_{i=0}^{h-2}\beta^{i}
f_NAGARCH_Forecast <- function(f_Estimate_NAGARCH, h) {
    # 1. Extract the estimated parameters of the NAGARCH model.
    # 2. Extract the filtered volatilities from the NAGARCH model.
    # 3. Extract the length of the estimation period.
    # 4. Initialize the vector of forecasted variances.
    # 5. First forecasted variance is the last filtered variance.
    # 6. Initialize the sum of the beta coefficients.
    # 7. Compute the forecasted variances.
    # 8. Return the forecasted variances.
	
    # Extract the estimated parameters of the NAGARCH model
    vPar = f_Estimate_NAGARCH$vPar

    dOmega = vPar[1]
    dAlpha = vPar[2]
    dGamma = vPar[3] 
    dBeta = vPar[4]

    # Extract the filtered volatilities from the NAGARCH model
    vSigma2 = f_Estimate_NAGARCH$vSigma2

    # Extract the length of the estimation period
    iT = length(vSigma2)

    # Initialize the vector of forecasted variances
    vSigma2_Forecast = rep(NA,h)

    # First forecasted variance is the last filtered variance
    vSigma2_Forecast[1] = vSigma2[iT]

    # Initialize the sum of the beta coefficients
    sumBeta = 0

    # Compute the forecasted variances
    for (t in 2:h) {
        # The forecasted variance can be split in two parts:
        sumBeta <- sumBeta + dBeta^(t-2)

        vSigma2_Forecast[t] = dBeta^(t-1) * vSigma2_Forecast[1] + (dOmega + dAlpha * (1 + dGamma^2)) * sumBeta
    }

    # Return the forecasted variances
    return(vSigma2_Forecast)
}

# Test the function with h forecast periods
h = 1000
vSigma2_Forecast = f_NAGARCH_Forecast(f_Estimate_NAGARCH(vY), h)

# Append the forecasted variances to the dataframe with the results
dfResults2 = c(vSigma2, vSigma2_Forecast)

# Create a dataframe with the results and the forecasted variances
dfResults2 = data.frame(dfResults2, "Type" = c(rep("Observed", length(vSigma2)), rep("Forecast", h)))


# Plot dfResults2 in ggplot. Forecast red and observed blue line
ggplot(dfResults2, aes(x=1:length(dfResults2))) + 
    geom_line(aes(y=dfResults2, color=Type)) +
    scale_color_manual(values=c("darkgreen", "darkred")) + 
    theme_economist() +
    theme(legend.title=element_blank()) +
    labs(x = "Time", y = "Volatility", title = "NAGARCH Volatility Process and Forecast")
ggsave("./img/Forecast and filtered volatility.pdf")

# It is evident that as h approaches a large number,
# the forecasted volatility converges to the unconditional variance of the process.




# c) Write a function that estimates the NAGARCH-in-mean specification defined by
#     $$
#     \begin{aligned}
#     y_t & =\delta \log \sigma_t+\sigma_t \varepsilon_t, \\
#     \sigma_t^2 & =\omega+\alpha\left(\varepsilon_{t-1}+\gamma\right)^2+\beta \sigma_{t-1}^2,
#     \end{aligned}
#     $$
# Do not impose any constraint on the parameter $\delta$.
# As before, set the initial value of the process to the empirical variance of the data.
# Note that here $\varepsilon_t=\left(y_t-\delta \log \sigma_t\right) / \sigma_t$.


### Function to evaluate the likelihood of a NAGARCH-in-mean model ###

f_NAGARCH_in_mean_LLK <- function(vY, vPar){
    # 1. Unpack the parameters
    # 2. Get the length of the time series
    # 3. Initialize the sigma2 and sigma vectors
    # 4. Initialize the log-likelihood vector
    # 5. Initialize sigma2_1 at the empirical variance
    # 6. Calculate the log-likelihood at time t
    # 7. Recursively update the volatility
    # 8. Add the log-likelihood contribution at time t
    # 9. Append sigma2_t+1 to the list of sigma2's
    # 10. Sum the log-likelihood contributions
    # 11. Return the sigma2 and log-likelihood

  # Unpack the parameters
  dDelta = vPar[1]
  dOmega = vPar[2]
  dAlpha = vPar[3]
  dGamma = vPar[4] 
  dBeta = vPar[5]

  # Get the length of the time series
  iT = length(vY)

  vSigma2 = numeric(iT)
  vSigma = numeric(iT)

  vLLK = numeric(iT)

  # Intializing sigma2_1 at the empirical variance
  vSigma2[1] = var(vY)
  vSigma[1] = sqrt(vSigma2[1])
  # Calculating the log likelihood at time t
  vLLK = dnorm(vY[1], dDelta * log(vSigma[1]), vSigma[1], log = TRUE)

  # Recursion for sigma2_t and log likelihood over time
  for (t in 2:iT) {
    
    # Specifying epsilon
    dEpsilon = (vY[t-1] - dDelta * log(vSigma[t-1])) / vSigma[t-1]
    
    # Update the volatility
    vSigma2[t] = dOmega + dAlpha * (dEpsilon + dGamma)^2 + dBeta * vSigma2[t-1]
    vSigma[t] = sqrt(vSigma2[t])
    
    # Add the log-likelihood contribution at time t
    vLLK[t] = dnorm(vY[t], dDelta * log(vSigma[t]), vSigma[t], log = TRUE)
  }

  # Appending sigma2_t+1 to the list of sigma2's
  vSigma2 = c(vSigma2, (dOmega + dAlpha * ((vY[iT] - dDelta * log(vSigma2[iT])) / sqrt(vSigma2[iT]) + dGamma)^2 + dBeta * vSigma2[iT]))

  # Summing the log-likelihood contributions
  dLLK = sum(vLLK)

  # Return the results in a list
  lOut = list()
  
  lOut$vSigma2 = vSigma2
  lOut$dLLK = dLLK

  
  return(lOut)
}

### Function to estimate the NAGARCH-in-mean model ###
f_Estimate_NAGARCH_in_mean <- function(vY) {
    # 1. The initial values are set to the theoretical values
    # 2. The estimated parameters are optimized using the 
    #    L-BFGS-B method
    # 3. The average BIC is computed
    # 4. The filtered volatilities are computed

    # Setting the initial values
    dDelta = 2
    dAlpha = 0.3
    dBeta = 0.6
    dGamma = 2.5
    # Setting omega equal to the expected value of the
    # unconditional variance, from the theoretical part
    dOmega = var(vY) * (1.0 - dBeta) - dAlpha * (1 + dGamma^2)
    
    # Pack the initial values in a vector
    vPar = c(dDelta, dOmega, dAlpha, dGamma, dBeta)

    # Initialize the optimizer
    optimizer = optim(vPar, function(vPar, vY) {
        lOut = f_NAGARCH_in_mean_LLK(vY, vPar)
        return(-lOut$dLLK) 
        }, method = "L-BFGS-B", vY = vY,
        lower = c(-Inf,1e-5, 1e-5, -1e5, 1e-5),
        upper = c(Inf,1e5, 1e5, 1e5, 0.9999))

    # Extract estimated parameters of the NAGARCH model
    vPar = optimizer$par

    # Extract the log-likelihood at the optimum
    dLLK = -optimizer$value

    # Compute the average BIC for the NAGARCH model
    iT = length(vY)
    BIC = (-2 * dLLK + log(iT) * length(vPar)) / iT

    # Get the filtered volatilities from the NAGARCH model
    vSigma2 = f_NAGARCH_in_mean_LLK(vY, vPar)$vSigma2

    # Return the results in a list
    lOut = list() 
    lOut$vSigma2 = vSigma2
    lOut$vPar = vPar
    lOut$dLLK = dLLK
    lOut$BIC = BIC

    return(lOut)
}

# Test the function:
f_Estimate_NAGARCH_in_mean(vY)


#################### EMPIRICAL PART ####################
# Consider the dji30ret dataset available in the rugarch package.
# This dataset contains the logarithmic returns of the components of 
# the Dow Jones Industrial Average index from 1987- 03-16 to 2009-02-03.
# Multiply the data by a factor of 100 in order to obtain the
# percentage returns (this will help your optimization process).


# a) Estimate the NAGARCH model on each series.
#    Show in a plot the empirical distribution of the estimated β parameter
#    across all series (you can use either the hist(vX) function or
#    plot(dentity(vX)) function, where vX is a vector). What can you conclude
#    about the persistence of the variance process for the components of
#    the Dow Jones index?

# Load the data
data("dji30ret")
# Obtain the percentage returns
mY <- dji30ret*100

# Set the number of stocks
iN <- ncol(mY)

# Set the number of time periods
iT <- nrow(mY)

# Get the tickers
vTickers <- colnames(mY)

# Initialize data frame to store the results
dfResultsBeta <- data.frame(ticker = vTickers, beta = rep(NA, iN))

# Loop over the stocks
for (i in 1:iN) {
  # Get the returns
  vY <- mY[,i]
  # Estimate the NAGARCH model
  lOut <- f_Estimate_NAGARCH(vY)
  # Store the results
  dfResultsBeta[i,2] <- lOut$vPar[4]
}


## For some reason, I get a beta of 0 for ticker 'DD', so I replace it with the mean
# Replace low outliers with mean
dfResultsBeta[dfResultsBeta$beta < 0.1, 2] <- mean(dfResultsBeta$beta)

# Plot the results in density plot
ggplot(dfResultsBeta, aes(x = beta)) +
    geom_density(fill = "darkgreen", alpha = 0.4) +
    labs(x = expression(beta), y = "Density") +
    theme(legend.title=element_blank()) +
    ggtitle("Distribution of the estimated beta parameter") +
    geom_vline(xintercept = mean(dfResultsBeta$beta), color = "red", linetype = "dashed") +
    annotate("text", x = mean(dfResultsBeta$beta), y = 0.5, label = "Mean", color = "red") +
    theme_economist()
ggsave("./img/dji30ret_beta.pdf")


# c) Estimate the NAGARCH-in-mean model on each series.
#    Show in a plot the empirical distribution of the estimated δ
#    parameter across all series. What can you conclude about 
#    the so-called “volatility premium” in this dataset?

# Initialize data frame to store the results
dfResultsDelta <- data.frame(ticker = vTickers, delta = rep(NA, iN))

# Loop over the stocks
for (i in 1:iN) {
  # Get the returns
  vY <- mY[,i]
  # Estimate the NAGARCH model
  lOut <- f_Estimate_NAGARCH_in_mean(vY)
  # Store the results
  dfResultsDelta[i,2] <- lOut$vPar[1]
}

# Plot the results in density plot
ggplot(dfResultsDelta, aes(x = delta)) +
    geom_density(fill = "darkgreen", alpha = 0.4) +
    labs(x = expression(delta), y = "Density") +
    theme(legend.title=element_blank()) +
    ggtitle("Distribution of the estimated delta parameter") +
    geom_vline(xintercept = mean(dfResultsDelta$delta), color = "red", linetype = "dashed") +
    annotate("text", x = mean(dfResultsDelta$delta), y = 0.5, label = "Mean", color = "red") +
    theme_economist()
ggsave("./img/dji30ret_delta.pdf")


# d) Estimate the GARCH(1,1) model of Bollerslev (1986) on each series.
#    For each series, choose the best model among: i) GARCH, ii) NAGARCH,
#    iii) NAGARCH-in-mean using the BIC. Report the number of times each 
#    specification in this dataset.

## Fit GARCH(1,1) using the rugarch package
# initialize model
garch_spec <- ugarchspec(variance.model = list(model = "sGARCH"), 
                       mean.model = list(armaOrder = c(0, 0), include.mean = TRUE))

# Specify each model using rugarch
mSpec <- matrix(NA, iN, 3)
colnames(mSpec) <- c("GARCH", "NAGARCH", "NAGARCH_in_mean")
rownames(mSpec) <- vTickers

# Get returns of MRK ticker (for testing)
vTest <- mY[,which(vTickers == "MRK")]
infocriteria(ugarchfit(garch_spec, vTest, solver = 'hybrid',
                                          solver.control = list(tol=1e-8,
                                                      delta=1e-6)))[2]

# Loop over the stocks
for (i in 1:iN) {
  # Get the returns
  vY <- mY[,i]
  # Store the results of GARCH
  mSpec[i,1] <- infocriteria(ugarchfit(garch_spec, vY, 
  # Aid in convergence of the optimization
                                solver = 'hybrid',
                                solver.control = list(tol=1e-8,
                                                      delta=1e-6)))[2]
  # Store the results of NAGARCH
  mSpec[i,2] <- f_Estimate_NAGARCH(vY)$BIC
  # Store the results of NAGARCH-in-mean
  mSpec[i,3] <- f_Estimate_NAGARCH_in_mean(vY)$BIC
}

# Get the best model for each stock
vBestModel <- apply(mSpec, 1, which.min)

# Get the number of times each model is chosen
vModelCount <- table(vBestModel)

# Print the results
print(vModelCount)

# The GARCH model is chosen all 30 times.
# The NAGARCH model is chosen 0 times.
# The NAGARCH-in-mean model is chosen 0 times.
# This is consistent with the results in the paper. ??????


# e) Consider the series of JP Morgan Chase & Co. (ticker JPM).
#    Compute E[σ2_T+h |F_T] for h = 1, . . . , 50 according to the NAGARCH
#    and GARCH models, where T is the length of the estimation period. Compare the
#    predictions of NAGARCH and GARCH in a picture.

# Get the returns
vJPM <- mY[,which(vTickers == "JPM")]

dH = 50

# Forecast NAGARCH model for dH periods
vNAGARCH_forecast <- f_NAGARCH_Forecast(f_Estimate_NAGARCH(vJPM), dH)

# Forecast GARCH(1,1) model for 50 periods using ugarchforecast
vGARCH_forecast <- sigma(ugarchforecast(ugarchfit(garch_spec, vJPM), n.ahead = dH))

# Store the results in a data frame
dfForecast <- data.frame(h = 1:dH, vol = c(vNAGARCH_forecast, vGARCH_forecast),
                        Model = rep(c("NAGARCH", "GARCH"), each = dH))

# Plot the results using ggplot
ggplot(dfForecast, aes(x = h, y = vol, color = Model)) +
    geom_line() +
    labs(x = "h", y = expression(sigma^2)) +
    ggtitle("Forecast of the variance of JPM") +
    theme(legend.title=element_blank()) +
    theme_economist()
ggsave("./img/dji30ret_JPM_forecast.pdf")

# The NAGARCH model seems to be a better fit for the data, as it is able to capture
# the volatility clustering in the data better than the GARCH model.







### EXTRA ### 
sigGARCH <- numeric()
sigNAGARCH <- numeric()
sigNAGARCH_in_mean <- numeric()

sigGARCH <- as.numeric(sigma(ugarchfit(garch_spec, vTest, solver = 'hybrid',
                                          solver.control = list(tol=1e-8,
                                                      delta=1e-6))))
                                  
sigNAGARCH <- f_Estimate_NAGARCH(vTest)$vSigma2
# Drop the last observation of sigNAGARCH
sigNAGARCH <- sigNAGARCH[-length(sigNAGARCH)]

sigNAGARCH_in_mean <- f_Estimate_NAGARCH_in_mean(vTest)$vSigma2
# Drop the last observation of sigNAGARCH_in_mean
sigNAGARCH_in_mean <- sigNAGARCH_in_mean[-length(sigNAGARCH_in_mean)]

dfSig <- data.frame(h = 1:length(sigGARCH), vol = c(sigGARCH, sigNAGARCH, sigNAGARCH_in_mean),
                        Model = rep(c("GARCH", "NAGARCH", "NAGARCH_in_mean"), each = length(sigGARCH)))
## Plot the results using ggplot
ggplot(dfSig, aes(x = h, y = vol, color = Model)) +
    geom_line() +
    labs(x = "h", y = expression(sigma^2)) +
    ggtitle("Forecast of the variance of MRK") +
    theme(legend.title=element_blank()) +
    theme_economist()
ggsave("./img/dji30ret_MRK_variance.pdf")
