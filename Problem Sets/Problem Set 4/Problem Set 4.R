rm(list=ls())
setwd("/Users/tobiasbrammer/Library/Mobile Documents/com~apple~CloudDocs/Documents/Aarhus Uni/7. semester/FinancialEconometrics/Problem Sets/Problem Set 4")

################################################################################
### Problem 1                                                                ###
################################################################################

### Point 1)
# Write a function to perform filtering of the volatility in the Stochastic
# Volatility (E[\exp\left(\alpha_{t}/2\right)\mid y_{1:t}) using the Bootstrap
# filter reported in slide 25 of Lecture 8.

bootstrap_filter <- function(returns, omega, phi, tau, N = 10000) {
    # function to perform bootstrap particle filter algorithm

    # number of observations
    obs <- length(returns)

    # matrix of bootstrapped particles
    alpha_bootstrap <- matrix(
        data = NA, nrow = obs, ncol = N
    )

    # vector of filtered volatility values
    volatility <- numeric(obs)

    # importance weight at each t
    weights <- numeric(N)


    # initialize values with draws from unconditional distribution
    # fill in first row
    alpha_bootstrap[1, ] <- rnorm(
        N,
        mean = omega / (1.0 - phi),
        sd = sqrt(tau**2 / (1.0 - phi**2))
    )
    # compute importance weights following notes in .pdf/.tex file
    # we draw `N` draws from gaussian distribution with varying scale
    weights <- dnorm(
        returns[1],
        mean = 0, sd = exp(alpha_bootstrap[1, ] / 2.0)
    )

    # normalize importance weights
    weights <- weights / sum(weights)

    # approximate \hat{\sigma}_{t}=E[\exp(\alpha_{t}/2)\mid y_{1:t}] for t = 1
    # this is a weighted average
    volatility[1] <- sum(exp(alpha_bootstrap[1, ] / 2) * weights)

    # draw `N` samples from `alpha_bootstrap` with replacement and
    # corresponding probabilities
    alpha_bootstrap[1, ] <- sample(
        alpha_bootstrap[1, ],
        size = N,
        replace = TRUE,
        prob = weights
    )

    for (t in seq(2, obs)) {
        # generate data
        alpha_bootstrap[t, ] <- omega + phi * alpha_bootstrap[t - 1, ] +
            rnorm(N, mean = 0, sd = tau)

        # draw weights
        weights <- dnorm(
            returns[t],
            mean = 0,
            sd = exp(alpha_bootstrap[t, ] / 2.0)
        )

        # normalize weights
        weights <- weights / sum(weights)

        # Approximate \hat{\sigma}_{t}=E[\exp(\alpha_{t}/2)\mid y_{1:t}] for t
        volatility[t] <- sum(exp(alpha_bootstrap[t, ] / 2) * weights)

        # resample
        alpha_bootstrap[t, ] <- sample(
            alpha_bootstrap[t, ],
            size = N,
            replace = TRUE,
            prob = weights
        )
    }


    return(
        list(
            volatility = volatility
        )
    )
}



### Point 2)
# Write a function like the one in point 1) but with resampling that occurs
# only when the Effective Sample Size (slide 26 of Lec. 8) is below the
# threshold gN. Let g be an argument of this function. Note that, when
# g = 1 resampling occurs at each iteration of the algorithm.

effective_sample_size <- function(weights, g, N) {
    #' returns boolean based on whether to perform resampling
    #' based on effective sample size.
    #' TRUE if ESS < gN
    #' FALSE if ESS >= gN
    #'
    #' weights:  `vector` of weights
    #' g:        `double` ESS threashold for effective sample size
    #' m:        `int` number og bootstraps

    ess <- 1.0 / sum(weights**2)

    return(
        ess < (g * N)
    )
}


bootstrap_filter <- function(returns, omega, phi, tau, ess_g, N = 10000) {
    #' #####################################################
    #' ### refactor of function in 1) to account for ESS ###
    #' #####################################################
    #'
    #' function to perform bootstrap particle filter algorithm
    #'
    #' returns:             `vector` of returns to model
    #' omega, phi, tau:     `double` parameters of SV model
    #' ess_g:               `double` threashold for effective sample size
    #' N = 10000:           `int` number of bootstraps
    #' 
    # number of observations
    obs <- length(returns)

    # matrix of bootstrapped particles
    alpha_bootstrap <- matrix(
        data = NA, nrow = obs, ncol = N
    )

    # vector of filtered volatility values
    volatility <- numeric(obs)

    # weights (tilde) for each bootstrap and time t
    weights <- matrix(NA, nrow = obs, ncol = N)

    # placeholder for normalized weights at each t
    norm_weights <- numeric(N)

    effective_sample_size <- function(weights, g, N) {
    #' returns boolean based on whether to perform resampling
    #' based on effective sample size.
    #' TRUE if ESS < gN
    #' FALSE if ESS >= gN
    #'
    #' weights:  `vector` of weights
    #' g:        `double` ESS threashold for effective sample size
    #' m:        `int` number og bootstraps

    ess <- 1.0 / sum(weights**2)

    return(
        ess < (g * N)
        )
    }


    # initialize values with draws from unconditional distribution
    # fill in first row
    alpha_bootstrap[1, ] <- rnorm(
        N,
        mean = omega / (1.0 - phi),
        sd = sqrt(tau**2 / (1.0 - phi**2))
    )
    # compute importance weights following notes in .pdf/.tex file
    # we draw `N` draws from gaussian distribution with varying scale
    weights[1, ] <- dnorm(
        returns[1],
        mean = 0,
        sd = exp(alpha_bootstrap[1, ] / 2.0)
    )

    # normalize importance weights
    norm_weights <- weights[1, ] / sum(weights[1, ])

    # approximate \hat{\sigma}_{t}=E[\exp(\alpha_{t}/2)\mid y_{1:t}] for t = 1
    # this is a weighted average
    volatility[1] <- sum(exp(alpha_bootstrap[1, ] / 2) * norm_weights)

    # draw `N` samples from `alpha_bootstrap` with replacement and
    # corresponding probabilities
    if (effective_sample_size(norm_weights, ess_g, N)) {
        alpha_bootstrap[1, ] <- sample(
            alpha_bootstrap[1, ],
            size = N,
            replace = TRUE,
            prob = norm_weights
        )
        # set weights equal to 1 for all
        weights[1, ] <- 1.0
    }

    for (t in seq(2, obs)) {
        # generate data
        alpha_bootstrap[t, ] <- omega + phi * alpha_bootstrap[t - 1, ] +
            rnorm(N, mean = 0, sd = tau)

        # draw weights
        weights[t, ] <- weights[t - 1, ] * dnorm(
            returns[t],
            mean = 0,
            sd = exp(alpha_bootstrap[t, ] / 2.0)
        )

        # normalize weights
        norm_weights <- weights[t, ] / sum(weights[t, ])

        # Approximate \hat{\sigma}_{t}=E[\exp(\alpha_{t}/2)\mid y_{1:t}] for t
        volatility[t] <- sum(exp(alpha_bootstrap[t, ] / 2) * norm_weights)

        # resample
        if (effective_sample_size(norm_weights, ess_g, N)) {
            alpha_bootstrap[t, ] <- sample(
                alpha_bootstrap[t, ],
                size = N,
                replace = TRUE,
                prob = norm_weights
            )

            # set weights equal to 1 for all
            weights[t, ] <- 1.0
        }
    }


    return(
        list(
            volatility = volatility
        )
    )
}


# verify my code with Leo's solution
# set.seed(69)
# my_boot <- bootstrap_filter(
#    simulated_sv$y,
#    omega = 0,
#    phi = 0.9,
#    tau = 0.5,
#    ess_g = 0.75,
#    N= 10
# )

# set.seed(69)
# leo_boot <- BootstrapFilter_ESS(
#    simulated_sv$y,
#    dOmega = 0,
#    dPhi = 0.9,
#    dTau = 0.5,
#    dG = 0.75,
#    iN = 10
# )


### Point 3)
# Simulate T = 1000 observations from the SV model reported in slide 12 of
# Lecture 8 with omega = 0, phi = 0.9, and tau = 0.5. Set the seed to 123.

sim_sv <- function(obs, omega, phi, tau) {
    #' function to simulate from SV model

    epsilon <- rnorm(obs, mean = 0, sd = 1)

    # placeholder for alpha
    alpha <- numeric(obs)

    # initialize alpha by drawing from unconditional distribution
    # derivation is available in .tex/.pdf file
    alpha[1] <- rnorm(
        1,
        mean = omega / (1.0 - phi),
        sd = sqrt(tau**2 / (1 - phi**2))
    )

    for (t in seq(2, obs)) {
        # draw eta on each iteration
        alpha[t] <- omega + phi * alpha[t - 1] + tau *
            rnorm(1, mean = 0, sd = 1)
    }

    y <- exp(alpha / 2) * epsilon

    return(
        list(
            "alpha" = alpha,
            "y" = y,
            "sigma2" = exp(alpha / 2) # calculate variance
        )
    )
}


# save the TRUE parameters in global scope
prior_omega <- 0.0
prior_phi <- 0.9
prior_tau <- 0.5


set.seed(123)
simulated_sv <- sim_sv(
    500,
    omega = prior_omega,
    phi = prior_phi,
    tau = prior_tau
)


### Point 4)
# Perform filtering of the volatility using the Bootstrap filter you
# derived in point 2) using N = 10000 particles and g = 1. Repeat also
# for g = 0.75 and g = 0.5. Is the quality of the estimate affected?
# Also play with N and see how the number of particles affect the
# precision of the estimate.


## Try with 10.000 particles
# number of particles
particles <- 10000

set.seed(123)
vol_bootstrapped_g1 <- bootstrap_filter(
    simulated_sv$y,
    omega = prior_omega,
    phi = prior_phi,
    tau = prior_tau,
    ess_g = 1,
    N = particles
)

set.seed(123)
vol_bootstrapped_g075 <- bootstrap_filter(
    simulated_sv$y,
    omega = prior_omega,
    phi = prior_phi,
    tau = prior_tau,
    ess_g = 0.75,
    N = particles
)

set.seed(123)
vol_bootstrapped_g050 <- bootstrap_filter(
    simulated_sv$y,
    omega = prior_omega,
    phi = prior_phi,
    tau = prior_tau,
    ess_g = 0.5,
    N = particles
)

# Plot the simulated volatility and the filtered volatility in one plot
if(!require(ggplot2)){install.packages('ggplot2')}
if(!require(ggthemes)){install.packages('ggthemes')}
# Create Dataframe for ggplot
df <- data.frame(
    time = 1:length(simulated_sv$sigma2),
    vol = c(simulated_sv$sigma2,
            vol_bootstrapped_g1$volatility,
            vol_bootstrapped_g075$volatility,
            vol_bootstrapped_g050$volatility),
    lab = rep(c("True","g = 1.0", "g = 0.75", "g = 0.50"),
            each = length(simulated_sv$sigma2))
    )

# Plot
ggplot(data = df, aes(x=time, y=vol, col = lab)) + 
                                    geom_line() +
                                    labs(title = paste0((particles), " Particles for Bootstrap")) +
                                    labs(x = "Time") +
                                    labs(y = "Volatility") +
                                    theme_economist() +
                                    theme(legend.title=element_blank())
ggsave("./img/10000 particles.pdf")


## Try with 25 particles
# number of particles
particles <- 25

set.seed(123)
vol_bootstrapped_g1 <- bootstrap_filter(
    simulated_sv$y,
    omega = prior_omega,
    phi = prior_phi,
    tau = prior_tau,
    ess_g = 1,
    N = particles
)

set.seed(123)
vol_bootstrapped_g075 <- bootstrap_filter(
    simulated_sv$y,
    omega = prior_omega,
    phi = prior_phi,
    tau = prior_tau,
    ess_g = 0.75,
    N = particles
)

set.seed(123)
vol_bootstrapped_g050 <- bootstrap_filter(
    simulated_sv$y,
    omega = prior_omega,
    phi = prior_phi,
    tau = prior_tau,
    ess_g = 0.5,
    N = particles
)

df <- data.frame(
    time = 1:length(simulated_sv$sigma2),
    vol = c(simulated_sv$sigma2,
            vol_bootstrapped_g1$volatility,
            vol_bootstrapped_g075$volatility,
            vol_bootstrapped_g050$volatility),
    lab = rep(c("True","g = 1.0", "g = 0.75", "g = 0.50"),
            each = length(simulated_sv$sigma2))
    )

# Plot
ggplot(data = df, aes(x=time, y=vol, col = lab)) + 
                                    geom_line() +
                                    labs(title = paste0((particles), " Particles for Bootstrap")) +
                                    labs(x = "Time") +
                                    labs(y = "Volatility") +
                                    theme_economist() +
                                    theme(legend.title=element_blank())
ggsave("./img/25 particles.pdf")

# The precision of the filtered volatility is affected by the number of particles.
# The more particles, the more precise the filtered volatility is.



### Point 5)
# Estimate the parameters omega, rho and tau using the QML estimator
# you derived in Exercise Set 3. Perform filtering via the Bootstrap
# filter using the estimated parameters. Is the precision of the
# filtered volatility affected?


# Source model from Problem Set 3

source("/Users/tobiasbrammer/Library/Mobile Documents/com~apple~CloudDocs/Documents/Aarhus Uni/7. semester/FinancialEconometrics/Functions/quasi_ml_sv.R")
source("/Users/tobiasbrammer/Library/Mobile Documents/com~apple~CloudDocs/Documents/Aarhus Uni/7. semester/FinancialEconometrics/Functions/Bootstrap.R")

# estimate model using QML
fit_qml <- quasi_ml_sv(simulated_sv$y)


# re-parameterize the models - following PS3 (B -> A)
omega_hat <- log(fit_qml$params["sigma"]**2) * (1.0 - fit_qml$params["rho"])
phi_hat <- fit_qml$params["rho"]
tau_hat <- sqrt(fit_qml$params["sigma2_eta"])

# re-parametrized parameters
omega_hat
phi_hat
tau_hat


# number of particles
particles <- 10000
# choose the "best" threshold for g
ess_threshold <- 0.75

set.seed(123)
vol_bootstrapped_true <- bootstrap_filter(
    simulated_sv$y,
    omega = prior_omega,
    phi = prior_phi,
    tau = prior_tau,
    ess_g = 0.75,
    N = particles
)

set.seed(123)
vol_bootstrapped_estimated <- bootstrap_filter(
    simulated_sv$y,
    omega = omega_hat,
    phi = phi_hat,
    tau = tau_hat,
    ess_g = 0.75,
    N = particles
)

# Create Dataframe for ggplot
df <- data.frame(
    time = 1:length(simulated_sv$sigma2),
    vol = c(simulated_sv$sigma2,
            vol_bootstrapped_true$volatility,
            vol_bootstrapped_estimated$volatility),
    lab = rep(c("True","True (Bootstrapped)", "Estimated (Bootstrapped)"),
            each = length(simulated_sv$sigma2))
    )

# Plot
ggplot(data = df, aes(x=time, y=vol, col = lab)) + 
                                    geom_line() +
                                    labs(title = paste0("True vs estimated vol")) +
                                    labs(x = "Time") +
                                    labs(y = "Volatility") +
                                    theme_economist() +
                                    theme(legend.title=element_blank())
ggsave("./img/true vs estimated vol.pdf")


# comments from solution
## Results are similar, and estimator error does not play a crucial role in our
## example. We are then satisfied with the results.
## Note that, this example shows us that we will never be able to recover
## the "true" volatility recursion.




################################################################################
### Problem 2                                                                ###
################################################################################
rm(list=ls())
source("/Users/tobiasbrammer/Library/Mobile Documents/com~apple~CloudDocs/Documents/Aarhus Uni/7. semester/FinancialEconometrics/Functions/Bootstrap.R")
source("/Users/tobiasbrammer/Library/Mobile Documents/com~apple~CloudDocs/Documents/Aarhus Uni/7. semester/FinancialEconometrics/Functions/quasi_ml_sv.R")

setwd("/Users/tobiasbrammer/Library/Mobile Documents/com~apple~CloudDocs/Documents/Aarhus Uni/7. semester/FinancialEconometrics/Problem Sets/Problem Set 4")

### Point 1)
# Download the time series of the S&P500 index from Yahoo finance from
# 2005-01-01 to 2018-01-01 and compute the percentage log returns. Replace the
# zero returns with their empirical mean.

library(quantmod)

price <- getSymbols(
    "^GSPC",
    from = "2005-01-01",
    to = "2018-01-01",
    auto.assign = FALSE
)

# extract and parse data
price <- price$GSPC.Adjusted
price <- as.numeric(price)

# auto-removes the NA observation
price <- diff(log(price)) * 100

# remove zeros
price[price == 0] <- mean(price)


### Point 2)
# Estimate the SV model by QML.

# estimate model using QML
fit_qml <- quasi_ml_sv(price)


# re-parameterize the models - following PS3 (B -> A)
omega_hat <- log(fit_qml$params["sigma"]**2) * (1.0 - fit_qml$params["rho"])
phi_hat <- fit_qml$params["rho"]
tau_hat <- sqrt(fit_qml$params["sigma2_eta"])

# re-parametrized parameters
omega_hat
phi_hat
tau_hat


### Point 3)

# Perform filtering using the Bootstrap filter with g = 1 and N = 10000.
sv_volatility <- bootstrap_filter(
    price,
    omega = omega_hat,
    phi = phi_hat,
    tau = tau_hat,
    ess_g = 1,
    N = 10000
)


### Point 4)
# Estimate a GARCH(1,1) with Gaussian shocks using the rugarch package.

library(rugarch)


spec = ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
                  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                  distribution.model = "norm")

fit_garch <- ugarchfit(spec, price)


### Point 5)
# Compare in a figure the series of filtered volatility from the SV model
# and the one obtained by the GARCH model.

# Visual inspection
df <- data.frame(
    time = 1:length(sv_volatility$volatility),
    vol = c(sv_volatility$volatility,
            as.numeric(sigma(fit_garch))),
    lab = rep(c("SV Model","GARCH(1,1)"),
            each = length(sv_volatility$volatility))
    )

# Plot
ggplot(data = df, aes(x=time, y=vol, col = lab)) + 
                                    geom_line() +
                                    labs(title = paste0("Volatility of SV vs GARCH(1,1)")) +
                                    labs(x = "Time") +
                                    labs(y = "Volatility") +
                                    theme_economist() +
                                    theme(legend.title=element_blank())
ggsave("./img/SV vs GARCH.pdf")
