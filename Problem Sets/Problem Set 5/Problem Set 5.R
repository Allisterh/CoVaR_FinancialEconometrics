library(Rsolnp)


################################################################################
### Problem 2                                                                ###
################################################################################

### Point 1)
# Create a function to evaluate the likelihood of the GAS model of the
# previous exercise.



tgas_log_likelihood <- function(y, params) {
    #' returns the log-likelihood of a t-distributed GAS model
    #'

    # unpack vector of parameters
    omega <- params[1]
    alpha <- params[2]
    beta <- params[3]
    nu <- params[4]

    # number of observations
    obs <- length(y)

    # placeholders for generated variables
    phi_tilde <- numeric(obs)
    phi <- numeric(obs)


    # initialize parameters
    phi_tilde[1] <- omega / (1.0 - beta)
    phi[1] <- exp(phi_tilde[1])

    # derive initial log-likelihood value
    # this is from the logged density from theoretical notes
    # ll <- dt(x = y[1] / phi[1], df = nu, log = TRUE) - log(phi[1])

    # alternative definition using the analytical expression from problem 1)
    ll <- log(gamma((nu + 1) / 2)) - log(sqrt(pi * nu) * gamma(nu / 2)) -
        log(phi[1]) - ((nu + 1) / 2) * log(1 + (y[1]**2 / (nu * phi[1]**2)))

    for (t in seq(2, obs)) {
        # z as helper variable - ratio of y and phi at time (t-1)
        z <- y[t - 1] / phi[t - 1]

        # calculate next value for phi
        phi_tilde[t] <- omega + alpha *
            ((((nu + 1.0) * z**2) / (nu + z**2)) - 1) +
            beta * phi_tilde[t - 1]

        # calculate phi from exponential mapping
        phi[t] <- exp(phi_tilde[t])

        # calculate log-likelihood
        # ll <- ll + dt(x = y[t] / phi[t], df = nu, log = TRUE) - log(phi[t])

        # again using the analytical solution
        ll <- ll + log(gamma((nu + 1) / 2)) - log(sqrt(pi * nu) *
            gamma(nu / 2)) - log(phi[t]) - ((nu + 1) / 2) *
            log(1 + (y[t]**2 / (nu * phi[t]**2)))
    }

    return(
        list(
            "neg_ll" = -ll,
            "phi" = phi
        )
    )
}


estimate_tgas <- function(returns) {
    #' TODO: add docstring

    require(Rsolnp)

    # set initial parameters for use in `params`
    initial_nu <- 5
    initial_beta <- 0.9

    # set starting parameters
    params <- c(
        omega = log(var(returns) * ((initial_nu - 2) / initial_nu)) *
            (1 / 2) * (1 - initial_beta),
        alpha = 0.05,
        beta = initial_beta,
        nu = initial_nu
    )

    # find optimal parameters using solnp solver
    # solnp allows for setting constraints
    optimizer <- solnp(
        pars = params,
        fun = function(...) {
            return(tgas_log_likelihood(...)$neg_ll)
        },

        # parameter constraints
        LB = c(-3.0, 0.0001, 0.00001, 2.01),
        UB = c(3.0, 2.0, 0.9999, 50),

        # returns vector - argument to neg ll function
        y = returns
    )

    # extract optimal parameters
    optimal_params <- optimizer$pars

    # estimate phi, sigma and log-likelihood for optimal parameters
    optimal_fit <- tgas_log_likelihood(returns, optimal_params)

    # deriving optimal volatility process
    sigma <- optimal_fit$phi * optimal_params["nu"] / (optimal_params["nu"] - 2)

    # return optimzer object with results
    return(
        list(
            "optimal_params" = optimal_params,
            "sigma" = sigma,
            "phi" = optimal_fit$phi
        )
    )
}


################################################################################
### Problem 2                                                                ###
################################################################################


### Point 1)
# Download the time series of the S&P500 index from Yahoo finance from
# 2005-01-01 to 2018-01-01 and compute the percentage log returns.
# Replace the zero returns with their empirical mean.

library(quantmod)

price <- getSymbols(
    "^GSPC",
    from = "2005-01-01",
    to = "2018-01-01",
    auto.assign = FALSE
)

# extract and parse data
price <- price$GSPC.Adjusted

# auto-removes the NA observation
price <- diff(log(price)) * 100

# save time-stamps
time_stamps <- index(price[-1])

# remove first observation (result of diff)
price <- as.numeric(price)[-1]

# remove zeros
price[price == 0] <- mean(price)


### Point 2)
# Estimate the GAS model you developed in point 1 using the code of point 2.
gas_fit <- estimate_tgas(price)

gas_fit$optimal_params

### Point 3)
# Compare the volatility implied by GAS and the ones from the GARCH and
# SV models you have estimated in the "Real Data" part of Exercise Set 4.

# sourcing this file does not work, so I'll just plot the process here
# source('Exam/Problem Sets/Problem Set 4/Problem Set 4.R')


plot(
    x = time_stamps,
    y = gas_fit$sigma,
    type = "l",
    xlab = "Time",
    ylab = "Volatility"
)