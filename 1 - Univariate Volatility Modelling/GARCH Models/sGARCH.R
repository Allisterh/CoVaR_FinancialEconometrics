library(rugarch)
library(Rsolnp)

# reproducability
set.seed(69)


################################################################################
### Simulate from a GARCH(1, 1) process                                      ###
################################################################################

simulate_garch <- function(obs, omega, alpha, beta) {
    #' Simulates `obs` number of GARCH(1,1) observations

    # placeholder for running variables
    y <- numeric(obs)
    sigma_2 <- numeric(obs)

    # variance initialized at unconditional value
    sigma_2[1] <- omega / (1.0 - alpha - beta)

    # sample first observation
    y[1] <- sqrt(sigma_2[1]) * rnorm(1, mean = 0, sd = 1)

    for (t in seq(2, obs)) {
        # sample volatility
        sigma_2[t] <- omega + alpha * y[t - 1]^2 + beta * sigma_2[t - 1]

        # sample new observation
        y[t] <- sqrt(sigma_2[t]) * rnorm(1, mean = 0, sd = 1)
    }

    # return simulated series
    return(
        list(
            "returns" = y,
            "sigma_2" = sigma_2
        )
    )
}


garch_sim <- simulate_garch(
    obs = 10000, omega = 0.1, alpha = 0.05, beta = 0.94
)

# plot simulated series
plot(
    garch_sim$sigma_2,
    type = "l", ylab = "Conditional variance", xlab = "Time"
)
plot(
    garch_sim$returns,
    type = "l", ylab = "Log returns", xlab = "Time"
)


################################################################################
### Manually implemented GARCH(1,1) model                                    ###
################################################################################

garch_log_likelihood <- function(params, returns, ...) {
    #' Calculates log-likelihood of an GARCH(1,1) model
    #'
    #' @param params vector. Vector of parameters to use for estimation
    #' contains omega, alpha and beta.
    #' @param returns vector. Vector of returns of an asset.

    # "unpack" parameter vector
    omega <- params[1]
    alpha <- params[2]
    beta <- params[3]

    # number of observations
    obs <- length(returns)

    # placeholder for variance
    sigma_2 <- numeric(obs)

    # variance initialized at unconditional value
    sigma_2[1] <- omega / (1.0 - alpha - beta)

    # initialize first log-likelihood value
    ll <- dnorm(returns[1], sd = sqrt(sigma_2[1]), log = TRUE)

    for (t in seq(2, obs)) {
        sigma_2[t] <- omega + alpha * returns[t - 1]^2 + beta * sigma_2[t - 1]
        ll <- ll + dnorm(returns[t], sd = sqrt(sigma_2[t]), log = TRUE)
    }

    # return negative log-likelihood
    return(
        list(
            sigma_2 = sigma_2,
            neg_ll = -ll,
            ll = ll
        )
    )
}

estimate_garch <- function(series, ...) {
    #' Calculates optimal parameters for GARCH(1,1) model using maximum
    #' likelihood estimation. Takes into account constraints on parameters
    #' series optimizer (solnp) object
    #'
    #' @param series vector. Vector of returns of an asset.
    #' @param ... any. Optional arguments to forward to likelihood function.

    require(Rsolnp)

    # convert data-type of input to numeric
    series <- as.numeric(series)

    # precision constants (lower and upper)
    pre_l <- 1e-4 # 0.0001
    pre_u <- 1 - pre_l # 0.9999

    # Set start value for alpha = 0.1, beta = 0.8, and chose omega to target
    # the empirical variance by targeting the unconditional variance of the
    # GARCH model
    alpha <- 0.1
    beta <- 0.8
    omega <- var(series) * (1.0 - alpha - beta)

    # "pack" parameters in vector
    params <- c(
        omega = omega,
        alpha = alpha,
        beta = beta
    )

    # find optimal parameters using solnp solver
    # solnp allows for setting constraints
    optimizer <- solnp(
        params,
        fun = function(...) {
            return(garch_log_likelihood(...)$neg_ll)
        },

        # weak stationarity condition: alpha + beta < 1
        ineqfun = function(params, ...) {
            return(params[2] + params[3])
        }, ineqLB = pre_l, ineqUB = pre_u,

        # positivity condition
        LB = c(pre_l, pre_l, pre_l),

        # upper conditions: 0 < (alpha, beta) < 1
        UB = c(10, pre_u, pre_u),

        # supress output (run quietly)
        control = list(trace = 0),

        # argument to neg ll function
        returns = series,

        # optional arguments parsed to log-likelihood function
        ...
    )

    # extract optimal parameters
    optimal_params <- optimizer$pars

    # estimate mu and calculate optimal log-likelihood
    optimal_fit <- garch_log_likelihood(
        params = optimal_params,
        returns = series,
        ...
    )

    # calculate Bayes Information Criteria
    # \mathrm{BIC}=k \ln (n)-2 \ln (\widehat{L})
    bic <- length(params) * log(length(series)) - 2 * optimal_fit$ll

    return(
        list(
            "sigma_2" = optimal_fit$sigma_2,
            "ll" = optimal_fit$ll,
            "bic" = bic,
            "params" = optimal_params
        )
    )
}

estimate_garch(garch_sim$returns)


################################################################################
### Estimate GARCH(1,1) model using rugarch                                  ###
################################################################################

## Fit GARCH(1,1) using the rugarch package
# initialize model
garch_spec <- ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
    mean.model = list(armaOrder = c(0, 0), include.mean = FALSE)
)

# fit model
ugarchfit(garch_spec, garch_sim$returns)
