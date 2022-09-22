library(Rsolnp)

# reproducability
set.seed(69)

################################################################################
### Simulate from a ARCH(1,1) process                                       ###
################################################################################

simulate_arch <- function(obs, omega, alpha) {
    #' Simulates `obs` number of ARCH(1,1) observations

    # placeholder for running variables
    y <- numeric(obs)
    sigma_2 <- numeric(obs)

    # variance initialized at unconditional value
    sigma_2[1] <- omega / (1.0 - alpha)

    # sample first observation
    y[1] <- sqrt(sigma_2[1]) * rnorm(1, mean = 0, sd = 1)

    for (t in seq(2, obs)) {
        # sample volatility
        sigma_2[t] <- omega + alpha * y[t - 1]^2

        # sample new observation
        y[t] <- sqrt(sigma_2[t]) * rnorm(1, mean = 0, sd = 1)
    }

    # return simulated series
    return(
        list(
            "returns" = y,
            "variance" = sigma_2
        )
    )
}


# simulate ARCH(1) by setting beta = 0
arch_sim <- simulate_arch(
    obs = 10000, omega = 0.3, alpha = 0.7
)

# plot simulated series
plot(
    arch_sim$variance,
    type = "l", ylab = "Conditional variance", xlab = "Time"
)
plot(
    arch_sim$returns,
    type = "l", ylab = "Log returns", xlab = "Time"
)

################################################################################
### Manually implemented ARCH(1,1) model                                     ###
################################################################################


arch_neg_log_likelihood <- function(params, returns) {
    #' Calculates negative log-likelihood of an ARCH(1,1) model
    #' Takes parameter "guesses" and return vector as argument

    # "unpack" parameter vector
    omega <- params[1]
    alpha <- params[2]

    # number of observations
    obs <- length(returns)

    # placeholder for variance
    sigma_2 <- numeric(obs)

    # variance initialized at unconditional value
    sigma_2[1] <- omega / (1.0 - alpha)

    # initialize first log-likelihood value
    ll <- dnorm(returns[1], sd = sqrt(sigma_2[1]), log = TRUE)

    for (t in seq(2, obs)) {
        sigma_2[t] <- omega + alpha * returns[t - 1]^2
        ll <- ll + dnorm(returns[t], sd = sqrt(sigma_2[t]), log = TRUE)
    }

    # return negative log-likelihood
    return(-ll)
}

estimate_arch <- function(returns) {
    #' Calculates optimal parameters for ARCH(1) model using maximum
    #' likelihood estimation. Takes into account constraints on parameters
    #' returns optimizer (solnp) object

    # precision constants (lower and upper)
    pre_l <- 1e-4 # 0.0001
    pre_u <- 1 - pre_l # 0.9999

    # Set start value for alpha=0.1, and chose omega to target
    # the empirical variance by targeting the unconditional variance of the
    # ARCH model
    alpha <- 0.1
    omega <- var(returns) * (1.0 - alpha)

    # "pack" parameters in vector
    params <- c(omega, alpha)

    # find optimal parameters using solnp solver
    # solnp allows for setting constraints
    optimizer <- solnp(
        params,
        fun = arch_neg_log_likelihood,

        # positivity condition
        LB = c(pre_l, pre_l),

        # upper conditions: alpha < 1, omega < 10
        UB = c(10, pre_u),

        # returns vector - argument to neg ll function
        returns = returns
    )

    # return optimzer object with results
    return(optimizer)
}

estimate_arch(arch_sim$returns)

