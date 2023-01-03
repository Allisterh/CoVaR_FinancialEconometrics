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
