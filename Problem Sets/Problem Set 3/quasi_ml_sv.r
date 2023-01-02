quasi_ml_sv <- function(y) {
    #' This function maximizes the Quasi Log Likelihood
    #' computed by the Kalman Filter for a log-linearized
    #' stochastic volatility model.
    #'
    #' The SV model is
    #' y_t = sigma*exp(w_t/2)*z_t, z_t ~ iid N(0,1)
    #' w_t = rho * w_{t-1} + eta_t, eta_t ~ iid N(0, sigma_eta^2)
    #' It's log linearized version is
    #' zeta_t = mu + w_t + eps_t,
    #' w_t = rho * w_{t-1} + eta_t, eta_t ~ iid N(0, sigma_eta^2)
    #' We have that E[eps_t|F_{t-1}] = 0, and E[eps_t^2|F_{t-1}] = pi^2/2
    #' and mu = log(sigma^2) + E[log(z_t^2)]
    #' where E[log(z_t^2)] = -1.270376
    #'
    #' y are the "level" non-transformed returns

    # number of observations
    obs <- length(y)

    if (any(y == 0)) {
        y[y == 0] <- mean(y)
    }

    # transform returns
    y_transform <- log(y^2)
    # demean returns
    y_transform <- y_transform - mean(y_transform)

    # starting parameters
    params <- c(
        cor(y_transform[-1], y_transform[-obs]),
        var(y_transform) * 0.1
    )

    # optimize the likelihood
    optimizer <- optim(
        params,
        function(params, y_transform, obs) {
            rho <- params[1]
            sigma2_eta <- params[2]

            params_tot <- c(
                rho, sqrt(pi**2 / 2), sqrt(sigma2_eta)
            )
            kalman <- KalmanFilter_AR1plusNoise(
                y_transform, params_tot,
                Smoothing = FALSE
            )

            return(
                -kalman$dLLK
            )
        },
        method = "L-BFGS-B", lower = c(0.001, 0.001), upper = c(0.99, 1.0),
        y_transform = y_transform, obs = obs
    )

    # extract estimated parameters
    # compute the vector of total parameters
    rho <- optimizer$par[1]
    sigma2_eta <- optimizer$par[2]
    params_tot <- c(optimizer$par[1], sqrt(pi^2 / 2), sqrt(optimizer$par[2]))

    # filtering and smoothing
    kalman <- KalmanFilter_AR1plusNoise(
        y_transform, params_tot,
        Smoothing = TRUE
    )

    model_params <- c(
        # I have no idea where 1.270376 comes from
        sigma = exp((mean(log(y^2)) + 1.270376) / 2),
        rho = rho,
        sigma2_eta = sigma2_eta
    )

    return(
        list(
            params = model_params,
            kalman = kalman
        )
    )
}

