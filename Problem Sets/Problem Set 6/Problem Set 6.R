library(rugarch)
library(Rsolnp)

################################################################################
### Problem 2                                                                ###
################################################################################

dcc_filter <- function(etas, alpha, beta, eta_correlation) {
    #' Function to do stuff...

    # dimensions
    n <- ncol(etas)
    t <- nrow(etas)

    # array placeholders for the Q matrices
    q <- array(NA, dim = c(n, n, t))

    # array placeholder for the correlations
    r <- array(NA, dim = c(n, n, t))


    # initialize these values at unconditional correlation
    q[, , 1] <- eta_correlation
    r[, , 1] <- eta_correlation

    # compute initial log-likelihood contribution
    ll <- etas[1, , drop = 0] %*% solve(r[, , 1]) %*% t(etas[1, , drop = 0]) -
        etas[1, , drop = 0] %*% t(etas[1, , drop = 0]) + log(det(r[, , 1]))

    for (i in seq(2, t)) {

        # update with new Q matrix
        q[, , i] <- eta_correlation * (1 - alpha - beta) + alpha *
            (t(etas[i - 1, , drop = 0]) %*% etas[i - 1, , drop = 0]) +
            beta * q[, , i - 1]

        # transform correlation to ensure P.D.'ness:
        # R_{t}=\widetilde{Q}_{t}^{-1/2} Q_{t} \widetilde{Q}_{t}^{-1/2}
        r[, , i] <- diag(sqrt(1 / diag(q[, , i]))) %*% q[, , i] %*%
            diag(sqrt(1 / diag(q[, , i])))

        # calculate updated ll
        ll <- ll + etas[i, , drop = 0] %*% solve(r[, , i]) %*%
            t(etas[i, , drop = 0]) - etas[i, , drop = 0] %*%
            t(etas[i, , drop = 0]) + log(det(r[, , i]))
    }

    return(
        list(
            "ll" = -(1 / 2) * ll,
            "r" = r
        )
    )
}

estimate_dcc <- function(returns) {
    #' Estimate DCC model on `returns` using QML

    # require external packages
    require(Rsolnp)
    require(rugarch)

    # set limits for optimizer
    lower_limit <- 1e-4
    upper_limit <- 1 - lower_limit

    # specify/initialize GARCH(1, 1) model
    garch_specification <- ugarchspec(
        variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
        mean.model = list(armaOrder = c(0, 0))
    )

    # create list of specifications for multifit
    garch_specification <- multispec(
        replicate(ncol(returns), garch_specification)
    )

    # fit all models in multispec
    garch_fit <- multifit(
        multispec = garch_specification, data = returns
    )

    # compute standardized residuals and store them in `n x T` matrix
    etas <- residuals(garch_fit, standardize = TRUE)

    # extract univariate volatility processes and store them in `n x T` matrix
    garch_sigma <- sigma(garch_fit)

    # extract estimated univariate parameters
    garch_params <- coef(garch_fit)

    # calculate unconditional correlation
    eta_correlation <- cor(etas)

    # set initial parameters `alpha, beta`
    params <- c(alpha = 0.04, beta = 0.9)

    # setup `Rsolnp` object to perform ML estimation
    optimizer <- solnp(
        params,
        fun = function(params, etas, eta_correlation) {
            filter <- dcc_filter(etas, params[1], params[2], eta_correlation)
            ll <- as.numeric(filter$ll)

            return(-ll)
        },

        # stationarity condition alpha + beta < 1
        ineqfun = function(params, ...) {
            params[1] + params[2]
        }, ineqLB = lower_limit, ineqUB = upper_limit,

        # 0 < (alpha, beta) < 1
        LB = rep(lower_limit, 2), UB = rep(upper_limit, 2),

        # supress output (run quietly)
        control = list(trace = 0),

        # additional arguments
        etas = etas, eta_correlation = eta_correlation
    )

    # extract estimated parameters
    optimal_params <- optimizer$pars

    # filter dynamic correlation using estimated parameters
    optimal_filter <- dcc_filter(
        etas = etas,
        alpha = optimal_params[1],
        beta = optimal_params[2],
        eta_correlation = eta_correlation
    )


    return(
        list(
            "params" = optimal_params,
            "garch_params" = garch_params,
            "garch_volatilites" = garch_sigma
        )
    )
}

dcc <- estimate_dcc(
    dj_returns
)


################################################################################
### Problem 3                                                                ###
################################################################################

# load dataset from rugarch package
data(dji30ret)

# save returns in vector (ticker AA and AXP first 1000 observations)
dj_returns <- dji30ret
dj_returns <- dj_returns[seq(1, 1000), c("AA", "AXP")]

# look at data
head(dj_returns)


estimate_dcc(
    dj_returns
)


Estimate_DCC(dj_returns)$mCoef