library(rugarch)
library(Rsolnp)

# replicability
set.seed(69)


################################################################################
### Simulate from a GARCH(1, 1) process                                      ###
################################################################################

f_SimGarch <- function(iT, dOmega, dAlpha, dBeta) {
    ## Simulates `iT` number of GARCH(1,1) observation

    # placeholder for running variables
    vY <- numeric(iT)
    vSigma2 <- numeric(iT)
    vZ <- rnorm(iT)

    # variance initialized at unconditional value
    vSigma2[1] <- dOmega / (1.0 - dAlpha - dBeta)

    # sample first observation
    vY[1] <- sqrt(vSigma2[1]) * rnorm(1, mean = 0, sd = 1)

    for (t in 2:iT) {
        # sample volatility
        vSigma2[t] <- dOmega + dAlpha * vY[t - 1]^2 + dBeta * vSigma2[t - 1]

        # sample new observation
        vY[t] <- sqrt(vSigma2[t]) * vZ[t]
    }

    # return simulated series
    return(
        list(
            "returns" = vY,
            "vSigma2" = vSigma2
        )
    )
}


lGARCH_SIM <- f_SimGarch(
    iT = 10000, dOmega = 0.1, dAlpha = 0.05, dBeta = 0.94
)

# plot simulated series
plot(
    lGARCH_SIM$vSigma2,
    type = "l", ylab = "Conditional variance", xlab = "Time"
)
plot(
    lGARCH_SIM$returns,
    type = "l", ylab = "Log returns", xlab = "Time"
)


################################################################################
### Manually implemented GARCH(1,1) model                                    ###
################################################################################

garch_log_likelihood <- function(params, returns, ...) {
    #' Calculates log-likelihood of an GARCH(1,1) model
    #'
    #' @param params vector. Vector of parameters to use for estimation
    #' contains dOmega, dAlpha and dBeta.
    #' @param returns vector. Vector of returns of an asset.

    # "unpack" parameter vector
    dOmega <- params[1]
    dAlpha <- params[2]
    dBeta <- params[3]

    # number of iTervations
    iT <- length(returns)

    # placeholder for variance
    vSigma2 <- numeric(iT)

    # variance initialized at unconditional value
    vSigma2[1] <- dOmega / (1.0 - dAlpha - dBeta)

    # initialize first log-likelihood value
    ll <- dnorm(returns[1], sd = sqrt(vSigma2[1]), log = TRUE)

    for (t in seq(2, iT)) {
        vSigma2[t] <- dOmega + dAlpha * returns[t - 1]^2 + dBeta * vSigma2[t - 1]
        ll <- ll + dnorm(returns[t], sd = sqrt(vSigma2[t]), log = TRUE)
    }

    # return negative log-likelihood
    return(
        list(
            vSigma2 = vSigma2,
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

    # Set start value for dAlpha = 0.1, dBeta = 0.8, and chose dOmega to target
    # the empirical variance by targeting the unconditional variance of the
    # GARCH model
    dAlpha <- 0.1
    dBeta <- 0.8
    dOmega <- var(series) * (1.0 - dAlpha - dBeta)

    # "pack" parameters in vector
    params <- c(
        dOmega = dOmega,
        dAlpha = dAlpha,
        dBeta = dBeta
    )

    # find optimal parameters using solnp solver
    # solnp allows for setting constraints
    optimizer <- solnp(
        params,
        fun = function(...) {
            return(garch_log_likelihood(...)$neg_ll)
        },

        # weak stationarity condition: dAlpha + dBeta < 1
        ineqfun = function(params, ...) {
            return(params[2] + params[3])
        }, ineqLB = pre_l, ineqUB = pre_u,

        # positivity condition
        LB = c(pre_l, pre_l, pre_l),

        # upper conditions: 0 < (dAlpha, dBeta) < 1
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
            "vSigma2" = optimal_fit$vSigma2,
            "ll" = optimal_fit$ll,
            "bic" = bic,
            "params" = optimal_params
        )
    )
}

estimate_garch(lGARCH_SIM$returns)


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
ugarchfit(garch_spec, lGARCH_SIM$returns)
