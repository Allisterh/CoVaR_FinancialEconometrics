library(rugarch)
library(Rsolnp)
library(copula)

################################################################################
### Problem 1                                                                ###
################################################################################

### Point 1)
# Write a function to estimate the bivariate Gaussian and Student’s t copula
# model of Patton (2006a,b) detailed in slide 32/33 of lecture 11 with M = 1.
# Each marginal model should be a Student’s t GARCH and can be estimated with
# the rugarch package (to compute the probability integral transform, u_{i,t},
# you can use the pit() function on an estimated uGARCHfit object, see
# help("uGARCHfit-class")). The function should return the total likelihood of
# the model, the filtered correlations and variances and estimated parameters.

patton_filter <- function(u, copula, omega, alpha, beta, nu) {
    lambda_function <- function(x) {
        #' Lambda mapping function (modified logistic function)
        #' returns \Lambda(x)=\frac{1-e^{-x}}{1+e^{-x}}

        return(
            (1 - exp(-x)) / (1 + exp(-x))
        )
    }

    if (copula == "norm") {
        u_transform <- qnorm(u)
    }
    if (coupula == "t") {
        u_transform <- qt(u, df = nu)
    }

    # get T observations
    obs <- nrow(u)

    
}

patton_filter(5)

estimate_patton <- function(copula, returns) {
    #' Main function to estimate the Patton model with t or normal copula
    #'
    #' @param copula character. Type of Copula to estimate ("norm" or "T")
    #' @param returns matrix. Matrix (T x N) of asset returns

    # require needed packages
    require(rugarch)
    require(copula)
    require(Rsolnp)

    # set float precision
    lower_limit <- 1e-4
    upper_limit <- 1 - lower_limit


    # specify univariate GARCH(1, 1) model(s)
    garch_specification <- ugarchspec(
        mean.model = list(armaOrder = c(0, 0)),
        # std = Student's T Distribution
        distribution.model = "std"
    )

    # estimate GARCH model for each asset (column) in returns
    garch_fit <- apply(returns, 2, function(asset) {
        ugarchfit(garch_specification, asset)
    })

    # extract probability integral transform
    u <- do.call(
        cbind,
        lapply(garch_fit, function(fit) {
            pit(fit)
        })
    )
    # preserve names of returns after running lapply
    names(u) <- names(returns)

    # make sure models are robust (lower_limit < u < upper_limit)
    u[u > upper_limit] <- upper_limit
    u[u < lower_limit] <- lower_limit

    ## set starting parameters for the Patton Filter
    # calculate approximate unconditional correlation
    approx_correlation <- cor(u[, 1], u[, 2]) * 0.16

    # map approximate unconditional correlation
    omega_start <- log(approx_correlation + 1) - log(1 - approx_correlation)

    # set start parameters in vector
    params <- c(
        omega = omega_start,
        alpha = 0.04,
        beta = 0.8
    )
    if (copula == "norm") {
        optimizer <- solnp(
            params,
            fun = function(params, etas, eta_correlation) {
                filter <- patton_filter(
                    u, params[1], params[2], params[3],
                    nu = NA
                )
                # negative log-likelihood
                ll <- -as.numeric(filter$ll)

                if (!is.finite(ll) || !is.numeric(ll)) {
                    # make sure the likelihood doesn't collapse
                    ll <- 1e4
                }

                return(ll)
            },

            # stationarity condition alpha + beta < 1
            ineqfun = function(params, ...) {
                params[2] + params[3]
            }, ineqLB = lower_limit, ineqUB = upper_limit,

            # 0 < beta < 1 && -1 < alpha < 1 && -3 < omega < 3
            LB = c(-3, -upper_limit, lower_limit),
            UB = c(3, upper_limit, upper_limit),

            # supress output (run quietly)
            control = list(trace = 0),

            # additional arguments
            u = u
        )
    }

    optimal_parameters <- optimizer$pars

    optimal_parameters
}


u <- estimate_patton(
    copula = "norm", returns = dj_returns
)
cor(u)


################################################################################
### Problem 2                                                                ###
################################################################################

### Point 1)
# Consider the same couple of assets you choose in the previous exercise set

# load dataset from rugarch package
data(dji30ret)

# save returns in vector (ticker AA and AXP first 1000 observations)
dj_returns <- dji30ret
dj_returns <- dj_returns[seq(1, 1000), c("AA", "AXP")]

# look at data
head(dj_returns)