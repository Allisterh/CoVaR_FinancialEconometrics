# load libraries
library(quantmod, quietly = TRUE)
library(Rsolnp)

# supress warnings in quantmod package
options("getSymbols.warning4.0" = FALSE)


################################################################################
### Question 1 - Computational Part                                          ###
################################################################################

### Point a)
# Write a function to estimate the GAMMA-GAS model of the previous point using
# the Maximum Likelihood estimator. The function should accept a vector of
# observations and return the estimated parameters, the filtered means μt, for
# t = 1,...,T, and the log likelihood evaluated at its maximum value. Assume
# that omega, alpha, and beta are the intercept, score coefficient and
# autoregressive coefficient of the GAS process. Impose the following
# constraints during the optimization: omega in [−0.5,0.5],
# alpha in [0.001,1.5], beta in [0.01,0.999], a in [0.1,300].


ggas_log_likelihood <- function(y, params, constraint = FALSE) {
    #' Evaluate log-likelihood of GAMMA-GAS model
    #'
    #' @description This function evaluates the log-likelihood given data and
    #' the vector of parameters. Is built to be a input to the solnp object
    #' in `estimate_gas()`
    #'
    #' @param y vector. Vector of time-series observations to evaluate
    #' likelihood on.
    #' @param params vector. Vector of parameters for the GAS process.
    #' Contains omega, alpha, beta and a.
    #' @param constraint boolean. Determines whether to constain the model
    #' by setting \mu_t = \mu. Defaults to FALSE.

    # unpack vector of parameters
    omega <- params[1]
    alpha <- params[2]
    beta <- params[3]
    a <- params[4]

    # number of observations
    obs <- length(y)

    # placeholders for processes
    mu_tilde <- numeric(obs)
    mu <- numeric(obs)

    # initialize model (mu_tilde = unconditional mean)
    mu_tilde[1] <- omega / (1.0 - beta)
    # apply exponential link function
    mu[1] <- exp(mu_tilde[1])

    # initial log-likelihood contribution
    ll <- -lgamma(a) + a * log(a) + (a - 1.0) * log(y[1]) - a * log(mu[1]) -
        a * (y[1] / mu[1])


    for (t in seq(2, obs)) {
        # run recursion and calculate the process for mu and mu_tilde along
        # with all log-likelihood contributions
        if (constraint) {
            mu_tilde[t] <- mu_tilde[1]
            mu[t] <- exp(mu_tilde[1])
        } else {
            mu_tilde[t] <- omega + alpha * (sqrt(a) * (y[t - 1] /
                exp(mu_tilde[t - 1]) - 1.0)) + beta * mu_tilde[t - 1]

            # apply exponential link function
            mu[t] <- exp(mu_tilde[t])
        }

        # initial log-likelihood contribution
        ll <- ll - lgamma(a) + a * log(a) + (a - 1.0) * log(y[t]) -
            a * log(mu[t]) - a * (y[t] / mu[t])
    }

    return(
        list(
            "mu" = mu,
            "ll" = ll,
            "neg_ll" = -ll
        )
    )
}


estimate_ggas <- function(series, ...) {
    #' Estimates GAMMA-GAS model using Maximum Likelihood estimator
    #'
    #' @description This function uses constaint optimization to estimate
    #' optimal paramters for the GAMMA-GAS model. It depends on the Rsolnp
    #' package to solve optimization problem.
    #'
    #' @param series vector. Vector of time-series observations to estimate GAS
    #' model on.
    #' @param ... any. Accepts optional arguments parsed to log-likelihood
    #' function. This is for example whether to constain the model.

    # require necessary packages
    require(Rsolnp)

    # convert input to numeric vector
    series <- as.numeric(series)

    # initialize parameters at reasonable values within constraints
    init_omega <- 0.0005
    init_alpha <- 0.5
    init_beta <- 0.5
    init_a <- 5

    # place in parameter vector
    params <- c(
        omega = init_omega,
        alpha = init_alpha,
        beta = init_beta,
        a = init_a
    )

    # optimal paramters using solnp solver
    # allows for setting constraints
    optimizer <- solnp(
        pars = params,
        fun = function(...) {
            return(ggas_log_likelihood(...)$neg_ll)
        },

        # parameter constraints
        LB = c(-0.5, 0.001, 0.01, 0.1),
        UB = c(0.5, 1.5, 0.999, 300),

        # supress output (run quietly)
        control = list(trace = 0),

        # default argument to log-likelihood function
        y = series,

        # optional arguments parsed to log-likelihood function
        ...
    )

    # extract optimal parameters
    optimal_params <- optimizer$pars

    # estimate mu and calculate optimal log-likelihood
    optimal_fit <- ggas_log_likelihood(
        series, optimal_params, ...
    )

    # calculate Bayes Information Criteria
    # \mathrm{BIC}=k \ln (n)-2 \ln (\widehat{L})
    bic <- length(params) * log(length(series)) - 2 * optimal_fit$ll

    return(
        list(
            "mu" = optimal_fit$mu,
            "ll" = optimal_fit$ll,
            "bic" = bic,
            "params" = optimal_params
        )
    )
}


################################################################################
### Question 1 - Empirical Analysis                                          ###
################################################################################

### Point a)
# Download the VIX index from Yahoo finance using the quantmod package.
# Consider the period from "2010-01-01" to "2019-01-01". Your series is the
# one reported in the column named "Adjusted".

# date range
date_min <- "2010-01-01"
date_max <- "2019-01-01"

# Create dummy dataframe for first ticker in ticker-list
vix <- getSymbols(
    Symbols = c("^VIX"),
    env = parent.frame(),
    reload.Symbols = FALSE,
    from = date_min,
    to = date_max,
    verbose = FALSE,
    warnings = TRUE,
    src = "yahoo",
    symbol.lookup = TRUE,
    auto.assign = FALSE
)[, paste("VIX", "Adjusted", sep = ".")]


### Point b)
# Estimate the GAS model you derived in the previous exercise.

ggas_unconstrained <- estimate_ggas(
    series = vix,
    constraint = FALSE
)

### Point c)
# Consider the constraint versions of the model with \mu_{t} = \mu.
# Using BIC, choose the best specification between the static
# and time–varying specifications.

ggas_constrained <- estimate_ggas(
    series = vix,
    constraint = TRUE
)

# print BIC for both models
print(
    paste(
        c("Constrained BIC:", round(ggas_constrained$bic, 2)),
        sep = ""
    )
)
print(
    paste(
        c("Unconstrained BIC:", round(ggas_unconstrained$bic, 2)),
        sep = ""
    )
)


# visual inspection plot model and VIX index
plot_colors <- c("black", "red", "blue")

plot(
    x = index(vix),
    y = as.numeric(vix),
    type = "l",
    col = plot_colors[1],
    ylim = c(min(vix) * 0.8, max(vix) * 1.2), # 20% +/- vix to y-axis
    xlab = "Time",
    ylab = "Volatility",
    lwd = 2
)
lines(
    x = index(vix),
    y = ggas_unconstrained$mu,
    col = plot_colors[2]
)
lines(
    x = index(vix),
    y = ggas_constrained$mu,
    col = plot_colors[3],
    lwd = 3
)

# increase legend scaling to 150%
op <- par(cex = 1.5)

legend(
    "topright",
    legend = c(
        "VIX index",
        "Unconstrained GAMMA-GAS",
        "Constrained GAMMA-GAS"
    ),
    col = plot_colors,
    lty = 1,
    fonts
)


### Point d)
# Write a function to estimate the MEM model of Engle and Gallo (2006).
# Impose these constraints on the parameters of the MEM model:
# kappa in [0.1,10], eta in [0.01,0.99], phi in [0.01,0.99], and a ∈[0.1,300].


mem_log_likelihood <- function(y, params) {
    #' Evaluate log-likelihood of MEM model
    #'
    #' @description This function evaluates the log-likelihood given data and
    #' the vector of parameters. Is built to be a input to the solnp object
    #' in `estimate_gas()`
    #'
    #' @param y vector. Vector of time-series observations to evaluate
    #' likelihood on.
    #' @param params vector. Vector of parameters for the MEM process.
    #' Contains kappa, eta, phi and a.

    # unpack vector of parameters
    kappa <- params[1]
    eta <- params[2]
    phi <- params[3]
    a <- params[4]

    # number of observations
    obs <- length(y)

    # placeholders for processes
    mu <- numeric(obs)

    # initialize model (mu = unconditional mean)
    mu[1] <- kappa / (1.0 - eta - phi)

    # initial log-likelihood contribution
    ll <- -lgamma(a) + a * log(a) + (a - 1.0) * log(y[1]) - a * log(mu[1]) -
        a * (y[1] / mu[1])


    for (t in seq(2, obs)) {
        # run recursion and calculate the process for mu along
        # with all log-likelihood contributions
        mu[t] <- kappa + eta * y[t - 1] + phi * mu[t - 1]

        # initial log-likelihood contribution
        ll <- ll - lgamma(a) + a * log(a) + (a - 1.0) * log(y[t]) -
            a * log(mu[t]) - a * (y[t] / mu[t])
    }

    return(
        list(
            "mu" = mu,
            "ll" = ll,
            "neg_ll" = -ll
        )
    )
}

estimate_mem <- function(series, ...) {
    #' Estimates MEM model using Maximum Likelihood estimator
    #'
    #' @description This function uses constaint optimization to estimate
    #' optimal paramters for the MEM model. It depends on the Rsolnp
    #' package to solve optimization problem.
    #'
    #' @param series vector. Vector of time-series observations to estimate GAS
    #' model on.
    #' @param ... any. Accepts optional arguments parsed to log-likelihood
    #' function.

    # require necessary packages
    require(Rsolnp)

    # convert input to numeric vector
    series <- as.numeric(series)

    # initialize parameters at reasonable values within constraints
    init_kappa <- 2
    init_eta <- 0.2
    init_phi <- 0.2
    init_a <- 5

    # place in parameter vector
    params <- c(
        kappa = init_kappa,
        eta = init_eta,
        phi = init_phi,
        a = init_a
    )

    # optimal paramters using solnp solver
    # allows for setting constraints
    optimizer <- solnp(
        pars = params,
        fun = function(...) {
            return(mem_log_likelihood(...)$neg_ll)
        },

        # parameter constraints
        LB = c(0.1, 0.01, 0.01, 0.1),
        UB = c(10, 0.99, 0.99, 300),

        # supress output (run quietly)
        control = list(trace = 0),

        # default argument to log-likelihood function
        y = series,

        # optional arguments parsed to log-likelihood function
        ...
    )

    # extract optimal parameters
    optimal_params <- optimizer$pars

    # estimate mu and calculate optimal log-likelihood
    optimal_fit <- mem_log_likelihood(
        series, optimal_params, ...
    )

    # calculate Bayes Information Criteria
    # \mathrm{BIC}=k \ln (n)-2 \ln (\widehat{L})
    bic <- length(params) * log(length(series)) - 2 * optimal_fit$ll

    return(
        list(
            "mu" = optimal_fit$mu,
            "ll" = optimal_fit$ll,
            "bic" = bic,
            "params" = optimal_params
        )
    )
}


mem <- estimate_mem(
    series = vix
)

### Point e)
# Compare the GAS and the MEM models: which one reports the highest likelihood?

# print likelihood for all three models
print(
    paste(
        c("Constrained GAS LL:", round(ggas_constrained$ll, 2)),
        sep = ""
    )
)
print(
    paste(
        c("Unconstrained GAS LL", round(ggas_unconstrained$ll, 2)),
        sep = ""
    )
)
print(
    paste(
        c("MEM LL:", round(mem$ll, 2)),
        sep = ""
    )
)


### Point f)
# Create a 3 ×1 figure containing: i) the original series, ii) the
# filtered conditional mean, and iii) the conditional variance of
# the selected model.

par(mfrow = c(3, 1), cex = 0.85)
plot(
    x = index(vix),
    y = as.numeric(vix),
    type = "l",
    col = plot_colors[1],
    ylim = c(min(vix) * 0.8, max(vix) * 1.2), # 20% +/- vix to y-axis
    xlab = "Time",
    ylab = "Volatility",
    lwd = 2
)
title(main = "Vix")

plot(
    x = index(vix),
    y = mem$mu,
    type = "l",
    col = plot_colors[2],
    ylim = c(min(mem$mu) * 0.8, max(mem$mu) * 1.2), # 20% +/- vix to y-axis
    xlab = "Time",
    ylab = "Volatility",
    lwd = 2
)
title(main = "MEM Model Conditional Mean")

plot(
    x = index(vix),
    y = mem$mu**2 / mem$params["a"],
    type = "l",
    col = plot_colors[3],
    xlab = "Time",
    ylab = "Volatility",
    lwd = 2
)
title(main = "MEM Model Conditional Variance")

# reset plotting settings
par(mfrow = c(1, 1))


################################################################################
### Question 2                                                               ###
################################################################################

### Point a)
# Estimate the DCC and CCC models on the series of GSPC and DJI returns. You
# cannot use the rugarch package for GARCH estimation. For both the DCC and CCC:

# start by downloading data
download_tickers <- function(tickers, date_min, date_max) {
    #' Function to download "Adjusted" returns from Yahoo finance

    # Create dummy dataframe for first ticker in ticker-list
    returns <- getSymbols(
        Symbols = tickers[1],
        env = parent.frame(),
        reload.Symbols = FALSE,
        from = date_min,
        to = date_max,
        verbose = FALSE,
        warnings = TRUE,
        src = "yahoo",
        symbol.lookup = TRUE,
        auto.assign = FALSE
    )[, paste(gsub("\\^", "", tickers[1]), "Adjusted", sep = ".")]


    # get remaining tickers in loop
    for (ticker in tickers[-1]) {
        returns_vec <- getSymbols(
            Symbols = c(ticker),
            env = parent.frame(),
            reload.Symbols = FALSE,
            from = date_min,
            to = date_max,
            verbose = FALSE,
            warnings = TRUE,
            src = "yahoo",
            symbol.lookup = TRUE,
            auto.assign = FALSE
        )[, paste(gsub("\\^", "", ticker), "Adjusted", sep = ".")]

        # append the column to the data-frame
        returns$ticker <- returns_vec
    }

    names(returns) <- gsub("\\^", "", tickers)

    return(
        returns
    )
}

pct_log_returns <- function(level_returns) {
    #' calculates % log returns and returns object of same shape as input
    #' math: y_{t}=(\ln(p_{t})-\ln(p_{t-1}))\cdot 100

    log_returns <- diff(log(level_returns)) * 100

    # drop NA values (first observation)
    log_returns <- na.omit(log_returns)

    return(
        log_returns
    )
}

# tickers to download
indices <- c("^DJI", "^GSPC")

# download data
returns <- download_tickers(
    tickers = indices,
    date_min = "2010-01-01",
    date_max = "2019-01-01"
)

# calculate percentage log returns
returns <- pct_log_returns(returns)

# inspect data
head(returns)


# source sGARCH file
source(
    "Exam/Old Exams/Exam 2019/MySolution/sGARCH.R"
)


dcc_filter <- function(etas, alpha, beta, eta_correlation, constraint = FALSE) {
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


estimate_dcc <- function(series) {
    #'
    #'

    # set limits for optimizer
    lower_limit <- 1e-4
    upper_limit <- 1 - lower_limit

    # require external packages
    require(Rsolnp)

    model_spec <- function(...) {
        #' Function to estimate GARCH(1, 1) model to extract:
        #'   - standardized residuals
        #'   - univariate volatility processes
        #'   - GARCH parameters
        #'
        #' @param ... any. Expects series vector from parent function

        # estimate GARCH(1,1) model for each asset in `series`
        garch_fit <- apply(series, 2, estimate_garch)

        # extract estimated parameters of GARCH models in matrix
        garch_params <- do.call(cbind, lapply(garch_fit, function(fit) {
            fit$params
        }))

        # extract univariate volatility processes
        garch_sigma <- do.call(cbind, lapply(garch_fit, function(fit) {
            fit$sigma_2
        }))

        # extract likelihood
        garch_likelihood <- do.call(cbind, lapply(garch_fit, function(fit) {
            fit$ll
        }))

        # calculate D matrix to derive eta and corr(eta)
        d_matrix <- apply(garch_sigma, 1, function(x) {
            diag(x)
        })

        # placeholder for eta time-series
        etas <- matrix(data = NA, nrow = nrow(series), ncol = ncol(series))
        colnames(etas) <- names(series)

        for (t in seq(1, nrow(series))) {
            etas[t, ] <- series[t, ] %*%
                solve(sqrt(matrix(d_matrix[, t], ncol = 2, nrow = 2)))
        }

        # calculate unconditional correlation
        eta_correlation <- cor(etas)

        return(
            list(
                garch_sigma = garch_sigma,
                garch_params = garch_params,
                garch_likelihood = garch_likelihood,
                etas = etas,
                eta_correlation = eta_correlation
            )
        )
    }

    # estimate GARCH model
    garch_fit <- model_spec()

    # set initial parameters
    params <- c(alpha = 0.04, beta = 0.9)

    # optimize log-likelihood function
    optimizer <- solnp(
        params,
        fun = function(params, etas, eta_correlation) {
            filter <- dcc_filter(
                etas = etas,
                alpha = params[1],
                beta = params[2],
                eta_correlation = eta_correlation
            )
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
        etas = garch_fit$etas, eta_correlation = garch_fit$eta_correlation
    )

    # extract estimated parameters
    optimal_params <- optimizer$pars

    # filter dynamic correlation using estimated parameters
    optimal_filter <- dcc_filter(
        etas = garch_fit$etas,
        alpha = optimal_params[1],
        beta = optimal_params[2],
        eta_correlation = garch_fit$eta_correlation
    )

    return(
        list(
            "params" = optimal_params,
            "garch_params" = garch_fit$garch_params,
            "garch_sigma" = garch_fit$garch_sigma,
            "garch_likelihood" = garch_fit$garch_likelihood
        )
    )
}

dcc <- estimate_dcc(returns)

dcc$garch_likelihood

dcc$params
