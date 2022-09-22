setwd("/Users/tobiasbrammer/Library/Mobile Documents/com~apple~CloudDocs/Documents/Aarhus Uni/7. semester/Financial Econometrics/Problem Sets/Problem Set 1")

################################################################################
### Problem 1                                                                ###
################################################################################

# load list of tickers
tickers <- read.csv("DJITicker.csv", sep = ";")

# delete unlisted firms
tickers <- tickers[!(tickers$Symbol == "UTX" | tickers$Symbol == "DWDP"), ]


################################################################################
### Problem 2                                                                ###
################################################################################

library(quantmod, quietly = TRUE)

date_min <- "2021-04-30"
date_max <- "2022-09-15"

# Create dummy dataframe for first ticker in ticker-list
returns_data <- getSymbols(
    Symbols = tickers$Symbol[1],
    env = parent.frame(),
    reload.Symbols = FALSE,
    from = date_min,
    to = date_max,
    verbose = FALSE,
    warnings = TRUE,
    src = "yahoo",
    symbol.lookup = TRUE,
    auto.assign = FALSE
)[, paste(tickers$Symbol[1], "Adjusted", sep = ".")]

# remove "Adjusted" from ticker column-name
names(returns_data)[1] <- tickers$Symbol[1]


# get remaining tickers in loop
for (ticker in tickers$Symbol[-1]) {
    returns <- getSymbols(
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
    )[, paste(ticker, "Adjusted", sep = ".")]
    # append the column to the data-frame
    returns_data$ticker <- returns
    # set name of column to the ticker-name
    names(returns_data)[ncol(returns_data)] <- ticker}

# print last return observations
tail(returns_data)

pct_log_returns <- function(level_returns) {
    #' calculates % log returns and returns object of same shape as input
    #' math: y_{t}=(\ln(p_{t})-\ln(p_{t-1}))\cdot 100
    return(
        diff(log(level_returns)) * 100
    )
}

returns_data <- pct_log_returns(returns_data)
# show returns of all tickers 
plot(returns_data, type = 'l', lwd = 0.75)

################################################################################
### Problem 3                                                                ###
################################################################################

library(moments)

# placeholder for moments
desc_stat <- matrix(data = NaN, nrow = dim(tickers)[1], ncol = 7)

# set names of columns and rows
colnames(desc_stat) <- c(
    "mean", "median", "variance", "kurtosis", "skewness", "rho", "rho2"
)
rownames(desc_stat) <- tickers$Symbol

# helper function to calculate acf
autocorrelation <- function(returns, exponent = 1) {
    #' returns the autocorrelation coefficient for lag = 1
    # `[2]` due to acf calculates both lag=0 and lag=1
    return(
        acf(
            returns^exponent,
            lag = 1,
            na.action = na.pass,
            plot = FALSE
        )$acf[2]
    )
}

desc_stat[, "mean"] <- colMeans(returns_data, na.rm = TRUE)
desc_stat[, "median"] <- apply(returns_data, 2, median, na.rm = TRUE)
desc_stat[, "variance"] <- apply(returns_data, 2, var, na.rm = TRUE)
desc_stat[, "kurtosis"] <- apply(returns_data, 2, kurtosis, na.rm = TRUE)
desc_stat[, "skewness"] <- apply(returns_data, 2, skewness, na.rm = TRUE)
desc_stat[, "rho"] <- apply(returns_data, 2, autocorrelation)
desc_stat[, "rho2"] <- apply(returns_data, 2, autocorrelation, exponent = 2)

# print descriptive statistics
desc_stat


################################################################################
### Problem 4                                                                ###
################################################################################

snp <- getSymbols(
    Symbols = c("^GSPC"),
    env = parent.frame(),
    reload.Symbols = FALSE,
    from = date_min,
    to = date_max,
    verbose = FALSE,
    warnings = TRUE,
    src = "yahoo",
    symbol.lookup = TRUE,
    auto.assign = FALSE
)[, paste("GSPC", "Adjusted", sep = ".")]

# remove "Adjusted" from column name
names(snp)[1] <- "GSPC"

# calculate log returns
snp_return <- pct_log_returns(snp)
plot(snp_return, type = 'l', lwd=0.75, col = 'cornflowerblue')


fit_capm <- function(asset, factor) {
    # clean data
    asset <- asset[!is.na(asset)]
    factor <- factor[!is.na(factor)]

    # estimate OLS
    model_fit <- lm(asset ~ factor)

    # return matrix of coefficents
    capm <- matrix(nrow = 3, ncol = 1)
    capm[1, 1] <- model_fit$coefficients[1]
    capm[2, 1] <- model_fit$coefficients[2]
    capm[3, 1] <- (summary(model_fit)$sigma)^2

    return(
        capm
    )
}

# calculate CAPM estimates
capm_stat <- t(apply(returns_data, 2, fit_capm, factor = snp_return))
colnames(capm_stat) <- c("alpha", "beta", "mse")

capm_stat


################################################################################
### Problem 5                                                                ###
################################################################################

library(rugarch)
library(Rsolnp)

# reproducability
set.seed(69)


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
            "variance" = sigma_2
        )
    )
}


# simulate ARCH(1) by setting beta = 0
garch_sim <- simulate_garch(
    obs = 10000, omega = 0.3, alpha = 0.7, beta = 0
)

# plot simulated series
plot(
    garch_sim$variance,
    type = "l", ylab = "Conditional variance", xlab = "Time"
)
plot(
    garch_sim$returns,
    type = "l", ylab = "Log returns", xlab = "Time"
)


################################################################################
### Problem 6                                                                ###
################################################################################

arch_neg_log_likelihood <- function(params, returns) {
    #' Calculates negative log-likelihood of an GARCH(1,1) model
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

        # upper conditions: 0 < (alpha, beta) < 1
        UB = c(10, pre_u),

        # returns vector - argument to neg ll function
        returns = returns
    )

    # return optimzer object with results
    return(optimizer)
}

estimate_arch(garch_sim$returns)


################################################################################
### Problem 7                                                                ###
################################################################################

set.seed(69)

bootstraps <- 500
obsercations <- c(200, 500, 1000)

for (b in seq(1, bootstraps)) {
    print(b)
}