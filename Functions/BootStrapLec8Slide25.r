# Write a function to perform filtering of the volatility in the Stochastic
# Volatility (E[\exp\left(\alpha_{t}/2\right)\mid y_{1:t}) using the Bootstrap
# filter reported in slide 25 of Lecture 8.

bootstrap_filter <- function(returns, omega, phi, tau, N = 10000) {
    # function to perform bootstrap particle filter algorithm

    # number of observations
    obs <- length(returns)

    # matrix of bootstrapped particles
    alpha_bootstrap <- matrix(
        data = NA, nrow = obs, ncol = N
    )

    # vector of filtered volatility values
    volatility <- numeric(obs)

    # importance weight at each t
    weights <- numeric(N)


    # initialize values with draws from unconditional distribution
    # fill in first row
    alpha_bootstrap[1, ] <- rnorm(
        N,
        mean = omega / (1.0 - phi),
        sd = sqrt(tau**2 / (1.0 - phi**2))
    )
    # compute importance weights following notes in .pdf/.tex file
    # we draw `N` draws from gaussian distribution with varying scale
    weights <- dnorm(
        returns[1],
        mean = 0, sd = exp(alpha_bootstrap[1, ] / 2.0)
    )

    # normalize importance weights
    weights <- weights / sum(weights)

    # approximate \hat{\sigma}_{t}=E[\exp(\alpha_{t}/2)\mid y_{1:t}] for t = 1
    # this is a weighted average
    volatility[1] <- sum(exp(alpha_bootstrap[1, ] / 2) * weights)

    # draw `N` samples from `alpha_bootstrap` with replacement and
    # corresponding probabilities
    alpha_bootstrap[1, ] <- sample(
        alpha_bootstrap[1, ],
        size = N,
        replace = TRUE,
        prob = weights
    )

    for (t in seq(2, obs)) {
        # generate data
        alpha_bootstrap[t, ] <- omega + phi * alpha_bootstrap[t - 1, ] +
            rnorm(N, mean = 0, sd = tau)

        # draw weights
        weights <- dnorm(
            returns[t],
            mean = 0,
            sd = exp(alpha_bootstrap[t, ] / 2.0)
        )

        # normalize weights
        weights <- weights / sum(weights)

        # Approximate \hat{\sigma}_{t}=E[\exp(\alpha_{t}/2)\mid y_{1:t}] for t
        volatility[t] <- sum(exp(alpha_bootstrap[t, ] / 2) * weights)

        # resample
        alpha_bootstrap[t, ] <- sample(
            alpha_bootstrap[t, ],
            size = N,
            replace = TRUE,
            prob = weights
        )
    }


    return(
        list(
            volatility = volatility
        )
    )
}

