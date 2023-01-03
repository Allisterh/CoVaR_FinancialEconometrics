# Kalman filter and smoother for the state space:
# Y_t         = Z * alpha_t + D*eps_t, eps_t ~ N(0, S)
# alpha_{t+1} = T * alpha_t + H*eta_t, eta_t ~N(0, Q)
#
# Y_t is p x 1
# Z   is p x m
# S   is p x p
# T   is m x m
# H   is m x l
# Q   is l x l
# a1 is the initialization for a
# P1 is the initialization for P

## Here is the explanation for the code above:
# The function kalman_filter() takes the following arguments:
#     1. mY: The observed data matrix of size p * n
#     2. mZ: The matrix of the deterministic terms of size p * m
#     3. mS: The matrix of the deterministic terms of size p * p
#     4. mT: The matrix of the deterministic terms of size m * m
#     5. mH: The matrix of the deterministic terms of size r_int * m
#     6. mQ: The matrix of the deterministic terms of size r_int * r_int
#     7. a1: The initial state vector of length m
#     8. P1: The initial covariance matrix of size m * m
#     9. mD: The matrix of the deterministic terms of size p * r_int, set to 1 by default
#     10. Smoothing: A boolean value that specifies if the smoothing should be performed, set to true by default

# The function kalman_filter() returns a list of the following elements:
#     1. v: The prediction error matrix of size p * n
#     2. a_filt: The filtered state mean matrix of size m * n
#     3. a_pred: The predicted state mean matrix of size m * n
#     4. P_filt: The filtered state variance matrix of size m * m * n
#     5. P_pred: The predicted state variance matrix of size m * m * n
#     6. F: The variance of the prediction error matrix of size p * p * n
#     7. K: The kalman gain matrix of size m * p * n
#     8. N: The error smoothing matrix of size m * m * n
#     9. a_smoot: The smoothed state mean matrix of size m * n
#     10. V: The smoothed state variance matrix of size m * m * n
#     11. L: The error smoothing matrix of size m * m * n
#     12. eps_smoot: The error smoothing matrix of size p * n
#     13. eta_smoot: The error smoothing matrix of size r_int * n
#     14. vLLK: The log likelihood contribution vector of length n

kalman_filter <- function(mY, mZ, mD = matrix(1,1,1),mS, mT, mH, mQ, a1, P1, Smoothing = TRUE) {
    n <- ncol(mY)
    p <- nrow(mY)
    r_int <- ncol(mH)
    m <- length(a1)
    
    v <- matrix(0, p, n)
    F <- array(0, dim = c(p, p, n))
    K <- array(0, dim = c(m, p, n))
    a_filt <- matrix(0, m, n)
    a_pred <- matrix(0, m, n)
    P_filt <- array(0, dim = c(m, m, n))
    P_pred <- array(0, dim = c(m, m, n))

    r <- matrix(0, m, n)
    N <- array(0, dim = c(m, m, n))
    a_smoot <- matrix(0, m, n)
    V <- array(0, dim = c(m, m, n))
    L <- array(0, dim = c(m, m, n))

    eps_smoot <- matrix(0, p, n)
    eta_smoot <- matrix(0, r_int, n)
    vLLK <- numeric(n)

    # initialise
    v[, 1] <- mY[, 1]
    a_filt[, 1] <- a1
    a_pred[, 1] <- a1
    P_filt[, , 1] <- P1
    P_pred[, , 1] <- P1

    HQH <- mH %*% mQ %*% t(mH)
    # Constant in the log likelihood contribution, from slide 23
    dC <- -0.5 * (n * p * 1.0) * log(pi * 2.0)

    # Filtering recursions
    for (t in 1:n) {
        # one step ahead prediction mean, from slide 23
        v[, t] <- mY[, t] - mZ %*% a_pred[, t]
        # one step ahead prediction variance, from slide 23
        F[, , t] <- mZ %*% P_pred[, , t] %*% t(mZ) + mS %*% mD %*% t(mD)

        # kalman gain, K_t = P_t * Z_t' * inv(F_t)
        K[, , t] <- mT %*% P_pred[, , t] %*% t(mZ) %*% solve(F[, , t])

        ## Calculate the posterior mean and variance
         # The posterior adjust the prior with the new information
         # filtered state mean E[alpha_t |Y_1:t], a_filt_t = a_pred_t + K_t * v_t
         # From slide 20
         a_filt[, t] <- a_pred[, t] + K[, , t] %*% v[, t]
         # filtered state variance Var[alpha_t |Y_{1:t}], P_filt_t = P_pred_t - K_t * P_pred_t
         P_filt[, , t] <- P_pred[, , t] %*% (1 - K[, , t]) 

        # likelihood contribution. 
        # This relation is not from slide 22 of the lecture notes, but from slide 23
        vLLK[t] <- (log(det(as.matrix(F[, , t]))) + c(v[, t] %*% solve(F[, , t]) * v[, t]))

        # Predicted states
        if (t < n) {
            # Predicted state mean E[alpha_{t+1}|Y_{1:t}], see slide 22
            a_pred[, t + 1] <- mT %*% a_pred[, t] + K[, , t] %*% v[, t]
            # Predicted state variance Var[alpha_{t+1}|Y_{1:t}], see slide 22
            P_pred[, , t + 1] <- mT %*% P_pred[, , t] %*% t((mT - K[, , t] %*% mZ)) + HQH
        }
    }
    # Smoothing recursion from 
    if (Smoothing) {
        for (t in n:2) {
            # Error smoothing matrix, L_t = P_t * Z_t' * inv(F_t) from slide 22
            L[, , t] <- mT - K[, , t] %*% mZ
            # Weighted sum of innovations from slide 27, r_t-1 = inv(F_t) * v_t + L_t * r_t
            r[, t - 1] <- solve(F[, , t]) %*% v[, t] + t(L[, , t]) %*% r[, t]
            # Variance of smoothed state recursion from slide 28
            N[, , t - 1] <- solve(F[, , t]) %*% mZ + t(L[, , t]) %*% N[, , t] %*% L[, , t]
            # Smoothed state mean E[alpha_t | Y_{1:n}], from slide 27
            a_smoot[, t] <- a_pred[, t] + P_pred[, , t] %*% r[, t - 1]
            # Smoothed state variance Var[alpha_t | Y_{1:n}], from slide 28
            V[, , t] <- P_pred[, , t] - P_pred[, , t] %*% N[, , t - 1] %*% P_pred[, , t]
            
            # Not that we have smoothed states for the entire model, we can solve for the shocks
            # eps_t = S * inv(F_t) * v_t - K_t' * r_t
            eps_smoot[, t] <- mS %*% (solve(F[, , t]) %*% v[, t] - t(K[, , t]) %*% r[, t])
            # eta_t = Q * inv(F_t) * v_t - K_t' * r_t
            eta_smoot[, t] <- mQ %*% t(mH) %*% r[, t]
        }
    }

    KF <- list()

    KF[["v"]] <- v
    KF[["a_filt"]] <- a_filt
    KF[["a_pred"]] <- a_pred
    KF[["P_filt"]] <- P_filt
    KF[["P_pred"]] <- P_pred
    KF[["F"]] <- F
    KF[["K"]] <- K
    KF[["N"]] <- N
    KF[["a_smoot"]] <- a_smoot
    KF[["V"]] <- V
    KF[["L"]] <- L
    KF[["eps_smoot"]] <- eps_smoot
    KF[["eta_smoot"]] <- eta_smoot
    KF[["vLLK"]] <- vLLK
    KF[["dLLK"]] <- -dC - 0.5 * sum(vLLK)

    return(KF)
}

# Create function to map AR(1) process to Kalman filter
# Kalman filter and smoother for the model
# y_t = alpha_t + sigma * eps_t
# alpha_{t+1} = phi * alpha_t + eta * xi_t
# vPar is the vector of parameters vPar = (phi, sigma, eta)
# note that the vPar contains the standard deviations (sigma and eta)
# and not the variances (sigma^2, eta^2).
# vY is the vector of observations
KalmanFilter_AR1plusNoise <- function(vY, vPar, Smoothing = FALSE) {
    dPhi <- vPar[1]
    dSigma <- vPar[2]
    dEta <- vPar[3]

    mY <- matrix(vY, nrow = 1)
    mZ <- matrix(1, 1, 1)
    mS <- matrix(dSigma^2, 1, 1)
    mT <- matrix(dPhi, 1, 1)
    mH <- matrix(1, 1, 1)
    mQ <- matrix(dEta^2, 1, 1)
    a1 <- 0
    mP1 <- matrix(dEta^2 / (1 - dPhi^2), 1, 1)

    return(
        kalman_filter(mY, mZ, mS, mT, mH, mQ, a1, mP1, Smoothing)
    )
}
