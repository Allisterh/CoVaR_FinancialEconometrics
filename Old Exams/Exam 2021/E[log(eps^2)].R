

# Computation of E[log(eps^2)], eps iid N(0,1)

# 1) Analytical solution
digamma(1/2) - log(1/2)

# 2) Direct Monte Carlo
iN = 500
vEps = rnorm(iN)
mean(log(vEps^2))

# 3) Solve the integral numerically
integrate(function(dEps) {
  log(dEps^2)*dnorm(dEps)
}, lower = -5.001, upper = 5)

# 4) Importance Sampling
iN = 50000
vEps = rt(iN, 4)

mean(log(vEps^2)*dnorm(vEps)/dt(vEps, 4))







