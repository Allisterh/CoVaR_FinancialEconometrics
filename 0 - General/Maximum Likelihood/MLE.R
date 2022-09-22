# Leopoldo's example from lecture


dLambda <- 5
vY <- rexp(10000, dLambda)

nLLK <- function(dLambda, vY) {
  dnLLK <- -sum(dexp(vY, dLambda, log = TRUE))
  return(dnLLK)
}

vLambda <- seq(1, 10, by = 0.1)

vnLLK <- sapply(vLambda, nLLK, vY = vY)

plot(vLambda, vnLLK, type = "l")
abline(v = 1 / mean(vY), col = "red")

mean(1 / mean(vY))

optimizer <- optim(2, nLLK,
  method = "L-BFGS-B",
  lower = 1e-6,
  upper = 10,
  vY = vY
)

options(digits = 10)
optimizer$par
1 / mean(vY)