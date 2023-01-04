## Compute the minimum variance portfolio
Compute_MVP <- function(mSigma) {

  iN = ncol(mSigma)
  vOnes = rep(1, iN)

  vOmega = vOnes %*% solve(mSigma)
  vOmega = vOmega/sum(vOmega)

  return(vOmega)
}