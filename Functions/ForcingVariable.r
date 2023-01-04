
# this function returns u_t
ForcingVariable <- function(dY, dMu, dA) {
  dU = sqrt(dA)*(dY/dMu - 1)
  return(dU)
}
