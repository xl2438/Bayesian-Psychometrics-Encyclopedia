compute_expected_numbers <- function(theta, psi) {
  n_cat <- length(psi)
  z <- outer(psi, theta, "+")
  p <- pnorm(z)
  x_hat <- 1 * (runif(length(p), 0, 1) < p)
  x_hat_sum <- apply(x_hat, 2, sum)
  out <- table(factor(x_hat_sum, levels = 0:n_cat))
  return(as.vector(out))
}