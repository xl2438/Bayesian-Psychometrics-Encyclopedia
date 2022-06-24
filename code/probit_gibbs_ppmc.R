compute_polycor <- function(ind, dat) {
  x_sum <- apply(dat[-ind, ], 2, sum)
  out <- polycor::polychor(dat[ind, ], x_sum)
  return(out)
}
compute_expected_numbers <- function(theta, psi) {
  n_cat <- length(psi)
  z <- outer(psi, theta, "+")
  p <- pnorm(z)
  x_hat <- 1 * (runif(length(p), 0, 1) < p)
  x_hat_sum <- apply(x_hat, 2, sum)
  mat_2_way <- table(x_hat[1, ], x_hat[2, ])
  out_1 <- as.vector(table(factor(x_hat_sum, levels = 0:n_cat)))
  out_2 <- var(sapply(c(1:5), FUN = compute_polycor, dat = x_hat))
  out_3 <- prod(mat_2_way[c(1, 4)]) / prod(mat_2_way[c(2, 3)])
  return(c(out_1, out_2, out_3))
}

irt_1pl_gibbs <- function(holder, dat, n_itr) {
  x <- dat
  n_person <- dim(x)[2]
  n_item <- dim(x)[1]
  a <- 1.0
  psi <- rep(0, n_item)
  theta <- rnorm(n_person, 0, a)
  post_psi <- matrix(NA, nrow = n_itr, ncol = n_item)
  post_a <- rep(NA, n_itr)
  post_ppmc <- matrix(NA, nrow = n_itr, ncol = n_item + 3)
  ######################################################
  temp <- outer(psi, theta, "+")
  y <- matrix(rnorm(length(temp), temp, 1), nrow = dim(temp)[1])
  mu_0 <- 0
  sigma_0 <- 10
  mu_0_mu <- 0
  mu_0_sigma <- 10
  sigma_0_alpha <- 0.001
  sigma_0_beta <- 0.001
  alpha <- 0.001
  beta <- 0.001
  lower_bound <- ifelse(x == 1, 0, -Inf)
  upper_bound <- ifelse(x == 0, 0, Inf)
  ####################################################
  for (i in 1:n_itr) {
    ############## update psi #################
    sum_theta <- sum(theta)
    sum_y <- apply(y, 1, sum)
    temp_mu <- (mu_0 + sigma_0 * (sum_y - sum_theta)) / (1 + sigma_0 * n_person)
    temp_var <- sigma_0 / (1 + sigma_0 * n_person)
    psi <- rnorm(n_item, temp_mu, sqrt(temp_var))
    ###########################################
    ############# update theta ################
    sum_y_j <- apply(y, 2, sum)
    sum_psi <- sum(psi)
    temp_mu <- (a * (sum_y_j - sum_psi)) / (1 + a * n_item)
    temp_var <- a / (1 + a * n_item)
    theta <- rnorm(n_person, temp_mu, sqrt(temp_var))
    ###########################################
    ############# update y ###################
    temp <- outer(psi, theta, "+")
    y <- EnvStats::rnormTrunc(length(temp), temp, 1, lower_bound, upper_bound)
    y <- matrix(y, nrow = n_item)
    #########################################
    ############ update a ##################
    temp_alpha <- alpha + n_person / 2
    temp_beta <- beta + sum(theta^2) / 2
    a <- invgamma::rinvgamma(1, shape = temp_alpha, rate = temp_beta)
    ########################################
    ############ update hyperparameters ########
    temp_mu <- (sigma_0 * mu_0_mu + mu_0_sigma * sum(psi)) / (sigma_0 + n_item
      * mu_0_sigma)
    temp_sigma <- 1 / (1 / mu_0_sigma + n_item / sigma_0)
    mu_0 <- rnorm(1, temp_mu, sqrt(temp_sigma))
    ###
    temp_alpha <- sigma_0_alpha + n_item / 2;
    temp_beta <- sigma_0_beta + sum((psi - mu_0)^2) / 2;
    sigma_0 <- invgamma::rinvgamma(1, shape = temp_alpha, rate = temp_beta)
    ###################
    post_psi[i, ] <- psi
    post_a[i] <- a
    ###################
    post_ppmc[i, ] <- compute_expected_numbers(theta, psi)
  }
  return(cbind(post_psi, post_a, post_ppmc))
}