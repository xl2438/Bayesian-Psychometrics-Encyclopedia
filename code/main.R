source("/home/xiang/Insync/xliu30@gmail.com/Google Drive/projects/Bayesian_Psychometrics/probit_gibbs_ppmc.R")
dat <- ltm::LSAT
dim(dat)
dat <- t(dat)

library(parallel)
fit <- mclapply(X = c(1:4), FUN = irt_1pl_gibbs, dat = dat, n_itr = 6000,
  mc.cores = 1)
temp <- fit
for (i in seq_len(length(temp))) {
  temp[[i]][, 6] <- sqrt(temp[[i]][, 6])
  temp[[i]][, -6] <- -temp[[i]][, -6] / temp[[i]][, 6]
}

str(temp)
library(coda)
chain1 <- mcmc(temp[[1]], start = 4001, thin = 1)
chain2 <- mcmc(temp[[2]], start = 4001, thin = 1)
chain3 <- mcmc(temp[[3]], start = 4001, thin = 1)
chain4 <- mcmc(temp[[4]], start = 4001, thin = 1)
fit_chains <- mcmc.list(chain1, chain2, chain3, chain4)
autocorr.plot(fit_chains)
plot(fit_chains)
gelman.diag(fit_chains)
summary(fit_chains)
ci_limits <- HPDinterval(fit_chains, prob = 0.95)[[1]]
str(summary(fit_chains))
post_tab <- cbind(summary(fit_chains)[[1]][, -c(3, 4)], ci_limits)
row.names(post_tab) <- c("b1", "b2", "b3", "b4", "b5", "a")
colnames(post_tab)
post_tab
stargazer::stargazer(as.matrix(post_tab))

dat_sum_score <- table(apply(dat, 2, sum))
png("raw_score_coverage.png", width = 800, height = 800)
plotrix::plotCI(x = c(0:5), y = as.numeric(dat_sum_score), li = ci_limits[-c
  (1:6), 1], ui = ci_limits[-c(1:6), 2], xlab = "raw sum score", ylab = "count")
dev.off()

compute_polycor <- function(ind, dat) {
  x_sum <- apply(dat[-ind, ], 2, sum)
  out <- polycor::polychor(dat[ind, ], x_sum)
  return(out)
}
var(sapply(c(1:5), FUN = compute_polycor, dat = dat))
compute_polycor(1, dat)

post_var <- temp[[1]][4001:8000, 13]
png("ppmc_polycor.png", width = 800, height = 800)
hist(post_var, freq = FALSE, xlab = "Var of polychoric correlations", main =
  NULL)
abline(v = var(sapply(c(1:5), FUN = compute_polycor, dat = dat)))
dev.off()

post_OR <- temp[[1]][4001:8000, 14]
mat_2_way <- table(dat[1, ], dat[2, ])
png("ppmc_or.png", width = 800, height = 800)
hist(post_OR, freq = FALSE, xlab = "odds ratio", main = NULL)
abline(v = prod(mat_2_way[c(1, 4)]) / prod(mat_2_way[c(2, 3)]))
dev.off()

plot(fit[[1]][, 1] ~ c(1:8000), type = "l")
for (i in 2:4) {
  lines(fit[[i]][, 1] ~ c(1:8000), col = i)
}

b <- c(-2, 0, 2)
theta <- seq(from = -4, to = 4, by = 0.01)
p <- pnorm(outer(theta, b, "-"))
dim(p)
png("ICC.png")
plot(p[, 1] ~ theta, type = "l", xlab = expression(theta), ylab = expression(P
  (X == 1)))
for (i in 2:3) {
  lines(p[, i] ~ theta, lty = i)
}
legend("bottomright", legend = c("b=-2", "b=0", "b=2"), lty = c(1:3))
dev.off()

x <- seq(from = -3, to = 3, by = 0.01)
png("prior_normal.png")
plot(dnorm(x, sd = 1) ~ x, type = "l", xlab = expression(psi[j]), ylab = "Prior
  Density")
lines(dnorm(x, sd = 2) ~ x, lty = 2)
lines(dnorm(x, sd = 10) ~ x, lty = 3)
legend("topright", legend = c(expression(sigma[psi[j]]^2 == 1), expression(sigma
  [psi[j]]^2 == 4), expression(sigma[psi[j]]^2 == 100)), lty = c(1, 2, 3))
dev.off()