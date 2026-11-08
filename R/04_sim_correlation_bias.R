# Impact of the correlation between variables and of a biased scale parameter
# on NB data thinning and Gaussian data fission (Appendix A6)
# -> Figure 1C, Supplementary Figures S2-S3; ~11 min on 12 cores
# p = 50 correlated count variables are simulated from a Poisson-log-normal
# distribution with the mean-variance relationship of a NB(theta = 10).

source("R/utils.R")
library(mvtnorm)
library(countsplit)

n <- 100
p <- 50
tau <- 0.5
theta <- rep(10, p)

sim_correlation_bias <- function(mu, rho, error) {
  # error: theta_hat / theta. Returns the Wilcoxon p-values of all variables
  # for NB thinning and Gaussian fission, with theta_hat ("Wrong") and theta
  # ("Oracle"), after k-means (2 clusters) on X^(1).
  R <- matrix(rho, p, p)
  diag(R) <- 1
  latent_cov <- function(th) diag(sqrt(log(1 + 1 / th))) %*% R %*% diag(sqrt(log(1 + 1 / th)))
  Sigma <- latent_cov(theta)
  Y <- matrix(rpois(n * p, sweep(exp(rmvnorm(n, sigma = Sigma)), 2, mu, "*")), n, p)
  Y_gauss <- rmvnorm(n, mean = mu, sigma = Sigma)
  test <- function(Xtrain, Xtest) {
    cl <- kmeans(Xtrain, centers = 2, nstart = 100)$cluster
    apply(Xtest, 2, wilcox_pvalue, groups = cl)
  }
  thin <- function(th) do.call(test, countsplit(Y, overdisps = th, epsilon = rep(tau, 2)))
  fission <- function(W) test(Y_gauss + tau * W, Y_gauss - (1 / tau) * W)
  data.frame(pval = c(thin(theta * error), thin(theta),
                      fission(rmvnorm(n, sigma = latent_cov(theta * error))), fission(rmvnorm(n, sigma = Sigma))),
             Variable = rep(paste0("Var", 1:p), 4),
             Method = rep(c("NB", "Gauss"), each = 2 * p),
             Estimation = rep(rep(c("Wrong", "Oracle"), each = p), 2),
             Error = error, Rho = rho)
}

grid <- expand.grid(rho = c(0, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9),
                    error = c(0.001, 0.01, 0.05, 0.1, 0.5, 0.75, 0.9, 0.95, 1, 1.01, 1.1, 1.25, 1.5, 2, 2.5, 5, 10))
res <- run_replicates(seeds = 5e6 + seq_len(n_rep), function() {
  mu <- runif(p, 2, 20)
  res <- do.call(rbind, Map(sim_correlation_bias, list(mu), grid$rho, grid$error))
  # Figure 1C: number of rejections among the p variables; Figures S2-S3:
  # p-values of the first and last variables in the displayed settings
  list(typeI = aggregate(cbind(n_reject = pval < 0.05, n_tests = 1) ~ Method + Estimation + Rho + Error, res, sum),
       pvalues = subset(res, Variable %in% c("Var1", "Var50") &
                          ((Estimation == "Oracle" & Error == 1) |
                             (Estimation == "Wrong" & Rho == 0 & Error %in% c(0.001, 0.5, 1, 1.5, 5)))))
})
saveRDS(list(typeI = rbind_rep(lapply(res, `[[`, "typeI")), pvalues = rbind_rep(lapply(res, `[[`, "pvalues"))),
        file.path(res_dir, "sim04_correlation_bias.rds"))
