# Marginal vs conditional Gaussian data fission in a two-component mixture
# -> Supplementary Figures S6 (ideal scenario) and S7 (adverse scenario)
# ~15 min on 12 cores

source("R/utils.R")
library(mvtnorm)

Sigma_marg <- (function() {  # marginal covariance of the mixture
  mu_bar <- Reduce("+", Map("*", pi_mix, mu_mix))
  Reduce("+", Map(function(p, m, S) p * (S + tcrossprod(m)), pi_mix, mu_mix, Sigma_mix)) - tcrossprod(mu_bar)
})()

sim_fission <- function(n, tau, adverse) {
  # Ideal scenario: fission with the true (marginal or component-specific)
  # covariances, k-means (2 clusters) on all observations of X^(1).
  # Adverse scenario: fission with empirical covariances, k-means on the first
  # component only (two spurious clusters).
  # t-tests between the estimated clusters on each variable of X^(2).
  sim <- simulate_mixture(n)
  Z <- sim$Z
  X <- sim$X
  keep <- if (adverse) Z == 1 else rep(TRUE, n)
  fission_test <- function(W, type) {
    cl <- kmeans((X + tau * W)[keep, ], centers = 2, nstart = 100)$cluster
    data.frame(pvalues = apply((X - (1 / tau) * W)[keep, ], 2, ttest_pvalue, groups = cl),
               Fission = type, Variable = c("X1", "X2"),
               ARI = if (adverse) NA else mclust::adjustedRandIndex(cl, Z), n = n, tau = tau)
  }
  res_marg <- fission_test(rmvnorm(n, sigma = if (adverse) cov(X) else Sigma_marg), "Marginal")
  W_cond <- matrix(NA, n, 2)
  for (g in 1:2) {
    W_cond[Z == g, ] <- rmvnorm(sum(Z == g), sigma = if (adverse) cov(X[Z == g, ]) else Sigma_mix[[g]])
  }
  rbind(res_marg, fission_test(W_cond, "Conditional"))
}

grid <- expand.grid(n = c(50, 100, 250, 500, 1000, 5000, 10000),
                    tau = c(seq(0.1, 1, by = 0.05), seq(1.5, 5, by = 0.5)))
for (adverse in c(FALSE, TRUE)) {
  res <- run_replicates(seeds = (1 + adverse) * 1e6 + seq_len(n_rep), function() {
    do.call(rbind, Map(sim_fission, grid$n, grid$tau, adverse))
  })
  saveRDS(rbind_rep(res), file.path(res_dir, if (adverse) "sim01_adverse.rds" else "sim01_ideal.rds"))
}
