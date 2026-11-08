# Type I error of the post-clustering t-test after univariate Gaussian data
# fission with a biased variance (Appendix A5)
# -> Supplementary Figure S1 (A: varying sigma, B: varying n); < 1 min on 12 cores

source("R/utils.R")

tau <- 0.4

sim_variance_bias <- function(n, sigma, sigma_hat) {
  X <- rnorm(n, mean = 0, sd = sigma)
  W <- rnorm(n, mean = 0, sd = sigma_hat)
  cl <- kmeans(X + tau * W, centers = 2)$cluster
  data.frame(pvalues = ttest_pvalue(X - (1 / tau) * W, cl), sigma = sigma, sigma_hat = sigma_hat, n = n)
}

grids <- list(
  S1A = expand.grid(sigma_hat = sort(c(seq(0, 4, length.out = 50), 2, 0.1)), sigma = c(0.1, 0.5, 1, 2), n = 100),
  S1B = expand.grid(sigma_hat = seq(0.8, 1.8, length = 50), n = c(50, 100, 200, 500, 1000), sigma = 1)
)
for (k in 1:2) {
  grid <- grids[[k]]
  res <- run_replicates(seeds = (2 + k) * 1e6 + seq_len(n_rep), function() {
    do.call(rbind, Map(sim_variance_bias, grid$n, grid$sigma, grid$sigma_hat))
  })
  saveRDS(rbind_rep(res), file.path(res_dir, paste0("sim02_variance_bias_", names(grids)[k], ".rds")))
}
