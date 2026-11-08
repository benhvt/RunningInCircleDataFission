# Marginal vs conditional NB data thinning in a two-component NB mixture where
# k-means (3 clusters) splits one component into two spurious clusters
# -> Figure 1B; < 1 min

source("R/utils.R")
library(countsplit)

n <- 100
Z <- rep(1:2, each = n / 2)
theta_mle <- function(X) apply(X, 2, function(x) npreg::theta.mle(x, mu = mean(x)))

test_spurious <- function(Xtrain, Xtest) {
  # Wilcoxon test on the first variable between the two clusters of X^(1)
  # with the closest means (i.e. originating from the same component)
  cl <- kmeans(Xtrain, centers = 3)$cluster
  m <- tapply(Xtrain[, 1], cl, mean)
  pairs <- combn(3, 2)
  keep <- cl %in% pairs[, which.min(abs(m[pairs[1, ]] - m[pairs[2, ]]))]
  wilcox_pvalue(Xtest[keep, 1], cl[keep])
}

sim_nb_mixture <- function() {
  X <- cbind(c(rnbinom(n / 2, prob = 0.5, size = 5), rnbinom(n / 2, prob = 0.4, size = 40)),
             rnbinom(n, prob = 0.5, size = 5))
  split <- lapply(countsplit(X, overdisps = theta_mle(X)), as.matrix)
  p_marg <- test_spurious(split[[1]], split[[2]])
  split <- lapply(1:2, function(z) lapply(countsplit(X[Z == z, ], overdisps = theta_mle(X[Z == z, ])), as.matrix))
  p_cond <- test_spurious(rbind(split[[1]][[1]], split[[2]][[1]]), rbind(split[[1]][[2]], split[[2]][[2]]))
  data.frame(pvalues = c(p_marg, p_cond), Fission = c("Marginal thinning", "Conditional thinning"))
}

# seeds of the original analysis
res <- run_replicates(seeds = 310123 * seq_len(n_rep), sim_nb_mixture)
saveRDS(rbind_rep(res), file.path(res_dir, "sim03_nb_thinning_mixture.rds"))
