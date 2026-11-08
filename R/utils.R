# Shared settings and functions. All paths are relative to the root folder of
# the archive, which must be the working directory.

n_cores <- getOption("mc.cores", max(1, parallel::detectCores() - 1))
n_rep <- getOption("ric.n_rep", 1000)               # Monte Carlo replicates
res_dir <- getOption("ric.results_dir", "results")  # intermediate results
dir.create(res_dir, showWarnings = FALSE)

run_replicates <- function(seeds, fun) {
  # Runs fun() once per seed, in parallel (forked processes, sequential on
  # Windows). Replicate r is run after set.seed(seeds[r]), so that each
  # replicate is reproducible on its own, whatever the number of cores.
  pbapply::pblapply(seq_along(seeds), function(r) {
    set.seed(seeds[r])
    fun()
  }, cl = n_cores)
}

rbind_rep <- function(res) do.call(rbind, Map(cbind, rep = seq_along(res), res))

ttest_pvalue <- function(x, groups) {
  # NA when an estimated cluster contains a single observation (rare, n = 50)
  tryCatch(t.test(x ~ groups)$p.value, error = function(e) NA)
}

wilcox_pvalue <- function(x, groups) {
  # Two-sided Wilcoxon rank-sum test p-value, computed as by wilcox.test() in
  # R < 4.5 (used for the original analyses), whatever the version of R:
  # R 4.6 uses exact p-values for small samples with ties and rounds asymptotic
  # p-values < 1e-16 to 0, and R 4.5 changed the underflow of pnorm().
  groups <- factor(groups)
  r <- rank(x)
  in_x <- groups == levels(groups)[1]
  n.x <- sum(in_x)
  n.y <- sum(!in_x)
  W <- sum(r[in_x]) - n.x * (n.x + 1) / 2
  if (n.x < 50 && n.y < 50 && !anyDuplicated(r)) {
    p <- if (W > n.x * n.y / 2) pwilcox(W - 1, n.x, n.y, lower.tail = FALSE) else pwilcox(W, n.x, n.y)
    return(min(2 * p, 1))
  }
  n_ties <- table(r)
  sigma <- sqrt(n.x * n.y / 12 * (n.x + n.y + 1 - sum(n_ties^3 - n_ties) / ((n.x + n.y) * (n.x + n.y - 1))))
  z <- (W - n.x * n.y / 2 - sign(W - n.x * n.y / 2) * 0.5) / sigma
  if (abs(z) >= 37.5193) return(0)
  2 * min(pnorm(z), pnorm(z, lower.tail = FALSE))
}

rho_colours <- c("#FFE5B4", "#FFBFA0", "#FF8CA0", "#C478B8", "#8050A0", "#50286F", "#2B1240")  # Figs 1C, S2

# Two-component bivariate Gaussian mixture (Supplementary Figures S6-S7)
pi_mix <- c(0.5, 0.5)
mu_mix <- list(c(0, 5), c(0, 0))
Sigma_mix <- list(diag(1, 2), cbind(c(1, 0.5), c(0.5, 1)))

simulate_mixture <- function(n) {
  Z <- sample(1:2, n, replace = TRUE, prob = pi_mix)
  X <- matrix(NA, n, 2)
  for (g in 1:2) X[Z == g, ] <- mvtnorm::rmvnorm(sum(Z == g), mu_mix[[g]], Sigma_mix[[g]])
  list(Z = Z, X = X)
}

load_TSP14 <- function() {
  # Tabula Sapiens bone marrow data (CZ CELLxGENE Discover, v5.0.0), donor TSP14
  file <- "data/TS_data.rds"
  if (!file.exists(file)) {
    dir.create("data", showWarnings = FALSE)
    options(timeout = 3600)
    download.file("https://datasets.cellxgene.cziscience.com/5e736dcd-01d8-4639-805a-31fea1528be0.rds",
                  file, mode = "wb")
  }
  stopifnot(tools::md5sum(file) == "3099d66fb374dc693a8e3486284ebcba")
  subset(readRDS(file), subset = donor_id == "TSP14")
}

counts_matrix <- function(seurat_object) {
  t(as.matrix(SeuratObject::GetAssayData(seurat_object, layer = "counts")))  # cells x genes
}

estimate_theta <- function(X) {
  # gene-wise maximum likelihood estimates of the NB overdispersion
  unlist(pbapply::pblapply(seq_len(ncol(X)), function(p) {
    tryCatch(npreg::theta.mle(X[, p], mu = mean(X[, p])), error = function(e) NA)
  }, cl = n_cores))
}
