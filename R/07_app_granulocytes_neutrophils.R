# Conditional vs marginal NB data thinning on a controlled mixture of two cell
# populations (granulocytes + CD24 neutrophils, Tabula Sapiens bone marrow
# data, donor TSP14), with 250 truly DE genes (H1) and 250 non-DE genes (H0)
# -> data and results for Figure 3 and Supplementary Figures S4-S5 (Appendix A7);
# < 1 min

source("R/utils.R")
library(Seurat)
library(countsplit)

TS_14 <- load_TSP14()
cells <- subset(TS_14, cell_ontology_class %in% c("granulocyte", "cd24 neutrophil"))
counts <- counts_matrix(cells)
counts <- counts[, apply(counts, 2, var) >= 1e-2]
pval <- unlist(pbapply::pblapply(1:ncol(counts), function(p) {
  wilcox_pvalue(counts[, p], cells$cell_ontology_class)
}, cl = n_cores))
counts <- counts[, c(order(pval)[1:250], order(pval, decreasing = TRUE)[1:250])]  # H1, then H0 genes

# cells with only zeros on the 500 genes are removed (none), ordered by cell type
cell_type <- cells$cell_ontology_class
keep <- rowMeans(counts == 0) < 1
counts <- counts[c(which(keep & cell_type == "granulocyte"), which(keep & cell_type == "cd24 neutrophil")), ]
true_cell <- factor(rep(c("granulocyte", "cd24 neutrophil"), c(sum(keep & cell_type == "granulocyte"),
  sum(keep & cell_type == "cd24 neutrophil"))))
true_pval <- apply(counts, 2, wilcox_pvalue, groups = true_cell)

cluster_and_test <- function(split, seed) {
  # k-means (2 clusters) on X^(1), multivariate or gene by gene (the k-means
  # of gene p being run after set.seed(seed + p)), Wilcoxon tests on X^(2)
  Xtrain <- as.matrix(split[[1]])
  Xtest <- as.matrix(split[[2]])
  cl <- kmeans(log2(Xtrain + 1), centers = 2, nstart = 100)$cluster
  univ <- do.call(rbind, pbapply::pblapply(1:ncol(Xtrain), function(p) {
    set.seed(seed + p)
    cl_p <- kmeans(log2(Xtrain[, p] + 1), centers = 2, nstart = 100)$cluster
    data.frame(pvalue = wilcox_pvalue(Xtest[, p], cl_p), ARI = mclust::adjustedRandIndex(cl_p, true_cell))
  }, cl = n_cores))
  rbind(data.frame(univ, Clustering = "Univariate"),
    data.frame(pvalue = apply(Xtest, 2, wilcox_pvalue, groups = cl),
      ARI = mclust::adjustedRandIndex(cl, true_cell), Clustering = "Multivariate"))
}

thin <- function(X) countsplit(X, folds = 2, epsilon = 0.5, overdisps = estimate_theta(X))
set.seed(202601)
marg <- cluster_and_test(thin(counts), seed = 202601)
set.seed(202601)
by_type <- lapply(levels(true_cell)[2:1], function(type) thin(counts[true_cell == type, ]))
cond <- cluster_and_test(lapply(1:2, function(f) rbind(by_type[[1]][[f]], by_type[[2]][[f]])), seed = 202602)

results <- data.frame(rbind(cbind(cond, Thinning = "Conditional"),
  cbind(marg, Thinning = "Marginal")),
Gene = colnames(counts),
Hyp = rep(c("H1", "H0"),
  each = 250),
GeneSymbol = TS_14[["RNA"]][[]][colnames(counts), "feature_name"],
TruePvalue = true_pval)

saveRDS(list(counts = counts, true_cell = true_cell, results = results,
  n_cells = c(table(as.character(cell_type)), kept = table(true_cell))),
file.path(res_dir, "app07_granulocytes_neutrophils.rds"))
