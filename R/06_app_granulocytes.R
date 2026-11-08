# NB data thinning on a homogeneous population (no true cluster): granulocytes
# of the Tabula Sapiens bone marrow data (donor TSP14), 500 most variable genes
# -> Figure 2, Type I error rates of Section 4.4, Supplementary Figure S9; < 1 min

source("R/utils.R")
library(Seurat)
library(countsplit)

TS_14 <- load_TSP14()
counts_all <- counts_matrix(subset(TS_14, cell_ontology_class == "granulocyte"))
counts <- counts_all[, order(apply(counts_all, 2, var), decreasing = TRUE)[1:500]]
counts <- counts[rowMeans(counts == 0) < 0.5, ]  # cells with >= 50% zeros removed

set.seed(20260122)
split <- lapply(countsplit(counts, folds = 2, epsilon = 0.5, overdisps = estimate_theta(counts)), as.matrix)
Xtrain <- split[[1]]
Xtest <- split[[2]]

# k-means (2 clusters) on X^(1), with all genes (multivariate) or gene by gene
# (univariate, the k-means of gene p being run after set.seed(20260122 + p)),
# and Wilcoxon tests on X^(2)
cluster <- kmeans(log2(Xtrain + 1), centers = 2, nstart = 100)$cluster
pval_univ <- unlist(pbapply::pblapply(1:ncol(Xtrain), function(p) {
  set.seed(20260122 + p)
  wilcox_pvalue(Xtest[, p], kmeans(log2(Xtrain[, p] + 1), centers = 2, nstart = 100)$cluster)
}, cl = n_cores))

genes <- data.frame(Gene = colnames(counts),
                    GeneSymbol = TS_14[["RNA"]][[]][colnames(counts), "feature_name"],
                    CorOriginal = cor(counts[, 1], counts)[1, ],     # with the first gene
                    CorThinning = cor(Xtrain[, 1], Xtest)[1, ],
                    PvalMultivariate = apply(Xtest, 2, wilcox_pvalue, groups = cluster),
                    PvalUnivariate = pval_univ)
saveRDS(list(counts = counts, Xtrain = Xtrain, Xtest = Xtest, cluster = cluster, genes = genes,
             n_cells = c(TSP14 = nrow(counts_all), kept = nrow(counts))),
        file.path(res_dir, "app06_granulocytes.rds"))
