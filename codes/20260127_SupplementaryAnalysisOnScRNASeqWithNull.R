################################################################################
# Applications of data thinning on scRNA-seq data (2)
# Homogeneous granulocyte populations – Tabula Sapiens Bone Marrow
#
# Benjamin Hivert – 21/01/2026
################################################################################

#------------------------------------------------------------------------------
# 0. Packages & global options
#------------------------------------------------------------------------------

library(Seurat)
library(SingleCellExperiment)
library(org.Hs.eg.db)
library(countsplit)
library(pbapply)
library(mclust)
library(ggplot2)
library(FactoMineR)
library(npreg)
library(umap)
library(dplyr)
library(latex2exp)
library(patchwork)
set.seed(2026)

#------------------------------------------------------------------------------
# 1. Load data
#------------------------------------------------------------------------------

TS_bone_marrow <- readRDS("data/TS_data.rds")

#------------------------------------------------------------------------------
# 2. Subsetting: donor & cell types
#------------------------------------------------------------------------------

TS_14 <- subset(TS_bone_marrow, subset = donor_id == "TSP14")

cells_of_interest <- c("granulocyte", "cd24 neutrophil")
TS_14_pos <- subset(TS_14, cell_ontology_class %in% cells_of_interest)

#------------------------------------------------------------------------------
# 3. Gene filtering: top variable genes & significant genes 
#------------------------------------------------------------------------------

raw_counts <- t(as.matrix(TS_14_pos@assays$RNA@counts))  # cells x genes
gene_var   <- apply(raw_counts, 2, var)
n_genes <- 250
top_genes <- order(gene_var, decreasing = TRUE)[1:n_genes]
counts_reduced <- raw_counts[, -which(gene_var < 1e-2)]
all_pval <- pbsapply(1:ncol(counts_reduced), function(p){
  wilcox.test(counts_reduced[,p]~TS_14_pos$cell_ontology_class)$p.value
  }, cl = 6)
idx <- which(all_pval > 0)
H0_genes <- order(all_pval, decreasing = TRUE)[1:n_genes]
H1_genes <- order(all_pval, decreasing = FALSE)[1:n_genes]
# H1_genes <- idx[order(all_pval[idx])][1:250]

counts_reduced2 <- counts_reduced[, c(H1_genes, H0_genes)]  

#-----------------------------------#-----------------------------------#------------------------------------------------------------------------------
# 4. Cell filtering: proportion of zeros
#------------------------------------------------------------------------------
split_by_type <- split(
  as.data.frame(counts_reduced2),
  TS_14_pos$cell_ontology_class
)

var_group <- lapply(split_by_type, function(x) {
  apply(x, 2, function(v) var(v))
})
lapply(var_group, boxplot)


prop_zero <- lapply(split_by_type, function(x) {
  apply(x, 1, function(v) mean(v == 0))
})
lapply(prop_zero, boxplot)

thresh_0 <- 1

counts_macro  <- split_by_type$`granulocyte`[prop_zero$`granulocyte` < thresh_0, ]
counts_neutro <- split_by_type$`cd24 neutrophil`[prop_zero$`cd24 neutrophil` < thresh_0, ]

counts_filt <- rbind(counts_macro, counts_neutro)

true_cell <- factor(
  c(rep("granulocyte", nrow(counts_macro)),
    rep("cd24 neutrophil", nrow(counts_neutro)))
)

## Save resulting data sets
counts_filt %>% 
  data.frame() %>%
  mutate(TrueCell = true_cell) %>%
  write.csv("results/ResultsOnScRNASeq2CellPopulations.csv",
            row.names = FALSE)
  
#------------------------------------------------------------------------------
# 5. Exploratory analysis (PCA / UMAP)
#------------------------------------------------------------------------------
pca_true <- PCA(log2(counts_filt + 1), 
                graph = FALSE,
                ncp = 50)
clust_true <- kmeans(log2(counts_filt + 1), centers = 2, nstart = 100)$cluster

ggplot(data.frame(pca_true$ind$coord,
                  Cell = true_cell,
                  Cluster = factor(clust_true))) +
  aes(Dim.1, Dim.2, colour = Cell) +
  geom_point()

umap_true <- umap(log2(counts_filt + 1))

ggplot(data.frame(umap_true$layout,
                  Cell = true_cell)) +
  aes(X1, X2, colour = Cell) +
  geom_point()


# true pvalues 
true_pval <- pbsapply(1:ncol(counts_filt), function(p){
  wilcox.test(counts_filt[,p]~true_cell)$p.value
}, cl = 6)

#------------------------------------------------------------------------------
# 6. Overdispersion estimation (NB MLE)
#------------------------------------------------------------------------------
estimate_theta <- function(X) {
  pbsapply(seq_len(ncol(X)), function(p) {
    res <- try(theta.mle(X[, p], mu = mean(X[, p])), silent = TRUE)
    if (inherits(res, "try-error")) NA else res
  }, cl = 6)
}

theta_marg   <- estimate_theta(counts_filt)
theta_macro  <- estimate_theta(counts_macro)
theta_neutro <- estimate_theta(counts_neutro)

#------------------------------------------------------------------------------
# 7. Marginal thinning
#------------------------------------------------------------------------------
set.seed(202601)
split_marg <- countsplit(
  X = counts_filt,
  folds = 2,
  epsilon = 0.5,
  overdisps = theta_marg
)

Xtrain_marg <- as.matrix(split_marg[[1]])
Xtest_marg  <- as.matrix(split_marg[[2]])

# Multivariate clustering
clust_marg <- kmeans(log2(Xtrain_marg + 1), centers = 2, nstart = 100)$cluster
ARI_marg <- adjustedRandIndex(clust_marg, true_cell)

# Multivariate testing
pval_multi_marg <- pbsapply(seq_len(ncol(Xtrain_marg)), function(p) {
  wilcox.test(Xtest_marg[, p] ~ clust_marg)$p.value
}, cl = 6)

# Univariate illustration
univ_marg <- pblapply(seq_len(ncol(Xtrain_marg)), function(p) {
  cl <- kmeans(log2(Xtrain_marg+1)[, p], centers = 2, nstart = 100)$cluster
  data.frame(
    pvalue = wilcox.test(Xtest_marg[, p] ~ cl)$p.value,
    ARI    = adjustedRandIndex(cl, true_cell),
    Gene = colnames(Xtrain_marg)[p],
    TruePvalue = true_pval[p]
  )
}, cl = 6)

univ_marg_df <- do.call(rbind, univ_marg) %>%
  mutate(Clustering = "Univariate")

multi_marg_df <- data.frame(ARI = ARI_marg,
                            pvalue = pval_multi_marg,
                            Gene = colnames(Xtrain_marg),
                            Clustering = "Multivariate",
                            TruePvalue = true_pval)
marg_df <- rbind.data.frame(univ_marg_df, 
                            multi_marg_df) %>%
  mutate(Thinning = "Marginal")

#------------------------------------------------------------------------------
# 8. Conditional thinning
#------------------------------------------------------------------------------
set.seed(202601)
split_macro <- countsplit(counts_macro,  folds = 2, epsilon = 0.5,
                          overdisps = theta_macro)
split_neutro <- countsplit(counts_neutro, folds = 2, epsilon = 0.5,
                           overdisps = theta_neutro)

Xtrain_cond <- rbind(split_macro[[1]], split_neutro[[1]])
Xtest_cond  <- rbind(split_macro[[2]], split_neutro[[2]])

# Multivariate clustering
clust_cond <- kmeans(log2(Xtrain_cond + 1), centers = 2, nstart = 100)$cluster
ARI_cond <- adjustedRandIndex(clust_cond, true_cell)

# Multivariate testing
pval_multi_cond <- pbsapply(seq_len(ncol(Xtrain_cond)), function(p) {
  wilcox.test(Xtest_cond[, p] ~ clust_cond)$p.value
}, cl = 6)

# Univariate illustration
univ_cond <- pblapply(seq_len(ncol(Xtrain_cond)), function(p) {
  cl <- kmeans(log2(Xtrain_cond[, p] + 1), centers = 2, nstart = 100)$cluster
  data.frame(
    pvalue = wilcox.test(Xtest_cond[, p] ~ cl)$p.value,
    ARI    = adjustedRandIndex(cl, true_cell),
    Gene = colnames(Xtrain_cond)[p],
    TruePvalue = true_pval[p]
  )
}, cl = 6)

univ_cond_df <- do.call(rbind, univ_cond) %>%
  mutate(Clustering = "Univariate") %>%
  mutate(Thinning = "Conditional")

#------------------------------------------------------------------------------
# 9. Bind results
#------------------------------------------------------------------------------

multi_cond_df <- data.frame(ARI = ARI_cond,
                            pvalue = pval_multi_cond,
                            Gene = colnames(Xtrain_cond),
                            Clustering = "Multivariate",
                            Thinning = "Conditional",
                            TruePvalue = true_pval)

cond_df <- rbind.data.frame(univ_cond_df, 
                            multi_cond_df)





results_df <- rbind.data.frame(cond_df, 
                               marg_df)

results_df <- results_df %>%
  mutate(GeneSymbol = mapIds(org.Hs.eg.db,
                             keys = Gene,
                             column = "SYMBOL",
                             keytype = "ENSEMBL",
                             multiVals = "first")) 

write.csv(results_df, file = "results/SupplementaryResultsApplicationOnScRNASeqWithNull.csv")
