################################################################################
#                                                                              #  
#            Applications of data thinning on scRNA-seq data (2)               #      
#                       Benjamin Hivert - 21/01/2026                           #  
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
#
# This script contains code to reproduce the analysis based on the homogeneous
# granulocytes populations from Tabula Sapiens bones Marrow 
#
#------------------------------------------------------------------------------#

#-- Packages 
library(Seurat)
library(SingleCellExperiment)
library(countsplit)
library(openxlsx)

# Load Data (data could be downloaded from cellxgenes)
TS_bone_marrow <-  readRDS(file = "data/TS_data.rds")

# Create an homogeneous data from one cell population of one donnor
TS_14 <- subset(TS_bone_marrow, subset = donor_id == "TSP14")
positive_TS_14 <- subset(TS_14, cell_ontology_class %in% c("granulocyte"))


# To help homogeneity, only the 500 highest variable genes were kept and cell
# with too much 0 were removed
raw_counts <- as.matrix(t(positive_TS_14@assays$RNA@counts))
var_counts <- apply(raw_counts, 2, var)

order_var <- sort(var_counts, index.return = TRUE, decreasing = TRUE)

nGene <- 500
raw_counts_reduced <- raw_counts[, order_var$ix[1:nGene]]


prop0_counts <- apply(raw_counts_reduced, 1, function(x){mean(x==0)})
boxplot(prop0_counts)

thresh_0 <- 0.5
raw_counts_reduced_filter <- raw_counts_reduced[which(prop0_counts<thresh_0),]
dim(raw_counts_reduced_filter)

## Estimataion by MLE
theta_mle_filter <- pbsapply(1:ncol(raw_counts_reduced_filter), function(p){
  theta.mle(raw_counts_reduced_filter[,p], mu = mean(raw_counts_reduced_filter[,p]))
}, cl = 5) #+ c(runif(250, 2.5, 10), runif(250, 0, 0.5))
names(theta_mle_filter) <- colnames(raw_counts_reduced_filter)

## Data thinning 
set.seed(20260122)
split_reduced <- countsplit(X=raw_counts_reduced_filter,
                            folds = 2, epsilon = 0.5, overdisps = theta_mle_filter)
Xtrain_reduced <- as.matrix(split_reduced[[1]])
Xtest_reduced <- as.matrix(split_reduced[[2]])

## Multivariate Clustering 
clust_train_reduced <- kmeans(log2(Xtrain_reduced +1), centers = 2, nstart = 100)$cluster

## P-values from wilcoxon test 
pval_multi <- pbsapply(1:ncol(Xtrain_reduced), function(p){
  wilcox.test(Xtest_reduced[,p]~clust_train_reduced)$p.value
}, cl = 6)

## Univariate appraoch to demonstrate impact of correlation between variables 
pval_univ <- pbsapply(1:ncol(Xtrain_reduced), function(p){
  train <- Xtrain_reduced[,p]
  test <- Xtest_reduced[,p]
  
  clust_train_univ <- kmeans(log2(train +1), centers = 2, nstart = 100)$cluster
  wilcox.test(test~clust_train_univ)$p.value
}, cl = 5)



cor_or <- cor(raw_counts_reduced_filter[,1], raw_counts_reduced_filter)
cor_thin <- cor(as.matrix(Xtrain_reduced)[,1], as.matrix(Xtest_reduced))
res_df <- data.frame(OriginalCorWithGene1 = as.numeric(cor_or),
                     ThinningCorWithGene1TrainAndTest = as.numeric(cor_thin),
                     Gene = colnames(cor_or),
                     PValuesAfterMultivariateClustering = pval_multi,
                     PvaluesAfterUnivariateClustering = pval_univ)
# Save results 
wb <- createWorkbook()
addWorksheet(wb, "RawCountsFilter")
addWorksheet(wb, "Xtrain")
addWorksheet(wb, "Xtest")
addWorksheet(wb, "Cluster")
addWorksheet(wb, "Results")

writeData(wb, sheet = "RawCountsFilter", raw_counts_reduced_filter)
writeData(wb, sheet = "Xtrain", Xtrain_reduced)
writeData(wb, sheet = "Xtest", Xtest_reduced)
writeData(wb, sheet = "Cluster", data.frame(Cluster = clust_train_reduced))
writeData(wb, sheet = "Results", res_df)

saveWorkbook(wb, "results/ResultsOnScRNASeq.xlsx", overwrite = TRUE)

