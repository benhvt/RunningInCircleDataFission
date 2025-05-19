################################################################################
#                                                                              #  
#                 Applications of data thinning on scRNA-seq data              #      
#                       Benjamin Hivert - 17/04/2025                           #  
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
#
# This script contains code to reproduce the comparison between overall and intra-
# cell population overdispersion based on real Tabula Sapiens datasets on Bone Marrow 
#
#------------------------------------------------------------------------------#


# Library

library(Seurat)
library(countsplit)
library(pbapply)
library(ggplot2)
library(patchwork)
library(dplyr)
library(latex2exp)

theme_set(theme_bw())

#------------------------------------------------------------------------------#
#                               Internal functions                             #
#------------------------------------------------------------------------------#

seurat_pipeline <- function(object, 
                            cell_pop,
                            nfeature = 200, 
                            ncomp = 10,
                            theta = NULL){
  # Function to apply all the Seurat Pipeline on a scRNA-seq datasets
  
  #Params:
  #cell_pop: a vector containing cell population to extract
  #nfeautre: number of features to keep
  #ncomp: number of PC to keep for (clustering is applied on PCA results)
  #theta: the overdispersion paramaters
  
  # Subset cell type
  print(cell_pop)
  print("Pre-processing")
  object_sub <- subset(object, cell_type %in% cell_pop)
  # object_sub <- FindVariableFeatures(object_sub, nfeature=5000)
  # object_sub <- subset(object_sub, features = VariableFeatures(object_sub))
  
  print("Overdispersion estimation and data thining")
  
  # Compute overdispersion for data thining
  if (is.null(theta)){
    overdisp <- sctransform::vst((object_sub@assays$RNA@counts), n_genes = NULL)
    theta <- overdisp$model_pars[,1]
  }
  split <- countsplit(as.matrix(t(as.matrix(object_sub@assays$RNA@counts))[, which(rownames(object_sub@assays$RNA@counts) %in% names(theta))]), overdisps = theta)
  Xtrain <- as.matrix(split[[1]])
  Xtest <- as.matrix(split[[2]])
  
  print("Clustering on the train")
  # Clustering
  sc.train <- CreateSeuratObject(counts = t(Xtrain))
  sc.train <- NormalizeData(sc.train, normalization.method = "LogNormalize", scale.factor = 10000)
  
  sc.train <- FindVariableFeatures(sc.train, selection.method = "vst", nfeatures = nfeature)
  all.genes.train <- rownames(sc.train)
  sc.train <- ScaleData(sc.train, features = all.genes.train)
  
  sc.train <- RunPCA(sc.train, features = VariableFeatures(object = sc.train))
  sc.train <- FindNeighbors(sc.train, dims = 1:ncomp)
  sc.train <- FindClusters(sc.train, resolution = 0.5)
  sc.train <- RunUMAP(sc.train, dims = 1:ncomp)
  umapTrain <- DimPlot(sc.train, reduction = "umap") + 
    xlab("UMAP1") + 
    ylab("UMAP2") + 
    ggtitle(cell_pop, subtitle = TeX(r'($X^{train}$)'))
  
  print("DEA on the test")
  #DEA 
  sc.test <- CreateSeuratObject(counts=t(Xtest))
  sc.test <- NormalizeData(sc.test)
  sc.test <- ScaleData(sc.test)
  sc.test <- RunPCA(sc.test, features = VariableFeatures(object = sc.train))
  sc.test <- RunUMAP(sc.test, dims = 1:ncomp)
  Idents(sc.test) <- Idents(sc.train)
  umapTest <- DimPlot(sc.test, reduction = "umap") + 
    xlab("UMAP1") + 
    ylab("UMAP2") +
    ggtitle(cell_pop, subtitle = TeX(r'($X^{test}$)'))
  
  K  <- length(table(Idents(sc.test)))
  if (K == 1){
    typeI <- data.frame(TypeI = NA, 
                        Pair = "Cluster 0",
                        CellType = cell_pop)
  }
  else{
    pair <- combn(0:(K-1), 2)
    marker.gene <- pblapply(1:ncol(pair), function(p){
      markers <- FindMarkers(sc.test, ident.1=pair[1,p], ident.2 = pair[2,p])
      return(data.frame(TypeI = mean(markers$p_val_adj<0.05),
                        Pair = paste("Cluster", pair[1,p], "vs", "Cluster", pair[2,p], sep = " ")))
    }, cl = 4)
    typeI <- do.call("rbind.data.frame", marker.gene)
    typeI$CellType <- cell_pop
  }
  umap_fin <- (umapTrain+umapTest + plot_layout(guides = "collect")) 
  return(list(theta = theta,
              umap = umap_fin,
              typeI = typeI))
}

#------------------------------------------------------------------------------#
#                               Application                                    #
#------------------------------------------------------------------------------#

#Load Data
TS <-  readRDS(file = "data/TS_data.rds")

#Subset only patient TSP14
TS_14 <- subset(TS, subset = donor_id == "TSP14")
dim(TS_14)

#------------------------------------------------------------------------------#
#                     Marginal overdispersion estimation                       #

overdisp_all <- sctransform::vst(TS_14@assays$RNA@counts, n_genes = NULL)
theta_all <- overdisp_all$model_pars[,1]

#------------------------------------------------------------------------------#
#           Conditionnal (intra cell-pop) overdispersion estimation            #


cell_pop_to_test <- c("neutrophil", 
                      "macrophage",
                      "monocyte",
                      "granulocyte",
                      "CD4-positive, alpha-beta T cell",
                      "memory B cell")

#Results when data thinning is performed using intra-cell type overdispersion
results <- pblapply(cell_pop_to_test, function(c){
  temp <- seurat_pipeline(TS_14, cell_pop = c, nfeature = 2000)
  return(temp)},
  cl = 5
)


#Note: In some cell populations, it is impossible to estimate overdispersion. So,
# comparision of gene-specific overdispersion between cell population will be based
# on the 8,333 genes where an overdispersion estimate is available 

name_theta_share <- Reduce(intersect, lapply(results, function(l){names(l$theta)}))

## Plot theta 
theta_share <- lapply(results, function(l){l$theta[name_theta_share]})
theta_share_df <- data.frame(Gene = theta_share[[1]],
                             neutrophil = theta_share[[1]],
                             macrophage = theta_share[[2]],
                             monocyte = theta_share[[3]],
                             granulocyte = theta_share[[4]],
                             CD4 = theta_share[[5]],
                             memory.b = theta_share[[6]])

write.csv(theta_share_df, file = "results/Application_CellPopulationOverdispersion.csv",
          row.names = FALSE)

