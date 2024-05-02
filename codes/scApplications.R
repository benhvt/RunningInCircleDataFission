library(Seurat)
library(countsplit)
library(pbapply)
library(ggplot2)
library(patchwork)
library(dplyr)
library(latex2exp)

theme_set(theme_bw())

seurat_pipeline <- function(object, 
                            cell_pop,
                            nfeature = 200, 
                            ncomp = 10,
                            theta = NULL){
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

#Load Data
TS <-  readRDS(file = "data/TS_data.rds")

#Subset only patient TSP14
TS_14 <- subset(TS, subset = donor_id == "TSP14")
dim(TS_14)

#Estimate Global overdispersion 
overdisp_all <- sctransform::vst(TS_14@assays$RNA@counts, n_genes = NULL)
theta_all <- overdisp_all$model_pars[,1]

#Test function
neutrophil <- seurat_pipeline(TS_14, cell_pop = c("neutrophil"))
neutrophil$typeI
neutrophil$umap

neutrophil_all <- seurat_pipeline(TS_14, cell_pop = c("neutrophil"),
                                  theta = theta_all)
neutrophil_all$typeI
neutrophil_all$umap

#Apply pipeline on 5 cell populations 
cell_pop_to_test <- c("neutrophil", 
                      "macrophage",
                      "monocyte",
                      "granulocyte",
                      "CD4-positive, alpha-beta T cell",
                      "memory B cell")

#Results when data thinning is performed using intra-cell type overdispersion
results <- lapply(cell_pop_to_test, function(c){
  temp <- seurat_pipeline(TS_14, cell_pop = c, nfeature = 2000)
  return(temp)}
)

#Plot Umap
lapply(1:length(results), function(l){results[[l]]$umap + ggtitle(cell_pop_to_test[l])})

#Plot type I 
typeI_full <- do.call("rbind.data.frame", lapply(results, function(l){l$typeI}))

typeI_full %>% 
  group_by(CellType) %>%
  summarise(TypeI = mean(TypeI))%>%
  ggplot() +
  aes(x=CellType, y = TypeI) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0.05, color = "darkred") +
  ylim(c(0,1))


# Plot theta
## How many theta can be estimated on all the cell population
ggvenn::ggvenn(data = list("neutrophil"=names(results[[1]]$theta), 
                           "macrophage"=names(results[[2]]$theta),
                           "monocyte" = names(results[[3]]$theta),
                           "granulocyte"=names(results[[4]]$theta),
                           "CD4-positive, alpha-beta T cell"=names(results[[5]]$theta),
                           "memory B cell"=names(results[[6]]$theta)))


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

pair_cellPop <- combn(2:(length(cell_pop_to_test)+1),2)
allPairplot <- lapply(1:ncol(pair_cellPop), function(p){
  if (colnames(theta_share_df)[pair_cellPop[1,p, drop = T]] == "memory.b"){
    lab_x <- "memory b cells"
  }
  else{
    lab_x <- colnames(theta_share_df)[pair_cellPop[1,p, drop = T]]
  }
  
  if (colnames(theta_share_df)[pair_cellPop[2,p, drop = T]] == "memory.b"){
    lab_y <- "memory b cells"
  }
  else{
    lab_y <- colnames(theta_share_df)[pair_cellPop[2,p, drop = T]]
  }
  df_temp <- data.frame(Gene = theta_share_df$Gene, 
                        Pop1 = theta_share_df[,pair_cellPop[1,p, drop = T]],
                        Pop2 = theta_share_df[,pair_cellPop[2,p, drop = T]])
  
  rmse <- round(sqrt(mean((df_temp$Pop1 - df_temp$Pop2)^2)), 2)
  df_temp$rmse <- paste0("RMSE=", rmse)
  ggplot(df_temp) + aes(x=Pop1, y = Pop2) +
    # geom_point(alpha = .5) +
    scattermore::geom_scattermore(pointsize = 4, alpha = .3) +
    geom_abline(slope = 1, intercept = 0, colour = "darkred", linetype = "dashed", linewidth = 1) +
    xlab(lab_x) +
    ylab(lab_y) +
    # geom_label(aes(x = Inf, y = Inf, label = paste("RMSE =", round(rmse, 2))), 
    #            hjust = 1, vjust = 1, size = 6, color = "white", fill = "darkred") +
    facet_grid(~rmse) +
    scale_x_log10() +
    scale_y_log10() +
    annotation_logticks(side = "bl") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    NULL
})

# do.call(cowplot::plot_grid, allPairplot)
plt <- ((allPairplot[[1]]+ ggtitle("C") + theme(plot.title = element_text(face = "bold", size = 20)))  + allPairplot[[2]] + allPairplot[[3]] + allPairplot[[4]]) + plot_layout(nrow = 1) 

source("codes/Figure4.R")

(plot_illustration) / (plt) &  #+ plot_annotation(tag_levels = "A") + plot_layout(tag_level = "new") & 
  theme(text = element_text(size = 20))
# ggsave("figures/figure4_C.pdf", width = 300, height = 100, units = "mm")
ggsave("figures/figure4.png", width = 300, height = 200, units = "mm")

#Result when data thinning is performed using global overdispersion
results_all <- lapply(cell_pop_to_test, function(c){
  temp <- seurat_pipeline(TS_14, cell_pop = c, nfeature = 2000, theta = theta_all)
  return(temp)}
)

#Plot Umap
lapply(1:length(results_all), function(l){results[[l]]$umap + ggtitle(cell_pop_to_test[l])})

#Plot type I 
typeI_full_all <- do.call("rbind.data.frame", lapply(results_all, function(l){l$typeI}))

typeI_full_all %>% 
  group_by(CellType) %>%
  summarise(TypeI = mean(TypeI))%>%
  ggplot() +
  aes(x=CellType, y = TypeI) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0.05, color = "darkred") +
  ylim(c(0,1))


# Plot theta
## How many theta can be estimated on all the cell population
ggvenn::ggvenn(data = list("neutrophil"=names(results_all[[1]]$theta), 
                           "macrophage"=names(results_all[[2]]$theta),
                           "monocyte" = names(results_all[[3]]$theta),
                           "granulocyte"=names(results_all[[4]]$theta),
                           "CD4-positive, alpha-beta T cell"=names(results_all[[5]]$theta),
                           "memory B cell"=names(results_all[[6]]$theta)))


name_theta_share_all <- Reduce(intersect, lapply(results_all, function(l){names(l$theta)}))

## Plot theta 
theta_share_all <- lapply(results_all, function(l){l$theta[name_theta_share]})
theta_share_df_all <- data.frame(Gene = theta_share_all[[1]],
                             neutrophil = theta_share_all[[1]],
                             macrophage = theta_share_all[[2]],
                             monocyte = theta_share_all[[3]],
                             granulocyte = theta_share_all[[4]],
                             CD4 = theta_share_all[[5]],
                             memory.b = theta_share_all[[6]])

pair_cellPop <- combn(2:(length(cell_pop_to_test)+1),2)
allPairplot_all <- lapply(1:ncol(pair_cellPop), function(p){
  df_temp <- data.frame(Gene = theta_share_df_all$Gene, 
                        Pop1 = theta_share_df_all[,pair_cellPop[1,p, drop = T]],
                        Pop2 = theta_share_df_all[,pair_cellPop[2,p, drop = T]])
  ggplot(df_temp) + aes(x=Pop1, y = Pop2) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, colour = "darkred") +
    xlab(colnames(theta_share_df)[pair_cellPop[1,p, drop = T]]) +
    ylab(colnames(theta_share_df)[pair_cellPop[2,p, drop = T]]) +
    scale_x_log10() +
    scale_y_log10()
  
})

do.call(cowplot::plot_grid, allPairplot_all)
