# Gene-wise NB overdispersion estimated within 6 cell populations of the
# Tabula Sapiens bone marrow data (donor TSP14) with sctransform::vst() (which
# uses a fixed internal seed) -> Supplementary Figure S8; ~1 min on 6 cores

source("R/utils.R")
library(Seurat)

TS_14 <- load_TSP14()
cell_pops <- c(neutrophil = "neutrophil", macrophage = "macrophage", monocyte = "monocyte",
               granulocyte = "granulocyte", CD4 = "CD4-positive, alpha-beta T cell", memory.b = "memory B cell")
theta <- pbapply::pblapply(cell_pops, function(pop) {
  counts <- GetAssayData(subset(TS_14, cell_type == pop), layer = "counts")
  sctransform::vst(counts, n_genes = NULL, verbosity = 0)$model_pars[, "theta"]
}, cl = n_cores)

# genes with an overdispersion estimate in all populations
genes <- Reduce(intersect, lapply(theta, names))
saveRDS(data.frame(Gene = genes, sapply(theta, `[`, genes)), file.path(res_dir, "app05_overdispersion.rds"))
