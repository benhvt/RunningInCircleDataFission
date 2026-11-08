# Reproduces all the data-based figures, tables and numbers of the article and
# of its Supplementary Materials (see README.md).
# Run from the root folder of the archive: `Rscript main.R`

RERUN <- FALSE  # TRUE: re-compute the intermediate results in results/ (~30 min
                # on 12 cores, downloads the Tabula Sapiens data); FALSE: ~1 min
options(mc.cores = max(1, parallel::detectCores() - 1))  # number of cores used

packages <- c("countsplit", "npreg", "mvtnorm", "mclust", "pbapply", "Seurat", "sctransform", "glmpca", "umap",
              "dplyr", "tidyr", "ggplot2", "patchwork", "latex2exp", "ggnewscale", "MetBrewer", "colorspace",
              "scattermore", "ggrepel")
missing <- setdiff(packages, rownames(installed.packages()))
if (length(missing) > 0) install.packages(missing)  # all on CRAN

scripts <- c("R/08_figures_main.R", "R/09_figures_supplementary.R")
if (RERUN) scripts <- c(sort(list.files("R", pattern = "^0[1-7]_", full.names = TRUE)), scripts)
for (s in scripts) {
  message("Running ", s)
  source(s, local = new.env())
}

invisible(lapply(packages, library, character.only = TRUE))
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
