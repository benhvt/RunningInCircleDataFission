# Practical limitations of information-splitting approaches for post-clustering inference

Code and intermediate results reproducing all the data-based figures, tables
and numbers of *Practical limitations of information-splitting approaches for
post-clustering inference, with implications for differential testing in
single-cell RNA-seq data* (B. Hivert, D. Agniel, R. Thiébaut, B. P. Hejblum)
and of its Supplementary Materials.
Zenodo: [10.5281/zenodo.19083400](https://doi.org/10.5281/zenodo.19083400). Licence: MIT.

## How to run

From the root folder of the archive (the working directory: all paths are
relative to it):

```bash
Rscript main.R
```

`main.R` installs the missing R packages (all on CRAN) and runs the scripts of
`R/` in order. By default, it produces the figures
(`figures/`) and tables (`tables/`) from the intermediate results provided in
`results/` (~1 min). With `RERUN <- TRUE`, it first re-computes all the
intermediate results (~30 min on 12 cores, ~2 h on 4 cores), downloading the
data if needed. Each script can also be run on its own (`Rscript R/<script>.R`).

| Script | Intermediate results / outputs | Article |
|---|---|---|
| `R/01_sim_gaussian_fission.R` | `results/sim01_ideal.rds`, `results/sim01_adverse.rds` | Figures S6, S7 |
| `R/02_sim_variance_bias.R` | `results/sim02_variance_bias_S1A.rds`, `..._S1B.rds` | Figure S1 |
| `R/03_sim_nb_thinning_mixture.R` | `results/sim03_nb_thinning_mixture.rds` | Figure 1B |
| `R/04_sim_correlation_bias.R` | `results/sim04_correlation_bias.rds` | Figures 1C, S2, S3 |
| `R/05_app_overdispersion.R` | `results/app05_overdispersion.rds` | Figure S8 |
| `R/06_app_granulocytes.R` | `results/app06_granulocytes.rds` | Figure 2, Section 4.4, Figure S9 |
| `R/07_app_granulocytes_neutrophils.R` | `results/app07_granulocytes_neutrophils.rds` |Figure 3, Figures S4, S5 (Appendix A7) |
| `R/08_figures_main.R` | `figures/Figure1.pdf`, `figures/Figure2.pdf`, `tables/Numbers_in_text.csv` | Figures 1-2, numbers of Sections 2.2.1 and 4.4 |
| `R/09_figures_supplementary.R` | `figures/FigureS1.pdf` to `figures/FigureS9.pdf`, `tables/Table_A7_performances.csv` | Figures S1-S9, Appendix A7 |

`R/utils.R` contains the shared functions. Figure S10 (`figures/FigureS10.png`)
is a drawing, and Tables S1-S2 are theoretical.

## Data

Tabula Sapiens bone marrow scRNA-seq data (v5.0.0, CZ CELLxGENE Discover,
CC BY 4.0), donor TSP14, from
<https://datasets.cellxgene.cziscience.com/5e736dcd-01d8-4639-805a-31fea1528be0.rds>
(277 MB, MD5 `3099d66fb374dc693a8e3486284ebcba`). This public file is not
included: it is downloaded automatically into `data/` when needed. The
intermediate results contain the subsets of the data used in the figures.

## Reproducibility

- **Seeds**: replicate `r` of each simulation is run after `set.seed(k * 1e6 + r)`
  (`k` = 1, 2 for Figures S6, S7; 3, 4 for Figures S1A, S1B; 5 for Figures 1C,
  S2, S3) or `set.seed(310123 * r)` (Figure 1B, seeds of the original
  analysis), so that each replicate is reproducible on its own. For instance,
  the following re-computes the first 10 replicates of `results/sim01_*.rds`
  (column `rep`):
  ```r
  options(ric.n_rep = 10, ric.results_dir = "results_check")
  source("R/01_sim_gaussian_fission.R")
  ```
  The real-data analyses use fixed seeds and are deterministic.
- **Parallelisation**: replicates (and gene-wise computations) are run in
  parallel on forked processes (`pbapply`; sequential on Windows), on
  `getOption("mc.cores")` cores (set in `main.R`, default: all cores but one).
  Results do not depend on the number of cores.
- **R version**: Wilcoxon tests use `wilcox_pvalue()` (`R/utils.R`), which
  reproduces `wilcox.test()` of R < 4.5 (used for the original analyses):
  changes in R 4.5 (`pnorm()` underflow) and R 4.6 (exact p-values with ties,
  asymptotic p-values < 1e-16 rounded to 0) otherwise alter Figure 1B and the
  genes selected in Appendix A7.
- **Run time and hardware**: MacBook Pro, Apple M4 Max (16 cores, 64 GB RAM),
  macOS 27.0.1, R 4.6.1: figures ~1 min; `RERUN <- TRUE` ~30 min with 12 cores
  (simulations ~5 CPU-hours: ~15 min for `R/01`, ~11 min for `R/04`, ~1 min or
  less for the other scripts).

## R session

Full output of `sessionInfo()` in `sessionInfo.txt`:

```
R version 4.6.1 (2026-06-24)
Platform: aarch64-apple-darwin23
Running under: macOS Golden Gate 27.0.1
other attached packages:
 [1] ggrepel_0.9.8      scattermore_1.2    colorspace_2.1-3   MetBrewer_0.2.0   
 [5] ggnewscale_0.5.2   tidyr_1.3.2        umap_0.2.10.0      glmpca_0.2.0      
 [9] sctransform_0.4.3  Seurat_5.5.1       SeuratObject_5.4.0 sp_2.2-3          
[13] pbapply_1.7-5      mclust_6.1.3       mvtnorm_1.4-2      npreg_1.1.1       
[17] countsplit_4.0.0   latex2exp_0.9.8    patchwork_1.3.2    ggplot2_4.0.3     
[21] dplyr_1.2.1       
```
