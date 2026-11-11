---
editor_options: 
  markdown: 
    wrap: 72
---

# Practical limitations of information-splitting approaches for post-clustering inference

Code and intermediate results reproducing all the data-based figures,
tables and numbers of *Practical limitations of information-splitting
approaches for post-clustering inference, with implications for
differential testing in single-cell RNA-seq data* (B. Hivert, D. Agniel,
R. Thiébaut, B. P. Hejblum) and of its Supplementary Materials. Zenodo:
[10.5281/zenodo.11207776](https://doi.org/10.5281/zenodo.11207776).
Licence: MIT.

## How to run

From the root folder of the archive (the working directory: all paths
are relative to it):

``` bash
Rscript main.R
```

`main.R` installs the missing R packages (all on CRAN) and runs the
scripts of `R/` in order. By default, it produces the figures
(`figures/`) and tables (`tables/`) from the intermediate results
provided in `results/` (\~1 min). With `RERUN <- TRUE`, it first
re-computes all the intermediate results (\~30 min on 12 cores, \~2 h on
4 cores), downloading the data if needed. Each script can also be run on
its own (`Rscript R/<script>.R`).

| Script | Intermediate results / outputs | Article |
|----|----|----|
| `R/01_sim_gaussian_fission.R` | `results/sim01_ideal.rds`, `results/sim01_adverse.rds` | Figures S6, S7 |
| `R/02_sim_variance_bias.R` | `results/sim02_variance_bias_S1A.rds`, `..._S1B.rds` | Figure S1 |
| `R/03_sim_nb_thinning_mixture.R` | `results/sim03_nb_thinning_mixture.rds` | Figure 1B |
| `R/04_sim_correlation_bias.R` | `results/sim04_correlation_bias.rds` | Figures 1C, S2, S3 |
| `R/05_app_overdispersion.R` | `results/app05_overdispersion.rds` | Figure S8 |
| `R/06_app_granulocytes.R` | `results/app06_granulocytes.rds` | Figure 2, Section 4.4, Figure S9 |
| `R/07_app_granulocytes_neutrophils.R` | `results/app07_granulocytes_neutrophils.rds` | Figure 3, Figures S4, S5 (Appendix A7) |
| `R/08_figures_main.R` | `figures/Figure1.pdf`, `figures/Figure2.pdf`, `tables/Numbers_in_text.csv` | Figures 1-2, numbers of Sections 2.2.1 and 4.4 |
| `R/09_figures_supplementary.R` | `figures/FigureS1.pdf` to `figures/FigureS9.pdf`, `tables/Table_A7_performances.csv` | Figures S1-S9, Appendix A7 |

`R/utils.R` contains the shared functions. Figure S10
(`figures/FigureS10.png`) is a drawing, and Tables S1-S2 are
theoretical.

## Data dictionary

### Input biological data

Tabula Sapiens bone marrow scRNA-seq data (v5.0.0, CZ CELLxGENE
Discover, CC BY 4.0), donor TSP14, from
<https://datasets.cellxgene.cziscience.com/5e736dcd-01d8-4639-805a-31fea1528be0.rds>
(277 MB, MD5 `3099d66fb374dc693a8e3486284ebcba`). This public file is
not included into the archive: `data/TS_data.rds` by `load_TSP14()`
(defined in `R/utils.R`): a Seurat object in which observations are
cells and variables are genes. Only donor TSP14 is kept.

| Field | Location in the object | Definition | Used in |
|----|----|----|----|
| `counts` | layer of the `RNA` assay | Raw UMI counts (non-negative integers), genes × cells; genes identified by their Ensembl ID | all `app` scripts |
| `feature_name` | gene metadata (`obj[["RNA"]][[]]`) | Gene symbol | `R/06`, `R/07` |
| `donor_id` | cell metadata | Donor identifier; data restricted to `"TSP14"` | all `app` scripts |
| `cell_ontology_class` | cell metadata | Cell-type annotation of the Tabula Sapiens consortium; `"granulocyte"` (454 cells) and `"cd24 neutrophil"` (2,106 cells) | `R/06`, `R/07` |
| `cell_type` | cell metadata | Standardized Cell Ontology label; `"neutrophil"`, `"macrophage"`, `"monocyte"`, `"granulocyte"`, `"CD4-positive, alpha-beta T cell"`, `"memory B cell"` | `R/05` |

### Intermediate results (in `results/`, read with `readRDS()`)

The intermediate results contain the subsets of the generated outputs
directly used in the figures.

| File | Observations (rows) | Variables (columns) |
|----|----|----|
| `sim01_ideal.rds`, `sim01_adverse.rds` | 756,000: one two-sample *t*-test per replicate × sample size × τ × fission scheme × variable | `rep`: replicate (1–1,000); `pvalues`: *t*-test p-value on X⁽²⁾ between the two estimated clusters (`NA` when a cluster has a single observation, 138 rows of the adverse file); `Fission`: `"Conditional"` or `"Marginal"`; `Variable`: `"X1"` (no mean difference between components) or `"X2"`; `ARI`: adjusted Rand index between estimated clusters and true components (`NA` in the adverse scenario); `n`: sample size (50–10,000); `tau`: fission parameter τ (0.1–5) |
| `sim02_variance_bias_S1A.rds` | 208,000: one *t*-test per replicate × σ × σ̂ | `rep`; `pvalues`; `sigma`: true standard deviation σ (0.1, 0.5, 1, 2, i.e. σ² = 0.01–4); `sigma_hat`: standard deviation used for fission (0–4); `n` = 100 |
| `sim02_variance_bias_S1B.rds` | 250,000: one *t*-test per replicate × n × σ̂ | same columns; `sigma` = 1; `sigma_hat` 0.8–1.8; `n` 50–1,000 |
| `sim03_nb_thinning_mixture.rds` | 2,000: one Wilcoxon test per replicate × thinning scheme | `rep`; `pvalues`: Wilcoxon p-value on X⁽²⁾ between the two spurious clusters; `Fission`: `"Conditional thinning"` or `"Marginal thinning"` |
| `sim04_correlation_bias.rds`, element `typeI` | 476,000: one setting per replicate × method × estimation × ρ × bias | `rep`; `Method`: `"NB"` (thinning) or `"Gauss"` (fission); `Estimation`: `"Oracle"` (true θ) or `"Wrong"` (θ̂ = θ × `Error`); `Rho`: correlation ρ (0–0.9); `Error`: ratio θ̂/θ (0.001–10); `n_reject`: number of the 50 variables with p \< 0.05; `n_tests` = 50 |
| `sim04_correlation_bias.rds`, element `pvalues` | 48,000: one Wilcoxon test (first and last variables, displayed settings) | `rep`; `pval`; `Variable`: `"Var1"` or `"Var50"`; `Method`; `Estimation`; `Error`; `Rho` |
| `app05_overdispersion.rds` | 8,333 genes | `Gene`: Ensembl ID; `neutrophil`, `macrophage`, `monocyte`, `granulocyte`, `CD4`, `memory.b`: NB overdispersion θ estimated by `sctransform::vst()` within each cell population |
| `app06_granulocytes.rds` | list (454 granulocytes × 500 most variable genes) | `counts`: UMI counts (cells × genes); `Xtrain`, `Xtest`: thinned counts X⁽¹⁾, X⁽²⁾; `cluster`: *k*-means cluster (1 or 2) of each cell; `genes` (one row per gene): `Gene`, `GeneSymbol`, `CorOriginal` (correlation with the first gene in `counts`), `CorThinning` (correlation between X⁽¹⁾ of the first gene and X⁽²⁾ of each gene), `PvalMultivariate` and `PvalUnivariate` (Wilcoxon p-values on X⁽²⁾ after multivariate or gene-wise clustering); `n_cells`: cells before / after filtering |
| `app07_granulocytes_neutrophils.rds` | list (2,560 cells × 500 genes) | `counts`: UMI counts, columns 1–250 = truly differentially expressed (H1) genes, 251–500 = genes without difference (H0), rows ordered as `true_cell`; `true_cell`: cell type of each row; `results` (2,000 rows = 500 genes × 2 thinnings × 2 clusterings): `pvalue` (Wilcoxon p-value on X⁽²⁾), `ARI` (between the clusters and `true_cell`), `Clustering` (`"Univariate"` or `"Multivariate"`), `Thinning` (`"Conditional"` or `"Marginal"`), `Gene`, `Hyp` (`"H1"` or `"H0"`), `GeneSymbol`, `TruePvalue` (Wilcoxon p-value between true cell types on the original counts; 0 below double precision); `n_cells`: cells before / after filtering |

### Output tables (`tables/`)

| File | Rows | Columns |
|----|----|----|
| `Numbers_in_text.csv` | one number quoted in Sections 2.2.1 and 4.4 | `Quantity`, `Value` |
| `Table_A7_performances.csv` | thinning × clustering (Appendix A7) | `Thinning`, `Clustering`, `TypeI` (proportion of the 250 H0 genes with p \< 0.05), `Power` (proportion of the 250 H1 genes with p \< 0.05) |
| `Figure3_numbers.tsv` | thinning × clustering (Figure 3) | `H0_significant`/`H0_total`, `H1_significant`/`H1_total`, and min, quartiles and max of the ARI over H1 genes |

## Reproducibility

- **Seeds**: replicate `r` of each simulation is run after
  `set.seed(k * 1e6 + r)` (`k` = 1, 2 for Figures S6, S7; 3, 4 for
  Figures S1A, S1B; 5 for Figures 1C, S2, S3) or `set.seed(310123 * r)`
  (Figure 1B, seeds of the original analysis), so that each replicate is
  reproducible on its own. For instance, the following re-computes the
  first 10 replicates of `results/sim01_*.rds` (column `rep`):

  ``` r
  options(ric.n_rep = 10, ric.results_dir = "results_check")
  source("R/01_sim_gaussian_fission.R")
  ```

  The real-data analyses use fixed seeds and are deterministic.

- **Parallelisation**: replicates (and gene-wise computations) are run
  in parallel on forked processes (`pbapply`; sequential on Windows), on
  `getOption("mc.cores")` cores (set in `main.R`, default: all cores but
  one). Results do not depend on the number of cores.

- **R version**: Wilcoxon tests use `wilcox_pvalue()` (`R/utils.R`),
  which reproduces `wilcox.test()` of R \< 4.5 (used for the original
  analyses): changes in R 4.5 (`pnorm()` underflow) and R 4.6 (exact
  p-values with ties, asymptotic p-values \< 1e-16 rounded to 0)
  otherwise alter Figure 1B and the genes selected in Appendix A7.

- **Run time and hardware**: MacBook Pro, Apple M4 Max (16 cores, 64 GB
  RAM), macOS 27.0.1, R 4.6.1: figures \~1 min; `RERUN <- TRUE` \~30 min
  with 12 cores (simulations \~5 CPU-hours: \~15 min for `R/01`, \~11
  min for `R/04`, \~1 min or less for the other scripts).

## R session

Full output of `sessionInfo()` in `sessionInfo.txt`:

```         
R version 4.6.0 (2026-04-24)
Platform: aarch64-apple-darwin23
Running under: macOS Sonoma 14.6.1
other attached packages:
 [1] ggrepel_0.9.8      scattermore_1.2    colorspace_2.1-3   MetBrewer_0.2.0   
 [5] ggnewscale_0.5.2   tidyr_1.3.2        umap_0.2.10.0      glmpca_0.2.0      
 [9] sctransform_0.4.3  Seurat_5.5.0       SeuratObject_5.4.0 sp_2.2-1          
[13] pbapply_1.7-4      mclust_6.1.2       mvtnorm_1.3-7      npreg_1.1.1       
[17] countsplit_4.0.0   latex2exp_0.9.8    patchwork_1.3.2    ggplot2_4.0.3     
[21] dplyr_1.2.1       
```
