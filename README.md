# Running in circles: practical limitations for real-life application of data fission and data thinning in post-clustering differential analysis

This repository contains code to reproduce simulations and Figure of the article *Running in circles: practical limitations for real-life application of data fission and data thinning in post-clustering differential analysis*.

## Organisation of the repos

```graphql
.
└── RevisionRunningInCircles/
    ├── codes      
    ├── data
    ├── raw_results
    ├── results
    ├── Figures
    ├── LICENSE
    ├── README.md
    └── RunningInCircleV2.Rproj
```

## Descriptions

### codes

This folder contains all the codes needed to reproduce simulations, applications and figures.

### data

This folder should contain the actual scRNA-seq data from**Tabula Sapiens**(in`.rds` format), which are used by the script `codes/20250417_BoneMarrowApplication.R`.

→ These data could be downloaded [here](https://datasets.cellxgene.cziscience.com/5e736dcd-01d8-4639-805a-31fea1528be0.rds).

### results

This folder contains all the pre-processed data, including results from simulations and scRNA-seq applications. The majority of these files were generated using the R script `20250417_PrepareResults.R`.

To optimize computation time, most of the simulations were run on [CURTA](https://www.mcia.fr). The raw simulation outputs were stored in the `raw_results` folder. The R script `20250417_PrepareResults.R` was then used to aggregate these outputs and prepare the final datasets for figure generation. No further processing or derivation of the simulation results was performed—only aggregation across the parallel computations.
**WARNNING: Some of the `.csv` result files have been compressed and must be unzipped before running `codes/MakeFigures.R`.**

### Figures

This folder contains all the figures generated with `R` script `codes/MakeFigures.R` .
