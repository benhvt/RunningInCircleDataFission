################################################################################
#                                                                              #  
#                         Main Make Supplementary Figures                      #      
#                         Benjamin Hivert - 23/01/2026                         #
#                                                                              #
################################################################################

# This files contains codes to generate all the main figures of the paper. 
# All results were pre-imported and prepared in 20250417_PrepareResults.R script

# Library 

library(dplyr)
library(ggplot2)
library(patchwork)
library(latex2exp)
library(mvtnorm)
library(readxl)
library(umap)
library(glmpca)

theme_set(theme_classic())

#------------------------------------------------------------------------------#
#                                   Supp Section 4                             #                                 
#               Impact of variance estimation on Type I error rate             # 
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#                           Internals functions                                #
#------------------------------------------------------------------------------#

# Marginal variance of a Gaussian Mixture
mixture_variance <- function(pi, mu_list, sigma_list) {
  # Function that computes the marginal variance of a Gaussian mixtures models of
  # parameters :
  #pi: a vector of length G that contains mixing proportion 
  # mu_list: a list of length G that contains vector of length p that are the means of each components
  #sigma_list: a list of G that contains matrix of dimension p x p that are the component-specific covariance
  #
  #return a matrix of dimension p x p that contains the marginal covariance of the mixture 
  
  G <- length(pi)  #Number of component in the mixture
  d <- length(mu_list[[1]])  # Dimension
  mu_Z <- Reduce("+", Map("*", pi, mu_list))  # Marginal mean of the mixture
  
  # Initialisation of the covariance matrix
  Sigma <- matrix(0, nrow = d, ncol = d)
  
  # Marginal variance computation
  for (g in 1:G) {
    Sigma <- Sigma + pi[g] * (sigma_list[[g]] + tcrossprod(mu_list[[g]]))
  }
  Sigma <- Sigma - tcrossprod(mu_Z)
  
  return(Sigma)
}

# Conditional covariance in marginal fission 
margFission_condCovariance <- function(g,
                                       pi,
                                       list_mu,
                                       list_Sigma){
  # Function that compute the conditional covariance between X1 and X2 when marginal
  # data fission is applied. This covariance is equale to Sigma_g - Sigma
  #
  # parameters :
  #g: The number of the component of the mixture where the conditional covariance must be computed
  #pi: a vector of length G that contains mixing proportion 
  # mu_list: a list of length G that contains vector of length p that are the means of each components
  #sigma_list: a list of G that contains matrix of dimension p x p that are the component-specific covariance
  #
  #return a matrix of dimension p x p that contains the conditional covariance 
  Sigma <- mixture_variance(pi, list_mu, list_Sigma)
  margFission <- Sigma - list_Sigma[[g]]
}


# Marginal covariance in conditional fission
condFission_margCovariance <- function(pi,
                                       list_mu,
                                       list_Sigma){
  # Function that compute the marginal covariance between X1 and X2 when conditional
  # data fission is applied.
  #
  # parameters :
  #pi: a vector of length G that contains mixing proportion 
  # mu_list: a list of length G that contains vector of length p that are the means of each components
  #sigma_list: a list of G that contains matrix of dimension p x p that are the component-specific covariance
  #
  #return a matrix of dimension p x p that contains the marginal covariance 
  
  # Marginal covariance in conditional fission
  list_pi <- list(pi[1], pi[2])
  margMu <- Reduce(`+`, Map(`*`,list_pi, list_mu))
  
  margSigma <- Reduce(`+`, Map(function(pi_k, mu_k) {
    diff <- matrix(mu_k - margMu, ncol = 1)
    pi_k * (diff %*% t(diff))  # Produit matriciel
  }, list_pi, list_mu))
  
}


compute_typeI <- function(n, tau, sigma, sigma_hat, alpha){
  # Function that computes theoritical type I error of data fission
  qu_t <- qt(alpha/2, df = n-2, lower.tail = FALSE)
  qu_n <- qnorm(alpha/2, lower.tail = F)
  cor_fg <- (sigma^2-sigma_hat^2)/sqrt((sigma^2 + (tau^2)*sigma_hat^2)*(sigma^2 + (1/tau^2)*sigma_hat^2))  
  mean_t <- sqrt(n)*sqrt((2/pi)*cor_fg^2)/sqrt(1-(2/pi)*cor_fg^2)
  
  return(pnorm(qu_n, mean = mean_t, sd= 1, lower.tail = FALSE) + pnorm(-qu_n, mean = mean_t, sd= 1, lower.tail = TRUE))
}


#------------------------------------------------------------------------------#
#                           Simulations setting                                #
rm("pi")
n <- 100
n_grid <- c(50, 100, 200, 500, 1000)
sigma <- c(0.1, 0.5, 1, 2)
sigma_grid <- sort(c(seq(0, 4, length.out = 50), 2, 0.1))
sigma_grid2 <- seq(0.8, 1.8, length = 50)
tau <- .4

#------------------------------------------------------------------------------#
#             First panel: Function of the original variance                   #


sigma_grid_plot <- seq(0, 6*max(sigma), length.out = 1000)

typeI_theo <- lapply(sigma_grid_plot, function(s){
  temp <- compute_typeI(n, tau , sigma, sigma_hat = s, alpha = .05)
  return(data.frame(TypeItheo = temp,
                    sigma_hat = s,
                    sigma = sigma))
})

typeI_df <- do.call(rbind.data.frame, typeI_theo) %>% mutate(sigma_name = paste0("sigma^2==", sigma^2))


var_bias <- read.csv(file = "results/VarianceEstimationAndTypeIError.csv")

var_bias_res <- var_bias %>%
  group_by(sigma, sigma_hat, n) %>%
  summarise(EmpTypeI = mean(pvalues < 0.05)) %>%
  ungroup() %>%
  # mutate(TheTypeI = compute_typeI(n = n, tau = tau, sigma = sigma, sigma_hat = sigma_hat, alpha = .05)) %>%
  mutate(sigma_name = paste0("sigma^2==", sigma^2)) %>%
  mutate(Ratio = (sigma^2-sigma_hat^2)/sigma^2) %>%
  reshape2::melt(id.vars = c("Ratio", "sigma_name", "sigma", "sigma_hat"), measure_vars = c("TheTypeI", "EmpTypeI"))

plot_sigma <-  ggplot(var_bias_res) + aes(x=Ratio, y = value, colour = sigma_name) +
  geom_point(data = subset(var_bias_res, variable == "EmpTypeI"), aes(shape = "Empirical"), size = 4) +
  scale_shape_manual(name = "Empirical", values = 2, labels = '', guide = "legend") +
  # geom_line(data = subset(df_res_sigma, variable == "TheTypeI"), aes(group = sigma_name, linetype = "Theoritical"), size = 1) +
  geom_line(data = typeI_df, aes(x= (sigma^2-sigma_hat^2)/sigma^2, y = TypeItheo, group = sigma_name, colour = sigma_name, linetype = "Theoritical"), size = 1) +
  scale_linetype_manual(name = "Theoritical", values = 1, labels = "", guide = "legend") +
  scale_colour_manual(name = '',
                      values = c("#93B5C6", "#DBC2CF", "#998bc0", "#BD4F6C"),
                      labels = c(TeX(r'($\sigma^2 = 0.01$)'),
                                 TeX(r'($\sigma^2 = 0.25$)'),
                                 TeX(r'($\sigma^2 = 1$)'),
                                 TeX(r'($\sigma^2 = 4$)'))) +
  xlab(TeX(r'($(\sigma^2 - \widehat{\sigma^2})/\sigma^2$)')) +
  xlim(c(-5,2)) +
  ylab("Type I error rate") +
  ggnewscale::new_scale_colour() +
  geom_hline(aes(yintercept = 0.05, 
                 colour = "5% nominal levels"),
             linetype = 2,
             size = 1.2) +
  scale_colour_manual(name = "",
                      values = "#6C0E23") +
  guides(
    shape = guide_legend(order = 1),
    linetype = guide_legend(order = 2),
    color = guide_legend(order = 3)
  ) +
  NULL


#------------------------------------------------------------------------------#
#               Second panel: Function of the sample size                      #


var_bias_sampsize <- read.csv(file = "results/VarianceEstimationAndTypeIErrorSampleSize.csv")

var_bias_sampsize_res <-  var_bias_sampsize %>%
  group_by(n, sigma_hat, sigma) %>%
  rename(SampSize = n) %>%
  summarise(EmpTypeI = mean(pvalues < 0.05)) %>%
  ungroup() %>%
  mutate(TheTypeI = compute_typeI(n = SampSize, 
                                  tau = tau, 
                                  sigma = sigma,
                                  sigma_hat = sigma_hat, 
                                  alpha = .05)) %>%
  mutate(sigma_name = paste0("sigma^2==", sigma^2)) %>%
  mutate(Ratio = (sigma^2-sigma_hat^2)/sigma^2) %>%
  reshape2::melt(id.vars = c("Ratio", "sigma_name", "SampSize", "sigma", "sigma_hat"), measure_vars = c("TheTypeI", "EmpTypeI"))

plot_sigma_n <- ggplot(var_bias_sampsize_res) + aes(x=Ratio, y = value, colour = as.factor(SampSize)) +
  geom_point(data = subset(var_bias_sampsize_res, variable == "EmpTypeI"),
             aes(shape = "Empirical"), 
             size = 4) +
  scale_shape_manual(name = "Empirical",
                     values = 2, labels = '',
                     guide = "legend") +
  geom_line(data = subset(var_bias_sampsize_res, variable == "TheTypeI"), 
            aes(group = SampSize, linetype = "Theoritical"), size = 1) +
  scale_linetype_manual(name = "Theoritical", 
                        values = 1,
                        labels = "", 
                        guide = "legend") +
  scale_colour_manual(name = "", 
                      values = colorRampPalette(c("#008154", "#0092a4", "#2a2956"))(length(n_grid)),
                      labels = paste0("n=", n_grid)) +
  xlab(TeX(r'($(\sigma^2 - \widehat{\sigma^2})/\sigma^2$)')) +
  ylab("Type I error rate") +
  ggnewscale::new_scale_colour() +
  geom_hline(aes(yintercept = 0.05, 
                 colour = "5% nominal levels"),
             linetype = 2,
             size = 1.2) +
  scale_colour_manual(name = "",
                      values = "#6C0E23") +
  guides(
    shape = guide_legend(order = 1),
    linetype = guide_legend(order = 2),
    color = guide_legend(order = 3)
  ) +
  NULL

(plot_sigma + plot_sigma_n)  +
  plot_annotation(tag_levels = "A") &
  theme_classic() &
  theme(axis.title = element_text(size = 24), 
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        plot.tag = element_text(face = "bold", size = 24),
        legend.spacing = unit(0.01, "cm"))

ggsave(filename = "Supplementary Figures/SuppFigureSection4.pdf",
       width = 350, 
       height = 120, 
       units = "mm",
       dpi = 600)


#------------------------------------------------------------------------------#
#                              Supp Section 5                                  #
#               Applications on real world data  sets of a                     #  
#                         two cell populations mixture                         #
#------------------------------------------------------------------------------#
H1_genes <- H0_genes <- 250

results_application <- read.csv("results/SupplementaryResultsApplicationOnScRNASeqWithNull.csv") %>%
  mutate(ClusteringName = paste(Clustering, "Clustering", sep = " "),
         ThinningName = paste(Thinning, "Thinning", sep = " ")) %>%
  mutate(Hyp = rep(c(rep("H1",H1_genes),
                     rep("H0",H0_genes)),4))

tagName <- results_application %>% 
  filter(Hyp == "H1") %>%
  group_by(ThinningName, ClusteringName) %>%
  arrange(pvalue) %>%
  slice_head(n = 5)

plt_scatter_pval_H0 <-results_application  %>%
  filter(Hyp == "H0") %>%
  ggplot() +
  aes(x=-log10(TruePvalue),
      y=-log10(pvalue), 
      colour = ARI) +
  geom_point(alpha = .5,
             size = 2) +
  # ggrepel::geom_label_repel(data = tagName, aes(label = GeneSymbol)) +
  geom_hline(yintercept = -log10(0.05)) +
  # geom_vline(xintercept = (0.05)) +
  facet_grid(ClusteringName~ThinningName, scales = "free") +
  xlab("-log10(True p-value)") +
  ylab("-log10(p-value)") +
  scale_colour_viridis_c(option = "plasma",
                         direction = -1, limits = c(min(results_application$ARI),1)) +
  ggtitle(TeX(r'(Distribution of p-values for genes under $H_0$)'))

plt_scatter_pval_H1 <-results_application  %>%
  filter(Hyp == "H1") %>%
  ggplot() +
  aes(x=-log10(TruePvalue),
      y=-log10(pvalue), 
      colour = ARI) +
  geom_point(alpha = .5,
             size = 2) +
  ggrepel::geom_label_repel(data = tagName, aes(label = GeneSymbol)) +
  geom_hline(yintercept = -log10(0.05)) +
  # geom_vline(xintercept = (0.05)) +
  facet_grid(ClusteringName~ThinningName, scales = "free") +
  xlab("-log10(True p-value)") +
  ylab("-log10(p-value)") +
  scale_colour_viridis_c(option = "plasma",
                         direction = -1, limits = c(min(results_application$ARI),1)) +
  ggtitle(TeX(r'(Distribution of p-values for genes under $H_1$)'))


suppFigureSec5 <- plt_scatter_pval_H0 / plt_scatter_pval_H1 + 
  plot_layout(guides = "collect") &
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 14)) 

ggsave(suppFigureSec5,
       filename = "Supplementary Figures/SuppFigureSection5.pdf",
       width = 150, 
       height = 250, 
       units = "mm")

indicator_performances <- results_application %>% 
  group_by(Thinning, Clustering) %>%
  summarise(Power = mean(pvalue < 0.05 & TruePvalue<0.05), 
            typeI = mean(pvalue < 0.05 & TruePvalue > 0.05))


#------------------------------------------------------------------------------#
#                 Sup. Figure 1: QQ-Plot of pvalues when 
#                         correlated data generation            
#------------------------------------------------------------------------------#

neg_bin_cor_res <- read.csv("results/CorrelationAndRelativeBiais.csv")

qqplot_cor <- neg_bin_cor_res %>% 
  filter(Error == 1) %>%
  filter(Variable %in% c("Var1", "Var50")) %>%
  filter(Estimation == "Oracle") %>%
  mutate(VarName = ifelse(Variable == "Var1", "X[1]", "X[50]")) %>%
  mutate(MethodName = ifelse(Method == "NB", "Negative~Binomial~Thinning", "Gaussian~Fisson")) %>%
  ggplot() +
  geom_abline(slope=1, intercept=0, col="red", size = 1.2, alpha = .7) + xlab("Theoretical Quantiles") + 
  stat_qq(aes(sample = pval, 
              colour = factor(Rho)),
          distribution = qunif, 
          size = 1,
          alpha = .5) +
  facet_grid(VarName~MethodName, labeller = label_parsed) +
  scale_colour_manual(name = TeX(r'(Correlation ($\rho$))'),
                      values = c(
                        "#FFE5B4",  # abricot clair
                        "#FFBFA0",  # corail doux
                        "#FF8CA0",  # rose saumon
                        "#C478B8",  # violet chaud
                        "#8050A0",  # violet saturé
                        "#50286F",  # prune foncé
                        "#2B1240"   # aubergine très foncé
                      )) +
  ylab("Empirical Quantiles") + 
  xlim(c(0, 1)) + ylim(c(0, 1)) + 
  theme_classic() +
  theme(text = element_text(size = 14))

ggsave(filename = "Supplementary Figures/SuppFigure1.pdf", 
       plot = qqplot_cor,
       width = 200,
       height = 150,
       units = "mm")

#------------------------------------------------------------------------------#
#                 Sup. Figure 2: QQ-Plot of pvalues when 
#                   correlated data generation as a
#                       function of relative biais 
#------------------------------------------------------------------------------#

neg_bin_cor_res <- read.csv("results/CorrelationAndRelativeBiais.csv")

qqplot_cor_biais <- neg_bin_cor_res %>% 
  mutate(Theta_hat = Theta*Error) %>% 
  filter(Rho == 0) %>%
  filter(Variable %in% c("Var1", "Var50")) %>%
  filter(Estimation == "Wrong") %>%
  filter(Error %in% c(0.001, 0.5, 1, 1.5, 5)) %>%
  group_by(Error) %>%
  mutate(RelativeBiais = (Theta_hat-Theta)/Theta) %>% 
  mutate(VarName = ifelse(Variable == "Var1", "X[1]", "X[50]")) %>%
  mutate(MethodName = ifelse(Method == "NB", "Negative~Binomial~Thinning", "Gaussian~Fisson")) %>%
  ggplot() +
  geom_abline(slope=1, 
              intercept=0, 
              col="red", 
              size = 1.2, 
              alpha = .7) + xlab("Theoretical Quantiles") + 
  stat_qq(aes(sample = pval, 
              colour = factor(RelativeBiais)),
          distribution = qunif, 
          size = 1, 
          alpha = .5) +
  facet_grid(VarName~MethodName, labeller = label_parsed) +
  MetBrewer::scale_color_met_d(name = "Hiroshige") +
  guides(colour = guide_legend("Relative Biais")) +
  # scale_colour_manual(name = TeX(r'(Correlation ($\rho$))'),
  #                     values = c(
  #                       "#FFE5B4",  # abricot clair
  #                       "#FFBFA0",  # corail doux
  #                       "#FF8CA0",  # rose saumon
  #                       "#C478B8",  # violet chaud
  #                       "#8050A0",  # violet saturé
  #                       "#50286F",  # prune foncé
  #                       "#2B1240"   # aubergine très foncé
  #                     )) +
  ylab("Empirical Quantiles") + 
  xlim(c(0, 1)) + ylim(c(0, 1)) + 
  theme_classic() +
  theme(text = element_text(size = 14))

ggsave(filename = "Supplementary Figures/SuppFigure2.pdf", 
       plot = qqplot_cor_biais,
       width = 200,
       height = 150,
       units = "mm")
#------------------------------------------------------------------------------#
#                               Sup. Figure 3                                  #
#                         Applications on Bonne marrow                         #
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#                     Gene-specific overdispersion parameters                  #

# Parameter
cell_pop_to_test <- c("neutrophil", 
                      "macrophage",
                      "monocyte",
                      "granulocyte",
                      "CD4-positive, alpha-beta T cell",
                      "memory B cell")

cell_theta <- read.csv(file = "results/Application_CellPopulationOverdispersion.csv")

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

pair_cellPop <- combn(2:(length(cell_pop_to_test)+1),2)
allPairplot <- lapply(1:ncol(pair_cellPop), function(p){
  if (colnames(cell_theta)[pair_cellPop[1,p, drop = T]] == "memory.b"){
    lab_x <- "memory b cells"
  }
  else{
    lab_x <- firstup(colnames(cell_theta)[pair_cellPop[1,p, drop = T]])
  }
  
  if (colnames(cell_theta)[pair_cellPop[2,p, drop = T]] == "memory.b"){
    lab_y <- "memory b cells"
  }
  else{
    lab_y <- firstup(colnames(cell_theta)[pair_cellPop[2,p, drop = T]])
  }
  df_temp <- data.frame(Gene = cell_theta$Gene, 
                        Pop1 = cell_theta[,pair_cellPop[1,p, drop = T]],
                        Pop2 = cell_theta[,pair_cellPop[2,p, drop = T]])
  
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
    scale_x_log10(breaks = c(0.01, 0.1, 1, 10), 
                  labels = c(0.01, 0.1, 1, 10)) +
    scale_y_log10() +
    annotation_logticks(side = "bl") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    NULL
})

# plt_application <- ((allPairplot[[1]]) + plot_spacer()  + allPairplot[[2]] + plot_spacer()  + allPairplot[[3]] + plot_spacer()  + allPairplot[[4]]) + plot_layout(nrow = 1, widths = c(4,.5,4,.5,4,.5,4))
plt_application <- ((allPairplot[[1]]) + allPairplot[[2]] + allPairplot[[3]] + allPairplot[[4]]) +
  plot_layout(nrow = 1) & 
  theme_classic() +
  theme(text = element_text(size = 14))

ggsave(plt_application,
       filename = "Supplementary Figures/SuppFigure3.pdf",
       width = 200,
       height = 75,
       units = "mm")


#------------------------------------------------------------------------------#
#                               Sup. Figure 4                                  #
#                 Visualisation of the homogenous population                   #
#------------------------------------------------------------------------------#

TS_granulo <- read_xlsx("results/ResultsOnScRNASeq.xlsx", 
                        sheet = "RawCountsFilter")

# GLM-PCA with L = 20 latent factors to estimate 
set.seed(20260311)
glmpca_granulo <- glmpca(Y = t(TS_granulo), 
                    L = 50, 
                    fam = "nb")

glmPCplot <- ggplot(glmpca_granulo$factors) +
  aes(x=dim1, 
      y = dim2,
      colour = "granulocyte") +
  geom_hline(yintercept = 0,
             colour = "black",
             linetype = "dashed") +
  geom_vline(xintercept = 0,
             colour = "black",
             linetype = "dashed") +
  geom_point(size = 1.5,
             alpha = .8) +
  scale_colour_manual(name = "True cells",
                      values = c("granulocyte"= "#294122",
                                 "cd24 neutrophil" = "#EB3D00")) +
  xlab("GLM-PC 1") +
  ylab("GLM-PC 2")

# UMAP on the L=20 GLM-PCs of the log2CPM 
## UMAP config
config <- umap::umap.defaults
config$n_neighbors <- 80
config$min_dist <- 0.6

## Run UMAP
umap_TS_granulo <- umap(
  glmpca_granulo$factors,
  config = config
)

## UMAP plot
umapPlot <- data.frame(umap_TS_granulo$layout) %>%
  ggplot() +
  aes(x = X1, 
      y = X2,,
      colour = "granulocyte") +
  geom_hline(yintercept = 0,
             colour = "black",
             linetype = "dashed") +
  geom_vline(xintercept = 0,
             colour = "black",
             linetype = "dashed") +
  scale_colour_manual(name = "True cells",
                      values = c("granulocyte"= "#294122",
                                "cd24 neutrophil" = "#EB3D00")) +
  geom_point(size = 1.5,
             alpha = .8) +
  xlab("UMAP 1") +
  ylab("UMAP 2") 

H0cellsPlot <- glmPCplot + umapPlot +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        text = element_text(size = 14))

ggsave(H0cellsPlot, 
       filename = "Supplementary Figures/SuppFigure4.pdf",
       width = 150, 
       height = 75,
       units = "mm")

#------------------------------------------------------------------------------#
#                        Sup. Figure 2 of Section5                             #
#              Visualisation of the 2 cells populations data                   #
#------------------------------------------------------------------------------#

TS_granulo_neutro <- read.csv("results/ResultsOnScRNASeq2CellPopulations.csv")
counts_granulo_neutro <- TS_granulo_neutro %>% dplyr::select(-c(TrueCell))
# GLM-PCA with L = 20 latent factors to estimate 
set.seed(20260311)

glmpca_granulo_neutro <- glmpca(Y = t(counts_granulo_neutro), 
                         L = 50, 
                         fam = "nb")

glmPCplot_H1 <- data.frame(glmpca_granulo_neutro$factors,
                           TrueCell = TS_granulo_neutro$TrueCell) %>%
  ggplot() +
  aes(x=dim1, y = dim2,
      colour = TrueCell) +
  geom_hline(yintercept = 0,
             colour = "black",
             linetype = "dashed") +
  geom_vline(xintercept = 0,
             colour = "black",
             linetype = "dashed") +
  geom_point(size = 1,
             alpha = .6) +
  scale_colour_manual(name = "True cells",
                      values = c("granulocyte"= "#294122",
                                 "cd24 neutrophil" = "#EB3D00")) +
  xlab("GLM-PC 1") +
  ylab("GLM-PC 2")

# UMAP on the L=20 GLM-PCs of the log2CPM 
## UMAP config
config <- umap::umap.defaults
config$n_neighbors <- 80
config$min_dist <- 0.6

## Run UMAP
umap_TS_granulo_neutro <- umap(
  glmpca_granulo_neutro$factors,
  config = config
)

## UMAP plot
umapPlot_H1 <- data.frame(umap_TS_granulo_neutro$layout,
                          TrueCell = TS_granulo_neutro$TrueCell) %>%
  ggplot() +
  aes(x = X1, 
      y = X2,
      colour = TrueCell) +
  geom_hline(yintercept = 0,
             colour = "black",
             linetype = "dashed") +
  geom_vline(xintercept = 0,
             colour = "black",
             linetype = "dashed") +
  geom_point(size = 1,
             alpha = .6) +
  scale_colour_manual(name = "True cells",
                      values = c("granulocyte"= "#294122",
                                 "cd24 neutrophil" = "#EB3D00")) +
  xlab("UMAP 1") +
  ylab("UMAP 2")

H1cellsPlot <- glmPCplot_H1 + umapPlot_H1 +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        text = element_text(size = 14))

ggsave(H1cellsPlot, 
       filename = "Supplementary Figures/SuppFigureSection52.pdf",
       width = 150, 
       height = 75,
       units = "mm")


