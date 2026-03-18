################################################################################
#                                                                              #  
#                         Main Make Figures Generations                        #      
#                         Benjamin Hivert - 17/04/2025                         #
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
library(scales)
theme_set(theme_classic())


#------------------------------------------------------------------------------#
#                                  Figures                                     #
#               Behaviour of marginal and conditional data fission             # 
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
#                           Simulations setting                                #


# Parameters of the mixtures (2 components)
pi <- c(0.5, 0.5)
mu <- list(c(0,5),
           c(0,0))
Sigma <- list(diag(1,2),
              cbind(c(1, 0.5),
                    c(0.5, 1)))

# Generation of the illustrative data example 
set.seed(20250404)

G <- length(pi)
sample_size <- c(50, 100, 250, 500, 1000, 5000, 10000)

Z <- sample(1:G, sample_size[4], replace = TRUE, prob = pi)
X <- matrix(NA, nrow = sample_size[4], ncol = 2)
for (g in 1:G){
  X[Z==g,] <- rmvnorm(n = sum(Z==g), 
                      mean = mu[[g]], 
                      sigma = Sigma[[g]])
}

#------------------------------------------------------------------------------#
#                        First panel of the Figure                             #
#                  Statistical power in the ideal scenario                     #


# Figure generation 
plt_illu_power <- data.frame(X1 = X[,1],
                             X2 = X[,2],
                             Cluster = as.factor(kmeans(X, centers = 2, nstart = 100)$cluster)) %>%
  ggplot() +
  aes(x=X1, 
      y = X2, 
      colour = Cluster) +
  geom_point(size = 2,
             alpha = .8) +
  scale_colour_manual(name = "Estimated clusters",
                      values = c("#294122", "#EB3D00")) +
  theme_classic() + 
  theme(legend.position = "top") +
  xlab(TeX(r'($X_1$)')) +
  ylab(TeX(r'($X_2$)'))

# Importation of results
ideal_sc <- read.csv(file = "results/IdealScenario_GaussianFission.csv")

# Derivation: compute Statistical power and ARI

powerResults <- ideal_sc %>%
  mutate(Fission_lab = paste(Fission, "fission", sep = " ")) %>%
  group_by(Fission_lab, Variable, tau, n) %>%
  summarise(Power = mean(pvalues < 0.05)) %>%
  mutate(Variable_lab = ifelse(Variable == "X1", "X[1]", "X[2]")) 

ariResults <- ideal_sc %>%
  mutate(Fission_lab = paste(Fission, "fission", sep = " ")) %>%
  group_by(Fission_lab, tau, n) %>%
  summarise(ARI_m = mean(ARI),
            ARI_sd = sd(ARI))

plt_power <- powerResults %>%
  mutate(Fission_lab = gsub(" ", "~", Fission_lab)) %>%
  ggplot() +
  aes(x=tau, y = Power, 
      colour = as.factor(n)) +
  geom_line(linewidth = .75,
            alpha = .8) +
  scale_colour_manual(name = "Sample Size",
                      values = MetBrewer::met.brewer("Derain", n=7)) +
  scale_linetype_manual(name = "Data fission",
                        values = c(1, 6)) +
  guides(colour = guide_legend(nrow = 2)) +
  ggnewscale::new_scale_colour() +
  geom_hline(aes(yintercept = 0.05, 
                 colour = "5% nominal level"),
             linewidth = .75, 
             linetype = "dashed",
             show.legend = FALSE) +
  scale_colour_manual(name = "",
                      values = "#DB2763") +
  facet_grid(Variable_lab~Fission_lab, 
             labeller = label_parsed) +
  xlab(TeX(r'( Tunning parameter $\tau$)')) +
  ylab("Statistical Power \n 5% level") +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal"
  ) 

plt_ari <- ggplot(ariResults) +
  aes(x=tau, y = ARI_m, 
      colour = as.factor(n)) +
  geom_line(linewidth = .75, 
            alpha = .8) +
  facet_grid(~Fission_lab) +
  scale_colour_manual(name = "Sample Size",
                      values = MetBrewer::met.brewer("Derain", n=7)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  
  xlab(TeX(r'( Tunning parameter $\tau$)')) +
  ylab("Adjusted Rand \n Index") +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal"
  ) +
  guides(color = guide_legend(nrow = 2))

figure1 <- ((plt_illu_power + theme(legend.position = "bottom")) +
               (plt_ari/plt_power + plot_layout(guides = "collect",
                                                heights = c(1,4)))) +
  plot_layout(widths = c(.8,1)) +
  plot_annotation(tag_levels = list(c("A", "B", "C"))) &
  theme(legend.position = "bottom",
        text = element_text(size = 12)) &
  theme(plot.tag = element_text(face = "bold"))

ggsave(figure1,
       filename = "Figures/figure1.pdf",
       width = 190,
       height = 120,
       units = "mm")

#------------------------------------------------------------------------------#
#                               Figure 2                                       #
#               Type I error rate in the adverse scenario                      #

# Clustering: Only the first component of the mixture is clusterized into two spurious clusters
clust_comp <- kmeans(X[Z==1,], centers = 2, nstart = 100)$cluster
cluster <- rep(NA, sample_size[4])
cluster[Z==1] <- clust_comp
cluster[Z==2] <- 3

# Figure generation 
plt_illu_typeI <- data.frame(X1 = X[,1],
                             X2 = X[,2],
                             Cluster = as.factor(cluster)) %>%
  ggplot() +
  aes(x=X1, 
      y = X2, 
      colour = Cluster) +
  geom_point(size = 2,
             alpha = .8) +
  scale_colour_manual(name = "Estimated clusters",
                      values = c("#294122", "#EB3D00", "#FFBBA6")) +
  theme_classic() + 
  theme(legend.position = "bottom") +
  xlab(TeX(r'($X_1$)')) +
  ylab(TeX(r'($X_2$)'))

# Import files of results
adverse_sc <- read.csv(file = "results/AdverseScenario_GaussianFission.csv")

plt_typeI <- adverse_sc %>%
  mutate(Fission_lab = paste(Fission, "fission", sep = "~")) %>%
  group_by(Fission_lab,Variable, tau, n) %>%
  summarise(TypeI = mean(pvalues < 0.05)) %>%
  mutate(Variable_lab = ifelse(Variable == "X1", "X[1]", "X[2]")) %>%
  ggplot() +
  aes(x=tau, y = TypeI, 
      colour = as.factor(n)) +
  geom_line(linewidth = .75, 
            alpha = .8) +
  scale_colour_manual(name = "Sample Size",
                      values = MetBrewer::met.brewer("Derain", n=7)) +
  guides(color = guide_legend(nrow = 2)) +
  ggnewscale::new_scale_colour() +
  geom_hline(aes(yintercept = 0.05, 
                 colour = "5% nominal level"),
             linewidth = .75, 
             linetype = "dashed",
             show.legend = FALSE) +
  scale_colour_manual(name = "",
                      values = "#DB2763") +
  facet_grid(Variable_lab~Fission_lab, 
             labeller = label_parsed) +
  xlab(TeX(r'( Tunning parameter $\tau$)')) +
  ylab("Type I error rate \n 5% nominal levels")  +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal"
  ) 


figure2 <- plt_illu_typeI + plt_typeI +
  plot_layout(widths = c(.8,1)) +
  plot_annotation(tag_levels = list(c("A", "B"))) &
  theme(text = element_text(size = 12)) &
  theme(plot.tag = element_text(face = "bold"))

ggsave(figure2,
       filename = "Figures/figure2.pdf",
       width = 190,
       height = 110,
       units = "mm")


#------------------------------------------------------------------------------#
#                                  Figure 3                                    #
#         Data fission for scRNA-seq (Binomial Negatives Simulations)          # 
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
#           Panel A and B: Simulations results under Mixture Settings          #

#-- Parameters 
n <- 100
Z <- rep(1:2, each = n/2)
probs <- c(0.5, 0.4)
size <- c(5, 40)
set.seed(09022024)

x1 <- rnbinom(n=n/2, prob = probs[1], size = size[1])
x2 <- rnbinom(n=n/2, prob = probs[2], size = size[2])
y <- rnbinom(n=n, prob = .5, size = 5)
X <- cbind.data.frame(X1=c(x1,x2), 
                      X2=y, 
                      TrueClasses=as.factor(Z))
cl_X <- kmeans(cbind(c(x1,x2), y), centers = 3)$cluster

plt_illu_negbin <- ggplot(X) + aes(x=X1, 
                                   y=X2, 
                                   colour = as.factor(cl_X)) + #, 
                                   #shape = as.factor(Z)) + 
  geom_point(size = 1.5) +
  # scale_shape_manual(name = "True classes",
  #                    values = c(15,16)) +
  scale_colour_manual(name = "Clusters", 
                      values = c("#294122", "#EB3D00", "#FFBBA6"),
                      labels = c(TeX(r'($C_1$)'),
                                 TeX(r'($C_2$)'),
                                 TeX(r'($C_3)'))) +
  xlab(TeX(r'($X_1$)')) +
  ylab(TeX(r'($X_2$)')) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.box = "vertical") +
  guides(color = guide_legend(nrow = 1))
  # theme(legend.position = "bottom") +
  NULL

neg_bin_typeI <- read.csv("results/TypeIThinningNegBin.csv")

plt_typeI_negbin <- neg_bin_typeI %>%
  ggplot() + 
  geom_abline(slope=1, intercept=0, col="red", size = 1.2, alpha = .7) + xlab("Theoretical Quantiles") + 
  stat_qq(aes(sample = pvalues, colour = factor(Fission)),
          distribution = qunif, size = 1.5) +
  scale_colour_manual(name = "", 
                      values = c("#334EAC", "#BAD6EB")) +
  ylab("Empirical Quantiles") + 
  xlim(c(0, 1)) + ylim(c(0, 1)) + 
  theme_classic()  +
  theme(legend.position = "bottom")

plt_negbin <- (plt_illu_negbin + ggtitle("A") + theme(plot.title = element_text(face = "bold", size = 20))) + 
  (plt_typeI_negbin + ggtitle("B") + theme(plot.title = element_text(face = "bold", size = 20)))


#------------------------------------------------------------------------------#
#            Panel C: Simulations results under Correlated setting             #


neg_bin_cor_res <- read.csv("results/CorrelationAndRelativeBiais.csv")

plt_neg_bin_cor <- neg_bin_cor_res %>% 
  mutate(Theta_hat = Theta*Error) %>% 
  mutate(RelativeBiais = (Theta_hat-Theta)/Theta) %>% 
  group_by(Method, Estimation, Rho, RelativeBiais, Error) %>% 
  filter(Rho != 0.01) %>%
  filter(!(RelativeBiais<0 & Method == "NB")) %>%
  filter(RelativeBiais < 5) %>%
  summarise(FDR = mean(p.adjust(pval, method = "BH") < 0.05),
            typeI = mean(pval < 0.05)) %>%
  mutate(Name = paste(Method, Estimation, collapse = "_")) %>%
  mutate(
    LabelName = case_when(
      Name == "Gauss Oracle" ~ "Gaussian~(theta)",
      Name == "Gauss Wrong"  ~ "Gaussian~(widehat(theta))",
      Name == "NB Oracle"    ~ "Negative~Binomial~(theta)",
      Name == "NB Wrong"     ~ "Negative~Binomial~(widehat(theta))",
      TRUE ~ Name
    )
  ) %>%
  ggplot() +
  aes(x=RelativeBiais,
      y = typeI, 
      colour = as.factor(Rho)) +
  geom_point(size = 2) +
  geom_line(linewidth = .9) +
  facet_grid(~LabelName, labeller = label_parsed, scales = "free_x") +
  geom_hline(yintercept = 0.05,
             color = "grey20",
             linewidth = 1.2,
             linetype = "dashed") +
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
  xlab(TeX(r'($\frac{\widehat{\theta} - \theta}{\theta}$)')) +
  ylab("Type I error") + 
  theme_classic() +
  # scale_x_continuous(breaks = c(round(seq(-0.99, 9, length.out = 5), 0) ,0)) +
  theme(legend.position = 'bottom') +
  guides(
    colour = guide_legend(nrow = 1))

figure3 <- ((plt_illu_negbin + ggtitle("A")) + 
  (plt_typeI_negbin + ggtitle("B"))) /
  (plt_neg_bin_cor + ggtitle("C")) +
  plot_layout(heights = c(1.2, 1)) &
  theme(plot.title = element_text(face = "bold", size = 12),
        text = element_text(size = 13))


ggsave(figure3,
       filename = "Figures/figure3.pdf", 
       width = 185,
       height = 130,
       units = "mm")


#------------------------------------------------------------------------------#
#                                 Figure 4                                     #
#                         Applications on Bonne marrow                         #
#------------------------------------------------------------------------------#


Xtrain <- read_xlsx("results/ResultsOnScRNASeq.xlsx", sheet = "Xtrain")
Xtest <- read_xlsx("results/ResultsOnScRNASeq.xlsx", sheet = "Xtest")
cluster <- read_xlsx("results/ResultsOnScRNASeq.xlsx", sheet = "Cluster")
results <- read_xlsx("results/ResultsOnScRNASeq.xlsx", sheet = "Results")


#------------------------------------------------------------------------------#
#                 Panel A Gene-specific overdispersion parameters              #

pval_toQQ <- results %>% dplyr::select(Gene,
                                PvaluesAfterUnivariateClustering,
                                PValuesAfterMultivariateClustering) %>%
  tidyr::pivot_longer(
    cols = c(PvaluesAfterUnivariateClustering, PValuesAfterMultivariateClustering),
    names_to = "Method",
    values_to = "Pvalue"
  ) %>%
  mutate(Method = case_when(
    Method == "PvaluesAfterUnivariateClustering" ~ "Univariate \n Clustering",
    Method == "PValuesAfterMultivariateClustering" ~ "Multivariate \n Clustering"
  ))

plt_qqApplication <- ggplot(pval_toQQ) +
  geom_abline(slope=1, intercept=0, col="red", size = 1.2, alpha = .7) + 
  stat_qq(aes(sample = Pvalue, colour = factor(Method)),
          distribution = qunif, size = 1.5) +
  scale_colour_manual(name = "Method",
                      values = c("#8C33FF",  # bleu
                                 "#2CA02C")) +
  xlab("Theoretical Quantiles") + 
  ylab("Empirical Quantiles") + 
  xlim(c(0, 1)) + ylim(c(0, 1)) + 
  theme_classic() +
  theme(text = element_text(size = 14))

#------------------------------------------------------------------------------#
#         Panel B-C Correlation between genes lead to problems even on         #
#                         homogeneous sub-population                           #

# Correlation plots 
plt_cor_orig_thin <- ggplot(results) +
  geom_abline(slope=1, intercept=0, col="red", size = 1.2, alpha = .7) + 
  
  aes(x=OriginalCorWithGene1, 
      y = ThinningCorWithGene1TrainAndTest) +
  geom_point(size = 1,
             alpha = .5) +
  ylab(TeX(r'(Cor$\left(X^{(1)}_1, X^{(2)}_j\right)$)')) +
  xlab(TeX(r'(Cor$\left(X_1, X_j\right)$)')) +
  xlim(c(-1,1)) +
  ylim(c(-1,1)) +
  theme_classic()


## Sélection du gène
GeneToSel <- sort(
  results$PValuesAfterMultivariateClustering,
  index.return = TRUE
)$ix[1]

## Clustering univarié (Train)
km_uni <- kmeans(
  log2(Xtrain[, GeneToSel] + 1),
  centers = 2,
  nstart  = 100
)

## Réplication Train + Test
cluster_uni <- c(km_uni$cluster, km_uni$cluster)

Xtrain1 <- tibble(
  Gene1 = c(
    Xtrain[, GeneToSel, drop = TRUE],
    Xtest[, GeneToSel, drop = TRUE],
    Xtrain[, GeneToSel, drop = TRUE],
    Xtest[, GeneToSel, drop = TRUE]
  ),
  Cluster = factor(c(
    rep(cluster$Cluster, 2), # Multivarié
    cluster_uni # Univarié
  )),
  Method = factor(
    rep(c("Multivariate", "Univariate"),
        each = 2 * nrow(Xtrain))
  ),
  Thinning = factor(
    rep(rep(c("Train", "Test"), each = nrow(Xtrain)), 2),
    levels = c("Train", "Test"))
) %>%
  mutate(
    MethodName = paste(Method, "Clustering", sep = " ")
  )


pvals <- Xtrain1 %>%
  group_by(MethodName, Thinning) %>%
  summarise(
    p_value = wilcox.test(log2(Gene1 + 1) ~ Cluster)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(
      "Wilcoxon~italic(p)~\"=\"~\"",
      scales::pvalue(p_value),
      "\""
    ),
    y_pos = max(log2(Xtrain1$Gene1 + 1), na.rm = TRUE) * 1.25
  )

plt_illuApplication <- 
  Xtrain1%>%
  ggplot(
    aes(
      x    = Thinning,
      y    = log2(Gene1 + 1),
      fill = Cluster
    )
  ) +
  introdataviz::geom_split_violin(alpha = .5, trim = FALSE) +
  geom_boxplot(width = .2, alpha = 1, show.legend = FALSE) +
  
  ## P-values (solution robuste)
  geom_text(
    data = pvals,
    aes(
      x     = Thinning,
      y     = y_pos,
      label = label
    ),
    inherit.aes = FALSE,
    parse = TRUE,
    size  = 4
  ) +
  
  scale_fill_manual(
    name   = "PGK1",
    values = rev(colorspace::lighten(c("#294122", "#EB3D00"), 0.25)),
    labels = c(TeX(r'($C_1$)'), TeX(r'($C_2$)'))
  ) +
  scale_colour_manual(
    name   = "PGK1",
    values = rev(colorspace::lighten(c("#294122", "#EB3D00"), 0.25)),
    labels = c(TeX(r'($C_1$)'), TeX(r'($C_2$)'))
  ) +
  facet_grid(~MethodName) +
  ylab("log2(counts + 1)") +
  xlab("") +
  scale_x_discrete(
    labels = c(
      "Train" = expression(X^{(1)}),
      "Test"  = expression(X^{(2)})
    )
  )

plt_illuApplication

figure4 <- ((plt_cor_orig_thin + plt_qqApplication ) / (plt_illuApplication)) +
  plot_annotation(tag_levels = "A") +
  plot_layout(height = c(2, 4)) &
  theme_classic() &
  theme(text = element_text(size = 16),
        plot.tag = element_text(face = "bold", size = 16))

ggsave("Figures/figure4.pdf",
       plot = figure4, 
       width = 250, 
       height = 200, 
       units = "mm")


# Table pvalues 
typeI_table <- data.frame(Method = c("Multivariate Clustering",
                                     "Univariate Clustering"),
                          `Type I` = c(mean(results$PValuesAfterMultivariateClustering < 0.05),
                                       mean(results$PvaluesAfterUnivariateClustering < 0.05))) 
xt <- xtable::xtable(typeI_table, 
             caption = "Type I error rates from Wilcoxon tests after Negative Binomial data thinning. 
  In the multivariate setting, k-means clustering is performed on all 500 genes in X^{(1)}, 
  and Wilcoxon tests are applied gene-wise on X^{(2)}. 
  In the univariate setting, clustering is performed individually for each gene in X^{(1)}."
)

xtable::print.xtable(xt, 
      type = "latex",           # génère du LaTeX
      file = "results/TypeIErrorRateNBOnScRNASeq.txt",    # nom du fichier
      include.rownames = FALSE) 


