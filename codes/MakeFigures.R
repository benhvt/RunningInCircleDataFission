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

theme_set(theme_classic())

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

plt_power <- ggplot(powerResults) +
  aes(x=tau, y = Power, 
      colour = as.factor(n),
      linetype = Fission_lab) +
  geom_line(linewidth = 1.1, 
            alpha = .8) +
  scale_colour_manual(name = "Sample Size",
                      values = MetBrewer::met.brewer("Derain", n=7)) +
  scale_linetype_manual(name = "Data fission",
                        values = c(1, 6)) +
  ggnewscale::new_scale_colour() +
  facet_grid(~Variable_lab, 
             labeller = label_parsed) +
  xlab(TeX(r'( Tunning parameter $\tau$)')) +
  ylab("Statistical Power \n 5% level")

plt_ari <- ggplot(ariResults) +
  aes(x=tau, y = ARI_m, 
      colour = as.factor(n),
      linetype = Fission_lab) +
  geom_line(linewidth = 1.1, 
            alpha = .8) +
  scale_colour_manual(name = "Sample Size",
                      values = MetBrewer::met.brewer("Derain", n=7)) +
  scale_linetype_manual(name = "Data fission",
                        values = c(1,6)) +
  ggnewscale::new_scale_colour() +
  xlab(TeX(r'( Tunning parameter $\tau$)')) +
  ylab("Adjusted Rand Index")



#------------------------------------------------------------------------------#
#                      Second panel of the Figure                              #
#                 Type I error rate in the worst scenario                      #


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
  xlab(TeX(r'($X_1$)')) +
  ylab(TeX(r'($X_2$)'))

# Import files of results
adverse_sc <- read.csv(file = "results/AdverseScenario_GaussianFission.csv")

plt_typeI <- adverse_sc %>%
  mutate(Fission_lab = paste(Fission, "fission", sep = " ")) %>%
  group_by(Fission_lab,Variable, tau, n) %>%
  summarise(TypeI = mean(pvalues < 0.05)) %>%
  mutate(Variable_lab = ifelse(Variable == "X1", "X[1]", "X[2]")) %>%
  ggplot() +
  aes(x=tau, y = TypeI, 
      colour = as.factor(n),
      linetype = Fission_lab) +
  geom_line(linewidth = 1.1, 
            alpha = .8) +
  scale_colour_manual(name = "Sample Size",
                      values = MetBrewer::met.brewer("Derain", n=7)) +
  scale_linetype_manual(name = "Data fission",
                        values = c(1, 6)) +
  ggnewscale::new_scale_colour() +
  geom_hline(aes(yintercept = 0.05, 
                 colour = "5% nominal level"),
             linewidth = .9, 
             linetype = "dashed") +
  scale_colour_manual(name = "",
                      values = "#DB2763") +
  facet_grid(~Variable_lab, 
             labeller = label_parsed) +
  xlab(TeX(r'( Tunning parameter $\tau$)')) +
  ylab("Type I error rate \n 5% nominal levels") 

figure2 <-  (plt_illu_power + ((plt_power / plt_ari) + plot_layout(heights = c(2, 1.5)))) /
  (plt_illu_typeI + plt_typeI) +
  plot_layout(guides = "collect", heights = c(2,1.5)) +
  plot_annotation(tag_levels = 'A') &
  theme_classic(base_size = 14) +
  theme(plot.tag = element_text(face = "bold"))
figure2
ggsave(plot = figure2,
       filename = "Figures/figure2.pdf",
       width = 300, 
       height = 200, 
       units = "mm")

#------------------------------------------------------------------------------#
#                                  Figures                                     #
#               Impact of variance estimation on Type I error rate             # 
#------------------------------------------------------------------------------#


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

ggsave(filename = "Figures/figure3.pdf",
       width = 350, 
       height = 120, 
       units = "mm",
       dpi = 600)

#------------------------------------------------------------------------------#
#                                  Figures                                     #
#         Data fission for scRNA-seq (Binomial Negatives Mixture)              # 
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
#                     Panel A and B: Simulations results                       #

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
                                   colour = as.factor(cl_X), 
                                   shape = as.factor(Z)) + 
  geom_point(size = 3) +
  scale_shape_manual(name = "True classes",
                     values = c(15,16)) +
  scale_colour_manual(name = "Clusters", 
                      values = c("#294122", "#EB3D00", "#FFBBA6"),
                      labels = c(TeX(r'($C_1$)'),
                                 TeX(r'($C_2$)'),
                                 TeX(r'($C_3)'))) +
  xlab(TeX(r'($X_1$)')) +
  ylab(TeX(r'($X_2$)')) +
  theme_classic() +
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
#                         Panel C: Applications results                        #

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

plt_application <- ((allPairplot[[1]]+ ggtitle("C") + theme(plot.title = element_text(face = "bold", size = 20))) + plot_spacer()  + allPairplot[[2]] + plot_spacer()  + allPairplot[[3]] + plot_spacer()  + allPairplot[[4]]) + plot_layout(nrow = 1, widths = c(4,.5,4,.5,4,.5,4)) 


(plt_negbin) / (plt_application) &  #+ plot_annotation(tag_levels = "A") + plot_layout(tag_level = "new") & 
  theme(text = element_text(size = 20))


ggsave("Figures/figure4.png", width = 315, height = 200, units = "mm")
