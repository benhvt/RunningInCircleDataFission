# --------------------------------- Figure 4 --------------------------------- #
# --------------------- Negative Binomial Simulations -------------------------# 
library(DataFission)
library(datathin)
library(ggplot2)
library(pbapply)
library(dplyr)
library(latex2exp)
library(patchwork)
theme_set(theme_bw())

# ------------------------------ Simulations study --------------------------- #

#-- Functions 
source("utils.R")

#-- Parameters 
n <- 100
true_class <- rep(1:2, each = n/2)
nsimu <- 1000

#-- Simulations with the wrong number of classes 
sim_res_wrong <- pblapply(1:nsimu, function(ns){
  set.seed(310123*ns)
  x1 <- rnbinom(n=n/2, prob = 0.5, size = 5)
  x2 <- rnbinom(n=n/2, prob = 0.4, size = 40)
  y <- rnbinom(n=n, prob = .5, size = 5)
  X <- cbind(c(x1,x2), 
             y)
  
  #--- Data thinning with the global parameter
  overdisp_global <- sapply(1:ncol(X), function(p){
    npreg::theta.mle(X[,p], mu = mean(X[,p]))
  })
  
  res <- countsplit::countsplit(X, overdisps = overdisp_global)
  clust <- km_fun(as.matrix(res[[1]]), K=3)
  clToTest <- order_cluster(as.matrix(res[[1]])[,1], clust)
  res_test <- wilcox.test(as.matrix(res[[2]])[clust == clToTest[1],1], 
                          as.matrix(res[[2]])[clust == clToTest[2],1])$p.value
  
  #--- Data thinning with the true intra-class parameter 
  overdisp_C1 <- sapply(1:ncol(X), function(p){
    npreg::theta.mle(X[true_class==1,p], mu = mean(X[true_class==1,p]))
  })
  
  overdisp_C2 <- sapply(1:ncol(X), function(p){
    npreg::theta.mle(X[true_class==2,p], mu = mean(X[true_class==2,p]))
  })
  overdisp_intra <- list(overdisp_C1, overdisp_C2)
  res_temp <- lapply(1:2, function(cl){
    res_datathin <- countsplit::countsplit(X[true_class==cl, ], overdisps = overdisp_intra[[cl]])
  })
  
  Xtrain <- rbind(as.matrix(res_temp[[1]][[1]]),
                  as.matrix(res_temp[[2]][[1]]))
  Xtest <- rbind(as.matrix(res_temp[[1]][[2]]),
                 as.matrix(res_temp[[2]][[2]]))
  
  clust_intra <- km_fun(Xtrain, K=3)
  clToTest_intra <- order_cluster(Xtrain[,1], clust_intra)
  res_test_intra <- wilcox.test(Xtest[clust_intra==clToTest_intra[1],1],
                          Xtest[clust_intra==clToTest_intra[2],1])$p.value  
  
  #--- Data thinning with the wrong intra-cluster parameter 
  cl_ref <- km_fun(X, K=3)
  
  overdisp_C1_hat <- sapply(1:ncol(X), function(p){
    npreg::theta.mle(X[cl_ref==1,p], mu = mean(X[cl_ref==1,p]))
  })
  
  overdisp_C2_hat <- sapply(1:ncol(X), function(p){
    npreg::theta.mle(X[cl_ref==2,p], mu = mean(X[cl_ref==2,p]))
  })
  
  overdisp_C3_hat <- sapply(1:ncol(X), function(p){
    npreg::theta.mle(X[cl_ref==3,p], mu = mean(X[cl_ref==3,p]))
  })
  
  overdisp_intra_hat <- list(overdisp_C1_hat, overdisp_C2_hat, overdisp_C3_hat)
  res_temp_hat <- lapply(1:3, function(cl){
    res_datathin <- countsplit::countsplit(X[cl_ref==cl, ], overdisps = overdisp_intra_hat[[cl]])
  })
  
  Xtrain_hat <- rbind(as.matrix(res_temp_hat[[1]][[1]]),
                  as.matrix(res_temp_hat[[2]][[1]]),
                  as.matrix(res_temp_hat[[3]][[1]]))
  Xtest_hat <- rbind(as.matrix(res_temp_hat[[1]][[2]]),
                 as.matrix(res_temp_hat[[2]][[2]]),
                 as.matrix(res_temp_hat[[3]][[2]]))
  
  
  clust_intra_hat <- km_fun(Xtrain_hat, K=3)
  clToTest_intra_hat <- order_cluster(Xtrain_hat[,1], clust_intra_hat)
  res_test_intra_hat <- wilcox.test(Xtest_hat[clust_intra_hat==clToTest_intra_hat[1],1],
                                Xtest_hat[clust_intra_hat==clToTest_intra_hat[2],1])$p.value  
  
  return(data.frame(pval = c(res_test, res_test_intra, res_test_intra_hat),
                    Overdispersion = c("Global", "Intra-comp", "Intra-cluster")))
}, cl = 6)


#-- Make figure
cluster_col <- c("#294122", "#EB3D00", "#FFBBA6")
results_col <- c("#334EAC", "#BAD6EB", "#6A2A23")
set.seed(09022024)

x1 <- rnbinom(n=n/2, prob = 0.5, size = 5)
x2 <- rnbinom(n=n/2, prob = 0.4, size = 40)
y <- rnbinom(n=n, prob = .5, size = 5)
X <- cbind.data.frame(X1=c(x1,x2), 
           X2=y, 
           TrueClasses=as.factor(true_class))
cl_X <- km_fun(cbind(c(x1,x2), y), K = 3)

plt1 <- ggplot(X) + aes(x=X1, y=X2, colour = cl_X, shape = as.factor(true_class)) + 
  geom_point(size = 3) +
  scale_shape_manual(name = "True classe",
                    values = c(15,16)) +
  scale_colour_manual(name = "Clusters", 
                      values = cluster_col,
                      labels = c(TeX(r'($C_1$)'),
                                 TeX(r'($C_2$)'),
                                 TeX(r'($C_3)'))) +
  xlab(TeX(r'($X_1$)')) +
  ylab(TeX(r'($X_2$)')) +
  theme_classic() +
  # theme(legend.position = "bottom") +
  NULL


plt2 <- do.call(rbind.data.frame, sim_res_wrong) %>% 
  ggplot() + 
  geom_abline(slope=1, intercept=0, col="red", size = 1.2, alpha = .7) + xlab("Theoretical Quantiles") + 
  stat_qq(aes(sample = pval, colour = factor(Overdispersion,
                                             levels = c('Intra-comp', 'Intra-cluster', 'Global'))),
              distribution = qunif, size = 1.5) +
  scale_colour_manual(name = "Overdispersion", 
                      values = results_col,
                      labels=c(TeX(r'($\hat{\theta}_{k}$)'),
                               TeX(r'($\hat{\theta}_{\hat{k}}$)'),
                               TeX(r'($\hat{\theta}$)'))) +
  ylab("Empirical Quantiles") + 
  xlim(c(0, 1)) + ylim(c(0, 1)) + 
  theme_classic()  +
  theme(legend.position = "bottom")
plt2

plot_illustration <- (plt1 + ggtitle("A") + theme(plot.title = element_text(face = "bold", size = 20))) + 
  (plt2 + ggtitle("B") + theme(plot.title = element_text(face = "bold", size = 20)))
  # theme(
  # #axis.title = element_text(size = 24),
  # #       axis.text = element_text(size = 18),
  # #       legend.text = element_text(size = 20),
  # #       legend.title = element_text(size = 24),
  #       plot.tag = element_text(face = "bold"))
# ------------------------------ Applications --------------------------- #
# 
# overdispersion_estimate <- read.csv(file ="results/overdispersion_estimate.csv")
# 
# R2 <- overdispersion_estimate  %>%
#   filter(Global < 4, Intra < 4) %>%
#   group_by(CellType) %>%
#   summarise(R2 = cor(Intra, Global))
# 
# plt_application <- ggplot(overdispersion_estimate) + 
#   aes(x=Intra, y = Global, colour = CellType) +
#   geom_point(size = 2, alpha = .5) + 
#   xlim(c(0, 4)) +
#   # scale_x_log10() +
#   ylim(c(0, 4)) +
#   # scale_y_log10() +
#   xlab(TeX(r'(Intra-cell population overdispersion $\hat{\theta}_k$)')) +
#   ylab(TeX(r'(Global overdispersion $\hat{\theta}$)')) +
#   scale_colour_manual(name = "Cell type", 
#                       values = c('#5C0029', 
#                                  "#E63946",
#                                  "#A8DADC",
#                                  "#457B9D")) +
#   facet_grid(~CellType) +
#   geom_text(data = R2, aes(x=1.5, y = 3.5, label = paste0("R2=", round(R2,3))), colour = "black", size = 8) +
#   geom_abline(slope = 1, intercept = 0, colour = "darkred", size = 1.2, linetype = "dashed") +
#   theme_classic() +
#   theme(legend.position = "bottom") +
#   NULL  
# 
# plt_final <- (plt1 + plt2) / plt_application +
#   plot_annotation(tag_levels = "A") &
#   theme(axis.title = element_text(size = 24), 
#         axis.text = element_text(size = 18),
#         legend.text = element_text(size = 20),
#         legend.title = element_text(size = 24),
#         plot.tag = element_text(face = "bold", size = 24))
# plt_final
# 
# ggsave(plt_final, filename = "figures/figure4.pdf",
#        width = 350, 
#        height = 240, 
#        units = "mm",
#         dpi = 600)
