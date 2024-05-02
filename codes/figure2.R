# --------------------------------- Figure 2 --------------------------------- #

library(DataFission)
library(datathin)
library(ggplot2)
library(pbapply)
library(dplyr)
library(latex2exp)
library(patchwork)
theme_set(theme_bw())

#-- Functions 
source("utils.R")
sim_fun <- function(n,sigma, sigma_hat, tau, nsim){
  res <- pbsapply(1:nsim, function(ns){
    X <- rnorm(n, mean = 0, sd = sigma)
    Z <- rnorm(n, mean = 0, sd = sigma_hat)
    X1 <- X + tau*Z
    X2 <- X - (1/tau)*Z
    
    cl <- km_fun(X1, K=2)
    test <- t.test(X2~cl)
    return(test$p.value)
  }, cl = parallel::detectCores()-1)
  return(mean(res < 0.05))
}

sim_fun_sigma_grid <- function(n,sigma, sigma_hat, tau, nsim){
  res <- pbsapply(sigma_hat, function(s){sim_fun(n, sigma, sigma_hat = s, tau, nsim)})
  return(res)
}

sim_fun_sigma <- function(n,sigma, sigma_hat, tau, nsim){
  res <- pblapply(sigma, function(s){sim_fun_sigma_grid(n, sigma = s, sigma_hat, tau, nsim)})
  return(res)
}

#-- Simulations  impact of sigma
n <- 100
sigma <- c(0.1, 0.5, 1, 2)
sigma_grid <- sort(c(seq(0, 4, length.out = 50), 2, 0.1))
nsimu <- 1000
tau <- .4

sim_res <- sim_fun_sigma(n=n, sigma = sigma, sigma_hat = sigma_grid, tau = tau, nsim = nsimu) 

df_res_sigma <- data.frame(EmpTypeI = unlist(sim_res),
                           sigma = rep(sigma, each = length(sigma_grid)),
                           sigma_hat = rep(sigma_grid, length(sigma))) %>%
  mutate(TheTypeI = compute_typeI(n = n, tau = tau, sigma = sigma, sigma_hat = sigma_hat, alpha = .05)) %>%
  mutate(sigma_name = paste0("sigma^2==", sigma^2)) %>%
  mutate(Ratio = (sigma^2-sigma_hat^2)/sigma^2) %>%
  reshape2::melt(id.vars = c("Ratio", "sigma_name", "sigma", "sigma_hat"), measure_vars = c("TheTypeI", "EmpTypeI"))

plot_sigma <-  ggplot(df_res_sigma) + aes(x=Ratio, y = value, colour = sigma_name) +
  geom_point(data = subset(df_res_sigma, variable == "EmpTypeI"), aes(shape = "Empirical"), size = 4) +
  scale_shape_manual(name = "Empirical", values = 2, labels = '', guide = "legend") +
  geom_line(data = subset(df_res_sigma, variable == "TheTypeI"), aes(group = sigma_name, linetype = "Theoritical"), size = 1) +
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



sim_fun_n <- function(n_grid,sigma, sigma_hat, tau, nsim){
  res <- pblapply(n_grid, function(n){sim_fun_sigma_grid(n=n, sigma = sigma, sigma_hat, tau, nsim)})
  return(res)
}

#-- Simulations: Impact of n

n <- c(50, 100, 200, 500, 1000)
sigma <- 1
sigma_grid <- seq(0.8, 1.8, length = 50)
nsimu <- 1000
tau <- .4

sim_res_n <- sim_fun_n(n_grid=n, sigma = sigma, sigma_hat = sigma_grid, tau = tau, nsim = nsimu) 

df_res_n <-  data.frame(EmpTypeI = unlist(sim_res_n),
                        sigma = sigma,
                        SampSize = rep(n, each = length(sigma_grid)),
                        sigma_hat = rep(sigma_grid, length(n))) %>%
  mutate(TheTypeI = compute_typeI(n = SampSize, tau = tau, sigma = sigma, sigma_hat = sigma_hat, alpha = .05)) %>%
  mutate(sigma_name = paste0("sigma^2==", sigma^2)) %>%
  mutate(Ratio = (sigma^2-sigma_hat^2)/sigma^2) %>%
  reshape2::melt(id.vars = c("Ratio", "sigma_name", "SampSize", "sigma", "sigma_hat"), measure_vars = c("TheTypeI", "EmpTypeI"))

#-- Plot 
plot_n <- ggplot(df_res_n) + aes(x=Ratio, y = value, colour = as.factor(SampSize)) +
  geom_point(data = subset(df_res_n, variable == "EmpTypeI"), aes(shape = "Empirical"), size = 4) +
  scale_shape_manual(name = "Empirical", values = 2, labels = '', guide = "legend") +
  geom_line(data = subset(df_res_n, variable == "TheTypeI"), aes(group = SampSize, linetype = "Theoritical"), size = 1) +
  scale_linetype_manual(name = "Theoritical", values = 1, labels = "", guide = "legend") +
  scale_colour_manual(name = "", 
                      values = colorRampPalette(c("#008154", "#0092a4", "#2a2956"))(length(n)),
                      labels = paste0("n=", n)) +
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

(plot_sigma + plot_n)  +
  plot_annotation(tag_levels = "A") &
  theme(axis.title = element_text(size = 24), 
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        plot.tag = element_text(face = "bold", size = 24),
        legend.spacing = unit(0.01, "cm"))
ggsave(filename = "figures/figure2.pdf",
       width = 350, 
       height = 120, 
       units = "mm",
       dpi = 600)
             