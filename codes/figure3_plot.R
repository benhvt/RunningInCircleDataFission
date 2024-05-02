# ------------------------------ Figure 3 FINALE ----------------------------- #

library(DataFission)
library(datathin)
library(ggplot2)
library(pbapply)
library(dplyr)
library(latex2exp)
library(patchwork)
theme_set(theme_bw())

################################################################################
#                                  Ratio                                       #
################################################################################

liste_csv <- list.files("results/figure3AB/", pattern = c(".csv"), full.names = TRUE)

all_csv <- lapply(liste_csv, function(l){read.csv(l, header = T)})
df_csv <- dplyr::bind_rows(all_csv)
df_ratio <- df_csv %>% 
  mutate(Method2 = ifelse(Method == "Individual hi", "Comaniciu et al.", Method)) %>%
  mutate(sigma = sigma) %>%
  group_by(delta, Method2, sigma) %>%
  # filter(Method == "Individual hi") %>%
  summarise(TypeI = mean(pval < 0.05), 
            MeanVariance= mean(Variance_hat),
            MeanMedian = mean(Median_Variance),
            sdVariance_hat = mean(sdVariance_hat),
            MeanQ1 = mean(Q1_Variance),
            MeanQ3 = mean(Q3_Variance))  %>%
  mutate(sigma_lab = paste0("sigma^2==", sigma^2)) %>%
  mutate(Ratio = delta/(sigma))  %>%
  filter(Method2 == "Change Point")

plot_typeI <- ggplot(df_ratio) +
  aes(x=Ratio, y = TypeI, colour = sigma_lab) +
  geom_line(size = 1.4) +
  scale_colour_manual(name = TeX(r'($\sigma^2$)'),
                      values = c("#93B5C6", "#DBC2CF","#998bc0", "#BD4F6C"),
                      labels = c(0.01, 0.5, 1, 4)) +
  ggnewscale::new_scale_colour() +
  geom_hline(aes(yintercept = 0.05, 
                 colour = "5% nominal level"),
             linetype = 2,
             size = 1.2) +
  facet_grid(~sigma_lab, scales = "free_x", labeller = labeller(sigma_lab = label_parsed)) +
  scale_colour_manual(name = "",
                      values = "#6C0E23") +
  scale_x_log10() +
  xlab(TeX(r'(Ratio $\delta/\sigma$)')) +
  annotation_logticks(sides = "b") +
  ylab("Empirical \n Type I error rate") +
  ylim(c(0,1)) +
  theme_classic() +
  NULL


plot_est_Med <- df_ratio %>%
  # filter(Method == "Change Point") %>%
  ggplot() +
  aes(x=Ratio, y = (MeanMedian - sigma^2)/sigma^2, colour = sigma_lab) +
  scale_colour_manual(name = TeX(r'($\sigma^2$)'),
                      values = c("#93B5C6", "#DBC2CF", "#998bc0", "#BD4F6C"),
                      labels = c(0.01, 0.5, 1, 4)) +
  geom_hline(aes(yintercept = 0), colour = "#6C0E23", size = 1, linetype = "dotted") + 
  geom_vline(aes(xintercept = 5), colour = "black", size = 1, linetype = 2) +
  geom_ribbon(aes(ymin = (MeanQ1 - sigma^2)/sigma^2,
                  ymax = (MeanQ3 - sigma^2)/sigma^2, fill = sigma_lab, colour = sigma_lab), alpha = .2) +
  scale_fill_manual(name = TeX(r'($\sigma^2$)'),
                    values = c("#93B5C6", "#DBC2CF", "#998bc0", "#BD4F6C"),
                    labels = c(0.01, 0.5, 1, 4)) +
  geom_line(size = 1.4) +
  scale_x_log10() +
  annotation_logticks(sides = "b") +
  xlab(TeX(r'(Ratio $\delta/\sigma$)')) +
  ylab("Relative bias") +
  facet_grid(~sigma_lab, scales = "free_x", labeller = labeller(sigma_lab = label_parsed)) +
  theme_classic() +
  NULL


################################################################################
#                               Sample Size                                    #
################################################################################


liste_csv_samp <- list.files("results/figure3CD/", pattern = c(".csv"), full.names = TRUE)

all_csv_samp <- lapply(liste_csv_samp, function(l){read.csv(l, header = T)})
df_csv_samp <- dplyr::bind_rows(all_csv_samp)

df_samp <- df_csv_samp %>% 
  mutate(Method2 = ifelse(Method == "Individual hi", "Comaniciu et al.", Method)) %>%
  mutate(sigma = sigma) %>%
  group_by(n, delta, Method2, sigma) %>%
  summarise(TypeI = mean(pval < 0.05), 
            MeanVariance= mean(Variance_hat),
            MeanMedian = mean(Median_Variance),
            sdVariance_hat = mean(sdVariance_hat),
            MeanQ1 = mean(Q1_Variance),
            MeanQ3 = mean(Q3_Variance))  %>%
  mutate(sigma_lab = paste0("sigma^2==", sigma^2)) %>%
  mutate(Ratio = delta/(sigma)) %>%
  filter(Method2 == "Change Point") %>%
  mutate(delta2 = paste('frac(delta, sigma)', delta, sep = "=="))

plot_typeI_samp <- df_samp %>%
  ggplot() + 
  aes(x=n, y = TypeI, colour = as.factor(delta)) +
  geom_line(size = 1.4) +
  scale_colour_manual(name = TeX(r'($\delta/\sigma$)'),
                      values = c("#B0DAF1", "#6E75A8", "#3E1929"), 
                      labels = c(0.5, 5, 10)) +
  ggnewscale::new_scale_colour() +
  geom_hline(aes(yintercept = 0.05, 
                 colour = "5% nominal level"),
             linetype = 2,
             size = 1.2) +
  scale_colour_manual(name = "",
                      values = "#6C0E23") +
  facet_grid(~factor(delta2, levels = paste("frac(delta, sigma)", c(0.5, 5, 10), sep = "==")), labeller = label_parsed) +
  xlab("Sample Size") +
  ylab("Empirical \n Type I error rate") +
  ylim(c(0,1)) +
  theme_classic() +
  NULL
plot_typeI_samp

plot_est_Med_samp <- df_samp %>%
  ggplot() +
  aes(x=n, y = (MeanMedian - sigma^2)/sigma^2, colour = as.factor(delta)) +
  scale_colour_manual(name = TeX(r'($\delta/\sigma$)'),
                      values = c("#B0DAF1", "#6E75A8", "#3E1929"), 
                      labels = c(0.5, 5, 10)) +
  geom_hline(aes(yintercept = 0), colour = "#6C0E23", size = 1, linetype = "dotted") + 
  geom_ribbon(aes(ymin = (MeanQ1- sigma^2)/sigma^2,
                  ymax = (MeanQ3 - sigma^2)/sigma^2, fill = as.factor(delta), colour = as.factor(delta)), alpha = .2) +
  scale_fill_manual(name = TeX(r'($\delta/\sigma$)'),
                    values = c("#B0DAF1", "#6E75A8", "#3E1929"), 
                    labels = c(0.5, 5, 10)) +
  geom_line(size = 1.4) +
  facet_grid(~factor(delta2, levels = paste("frac(delta, sigma)", c(0.5, 5, 10), sep = "==")), labeller = label_parsed) +
  xlab("Sample Size") +
  ylab("Relative bias") +
  theme_classic() +
  NULL
(plot_est_Med + plot_typeI)/(plot_est_Med_samp + plot_typeI_samp) + 
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(face = "bold", size = 24),
        text = element_text(size = 22),
        axis.text.x = element_text(size = 18, angle = 45, hjust =1),
        legend.position = "bottom")


ggsave(filename = "figures/figure3.pdf", 
       width = 400, 
       height = 250, 
       units = "mm", 
       dpi = 600) 
