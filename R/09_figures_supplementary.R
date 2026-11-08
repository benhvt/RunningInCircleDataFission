# Supplementary Figures S1-S9 and table of Appendix A7 (from the intermediate
# results in results/). Figure S10 is a drawing (figures/FigureS10.png).

source("R/utils.R")
library(dplyr)
library(ggplot2)
library(patchwork)
library(latex2exp)
theme_set(theme_classic())
dir.create("figures", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)
read_res <- function(file) readRDS(file.path(res_dir, file))
save_fig <- function(plot, name, width, height) ggsave(paste0("figures/", name, ".pdf"), plot, width = width, height = height, units = "mm")

#------------------------------------------------------------------------------#
# Figure S1: Type I error rate after data fission with a biased variance

typeI_theoretical <- function(n, sigma, sigma_hat, tau = 0.4) {
  # asymptotic Type I error rate at the 5% level (Appendix A5, equation (7))
  rho <- (sigma^2 - sigma_hat^2) / sqrt((sigma^2 + tau^2 * sigma_hat^2) * (sigma^2 + sigma_hat^2 / tau^2))
  m <- sqrt(n) * sqrt(2 / pi * rho^2) / sqrt(1 - 2 / pi * rho^2)
  pnorm(qnorm(0.975), m, lower.tail = FALSE) + pnorm(-qnorm(0.975), m)
}

plot_S1 <- function(file, theo, group, colours, labels) {
  # empirical (points) and theoretical (lines, on the `theo` grid) Type I error
  emp <- read_res(file) %>%
    group_by(n, sigma, sigma_hat) %>%
    summarise(typeI = mean(pvalues < 0.05, na.rm = TRUE), .groups = "drop")
  ggplot(emp, aes((sigma^2 - sigma_hat^2) / sigma^2, typeI, colour = factor(.data[[group]]))) +
    geom_point(aes(shape = "Empirical"), size = 4) +
    geom_line(data = if (is.null(theo)) emp else theo,
              aes(y = typeI_theoretical(n, sigma, sigma_hat), linetype = "Theoretical"), linewidth = 1) +
    scale_shape_manual(name = "Empirical", values = 2, labels = "") +
    scale_linetype_manual(name = "Theoretical", values = 1, labels = "") +
    scale_colour_manual(name = "", values = colours, labels = labels) +
    ggnewscale::new_scale_colour() +
    geom_hline(aes(yintercept = 0.05, colour = "5% nominal levels"), linetype = 2, linewidth = 1.2) +
    scale_colour_manual(name = "", values = "#6C0E23") +
    labs(x = TeX("$(\\sigma^2 - \\widehat{\\sigma^2})/\\sigma^2$"), y = "Type I error rate") +
    guides(shape = guide_legend(order = 1), linetype = guide_legend(order = 2), colour = guide_legend(order = 3))
}

S1A <- plot_S1("sim02_variance_bias_S1A.rds", group = "sigma",
               theo = expand.grid(sigma_hat = seq(0, 12, length.out = 1000), sigma = c(0.1, 0.5, 1, 2), n = 100),
               colours = c("#93B5C6", "#DBC2CF", "#998bc0", "#BD4F6C"),
               labels = parse(text = paste("sigma^2 ==", c(0.01, 0.25, 1, 4)))) + xlim(-5, 2)
S1B <- plot_S1("sim02_variance_bias_S1B.rds", theo = NULL, group = "n",
               colours = colorRampPalette(c("#008154", "#0092a4", "#2a2956"))(5),
               labels = paste0("n=", c(50, 100, 200, 500, 1000)))
save_fig(S1A + S1B + plot_annotation(tag_levels = "A") &
           theme(axis.title = element_text(size = 24), axis.text = element_text(size = 18),
                 legend.text = element_text(size = 16), legend.title = element_text(size = 18),
                 plot.tag = element_text(face = "bold", size = 24), legend.spacing = unit(0.01, "cm")),
         "FigureS1", 350, 120)

#------------------------------------------------------------------------------#
# Figures S2-S3: QQ-plots of the Wilcoxon p-values of the first and last
# variables (Appendix A6), with oracle overdispersion (S2) or without
# correlation (S3)

pvals <- read_res("sim04_correlation_bias.rds")$pvalues %>%
  mutate(VarName = ifelse(Variable == "Var1", "X[1]", "X[50]"),
         MethodName = ifelse(Method == "NB", "Negative~Binomial~Thinning", "Gaussian~Fission"))
plot_qq <- function(df, colour) {
  ggplot(df) +
    geom_abline(slope = 1, intercept = 0, col = "red", linewidth = 1.2, alpha = .7) +
    stat_qq(aes(sample = pval, colour = factor(.data[[colour]])), distribution = qunif, size = 1, alpha = .5) +
    facet_grid(VarName ~ MethodName, labeller = label_parsed) +
    labs(x = "Theoretical Quantiles", y = "Empirical Quantiles") +
    lims(x = c(0, 1), y = c(0, 1)) +
    theme(text = element_text(size = 14))
}
save_fig(plot_qq(filter(pvals, Estimation == "Oracle", Error == 1), "Rho") +
           scale_colour_manual(name = TeX("Correlation ($\\rho$)"), values = rho_colours),
         "FigureS2", 200, 150)
save_fig(plot_qq(filter(pvals, Estimation == "Wrong", Rho == 0, Error %in% c(0.001, 0.5, 1, 1.5, 5)) %>%
                   mutate(RelativeBias = Error - 1), "RelativeBias") +
           MetBrewer::scale_color_met_d("Hiroshige") + labs(colour = "Relative Bias"),
         "FigureS3", 200, 150)

#------------------------------------------------------------------------------#
# Figures S4 and S9: GLM-PCA and UMAP of the Tabula Sapiens datasets

plot_glmpca_umap <- function(counts, true_cell, size, alpha) {
  set.seed(20260311)
  factors <- glmpca::glmpca(Y = t(counts), L = 50, fam = "nb")$factors
  config <- umap::umap.defaults
  config$n_neighbors <- 80
  config$min_dist <- 0.6
  layout <- umap::umap(factors, config = config)$layout
  panel <- function(x, y, lab) {
    ggplot(data.frame(x, y, TrueCell = true_cell), aes(x, y, colour = TrueCell)) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_vline(xintercept = 0, linetype = "dashed") +
      geom_point(size = size, alpha = alpha) +
      scale_colour_manual(name = "True cells", values = c(granulocyte = "#294122", `cd24 neutrophil` = "#EB3D00")) +
      labs(x = paste(lab, 1), y = paste(lab, 2))
  }
  panel(factors[, 1], factors[, 2], "GLM-PC") + panel(layout[, 1], layout[, 2], "UMAP") +
    plot_layout(guides = "collect") & theme(legend.position = "bottom", text = element_text(size = 14))
}
app6 <- read_res("app06_granulocytes.rds")
app7 <- read_res("app07_granulocytes_neutrophils.rds")
save_fig(plot_glmpca_umap(app7$counts, app7$true_cell, size = 1, alpha = .6), "FigureS4", 150, 75)
save_fig(plot_glmpca_umap(app6$counts, "granulocyte", size = 1.5, alpha = .8), "FigureS9", 150, 75)

#------------------------------------------------------------------------------#
# Figure S5 and table: conditional vs marginal NB data thinning on the mixture
# of granulocytes and CD24 neutrophils (Appendix A7)

res <- app7$results %>%
  mutate(ClusteringName = paste(Clustering, "Clustering"), ThinningName = paste(Thinning, "Thinning"))
plot_S5 <- function(hyp) {
  ggplot(filter(res, Hyp == hyp), aes(-log10(TruePvalue), -log10(pvalue), colour = ARI)) +
    geom_point(alpha = .5, size = 2) +
    geom_hline(yintercept = -log10(0.05)) +
    facet_grid(ClusteringName ~ ThinningName, scales = "free") +
    scale_colour_viridis_c(option = "plasma", direction = -1, limits = c(min(res$ARI), 1),
                           guide = guide_colourbar(barwidth = unit(5, "cm"))) +
    labs(x = "-log10(True p-value)", y = "-log10(p-value)",
         title = TeX(paste0("Distribution of p-values for genes under $H_", substr(hyp, 2, 2), "$")))
}
top5 <- res %>% filter(Hyp == "H1") %>% group_by(ThinningName, ClusteringName) %>% arrange(pvalue) %>% slice_head(n = 5)
save_fig(plot_S5("H0") / (plot_S5("H1") + ggrepel::geom_label_repel(data = top5, aes(label = GeneSymbol), seed = 1)) +
           plot_layout(guides = "collect") & theme_bw() + theme(legend.position = "bottom", text = element_text(size = 14)),
         "FigureS5", 150, 250)

performances <- res %>%
  group_by(Thinning, Clustering) %>%
  summarise(TypeI = mean(pvalue[Hyp == "H0"] < 0.05), Power = mean(pvalue[Hyp == "H1"] < 0.05), .groups = "drop")
write.csv(performances, "tables/Table_A7_performances.csv", row.names = FALSE)

#------------------------------------------------------------------------------#
# Figures S6-S7: marginal vs conditional Gaussian data fission, in the ideal
# (S6) and adverse (S7) scenarios

set.seed(20250404)  # illustrative example (panels A)
illu <- simulate_mixture(500)
cl_ideal <- kmeans(illu$X, centers = 2, nstart = 100)$cluster
cl_adverse <- rep(3, 500)
cl_adverse[illu$Z == 1] <- kmeans(illu$X[illu$Z == 1, ], centers = 2, nstart = 100)$cluster

plot_illustration <- function(cl, colours) {
  ggplot(data.frame(X1 = illu$X[, 1], X2 = illu$X[, 2], Cluster = factor(cl)), aes(X1, X2, colour = Cluster)) +
    geom_point(size = 2, alpha = .8) +
    scale_colour_manual(name = "Estimated clusters", values = colours) +
    labs(x = TeX("$X_1$"), y = TeX("$X_2$")) +
    theme(legend.position = "bottom")
}
size_colours <- scale_colour_manual(name = "Sample Size", values = MetBrewer::met.brewer("Derain", n = 7))
tau_lab <- labs(x = TeX("Tuning parameter $\\tau$"))
plot_rejection <- function(df, ylab) {
  df %>%
    group_by(Fission, Variable, tau, n) %>%
    summarise(rate = mean(pvalues < 0.05, na.rm = TRUE), .groups = "drop") %>%
    mutate(Fission = paste0(Fission, "~fission"), Variable = ifelse(Variable == "X1", "X[1]", "X[2]")) %>%
    ggplot(aes(tau, rate, colour = factor(n))) +
    geom_line(linewidth = .75, alpha = .8) +
    size_colours +
    guides(colour = guide_legend(nrow = 2)) +
    ggnewscale::new_scale_colour() +
    geom_hline(aes(yintercept = 0.05, colour = "5% nominal level"), linewidth = .75, linetype = "dashed", show.legend = FALSE) +
    scale_colour_manual(name = "", values = "#DB2763") +
    facet_grid(Variable ~ Fission, labeller = label_parsed) +
    tau_lab + labs(y = ylab) +
    theme(legend.position = "bottom")
}

ideal <- read_res("sim01_ideal.rds")
ari <- ideal %>%
  filter(Variable == "X1") %>%
  group_by(Fission = paste(Fission, "fission"), tau, n) %>%
  summarise(ARI = mean(ARI), .groups = "drop") %>%
  ggplot(aes(tau, ARI, colour = factor(n))) +
  geom_line(linewidth = .75, alpha = .8) +
  facet_grid(~Fission) +
  size_colours +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  tau_lab + labs(y = "Adjusted Rand \n Index") +
  guides(colour = guide_legend(nrow = 2))
save_fig(plot_illustration(cl_ideal, c("#294122", "#EB3D00")) +
           (ari / plot_rejection(ideal, "Statistical Power \n 5% level") +
              plot_layout(guides = "collect", heights = c(1, 4))) +
           plot_layout(widths = c(.8, 1)) + plot_annotation(tag_levels = "A") &
           theme(legend.position = "bottom", text = element_text(size = 12), plot.tag = element_text(face = "bold")),
         "FigureS6", 190, 120)
save_fig(plot_illustration(cl_adverse, c("#294122", "#EB3D00", "#FFBBA6")) +
           plot_rejection(read_res("sim01_adverse.rds"), "Type I error rate \n 5% nominal levels") +
           plot_layout(widths = c(.8, 1)) + plot_annotation(tag_levels = "A") &
           theme(text = element_text(size = 12), plot.tag = element_text(face = "bold")),
         "FigureS7", 190, 110)

#------------------------------------------------------------------------------#
# Figure S8: overdispersions estimated within 4 cell populations vs neutrophils

theta <- read_res("app05_overdispersion.rds")
plot_S8 <- function(pop) {
  rmse <- round(sqrt(mean((theta$neutrophil - theta[[pop]])^2)), 2)
  ggplot(theta, 
         aes(neutrophil, .data[[pop]])) +
    scattermore::geom_scattermore(pointsize = 4, alpha = .3) +
    geom_abline(slope = 1, 
                intercept = 0, 
                colour = "darkred", 
                linetype = "dashed", 
                linewidth = 1) +
    facet_grid(~ paste0("RMSE=", rmse)) +
    scale_x_log10(breaks = c(0.01, 0.1, 1, 10), 
                  labels = c(0.01, 0.1, 1, 10)) +
    scale_y_log10() +
    annotation_logticks(sides = "bl") +
    labs(x = "Neutrophil", y = sub("^(.)", "\\U\\1", pop, perl = TRUE))
}
save_fig(wrap_plots(lapply(c("macrophage", "monocyte", "granulocyte", "CD4"), plot_S8),
                    nrow = 1) +
           plot_annotation(tag_levels = "A") &
           theme(text = element_text(size = 14),
                 plot.tag = element_text(face = "bold")),
         "FigureS8", 200, 75)
