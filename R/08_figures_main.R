# Figures 1,2 and 3 of the main manuscript and numbers reported in the text
# (from the intermediate results in results/)

source("R/utils.R")
library(dplyr)
library(ggplot2)
library(patchwork)
library(latex2exp)
theme_set(theme_classic())
dir.create("figures", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)
read_res <- function(file) readRDS(file.path(res_dir, file))

#------------------------------------------------------------------------------#
# Figure 1: NB data thinning in post-clustering inference

# A: illustrative two-component NB mixture, split into 3 clusters by k-means
set.seed(9022024)
X <- cbind(c(rnbinom(50, prob = 0.5, size = 5), rnbinom(50, prob = 0.4, size = 40)),
  rnbinom(100, prob = 0.5, size = 5))
fig1A <- ggplot(data.frame(X1 = X[, 1], X2 = X[, 2], cl = factor(kmeans(X, centers = 3)$cluster))) +
  aes(X1, X2, colour = cl) +
  geom_point(size = 1.5) +
  scale_colour_manual(name = "Clusters", values = c("#294122", "#EB3D00", "#FFBBA6"),
    labels = expression(C[1], C[2], C[3])) +
  labs(x = TeX("$X_1$"), y = TeX("$X_2$"), title = "A") +
  theme(legend.position = "bottom")

# B: QQ-plots of the Wilcoxon p-values on X^(2)
fig1B <- ggplot(readRDS(file.path(res_dir, "sim03_nb_thinning_mixture.rds"))) +
  geom_abline(slope = 1, intercept = 0, col = "red", linewidth = 1.2, alpha = .7) +
  stat_qq(aes(sample = pvalues, colour = Fission), distribution = qunif, size = 1.5) +
  scale_colour_manual(name = "", values = c("#334EAC", "#BAD6EB")) +
  labs(x = "Theoretical Quantiles", y = "Empirical Quantiles", title = "B") +
  lims(x = c(0, 1), y = c(0, 1)) +
  theme(legend.position = "bottom")

# C: Type I error rate as a function of the relative bias in theta
typeI <- readRDS(file.path(res_dir, "sim04_correlation_bias.rds"))$typeI %>%
  mutate(RelativeBias = Error - 1) %>%
  filter(Rho != 0.01, RelativeBias < 5, !(RelativeBias < 0 & Method == "NB")) %>%
  group_by(Method, Estimation, Rho, RelativeBias) %>%
  summarise(typeI = sum(n_reject) / sum(n_tests), .groups = "drop") %>%
  mutate(Label = paste0(ifelse(Method == "NB", "Negative~Binomial", "Gaussian"),
    ifelse(Estimation == "Oracle", "~(theta)", "~(widehat(theta))")))
fig1C <- ggplot(typeI) +
  aes(RelativeBias, typeI, colour = factor(Rho)) +
  geom_point(size = 2) +
  geom_line(linewidth = .9) +
  facet_grid(~Label, labeller = label_parsed, scales = "free_x") +
  geom_hline(yintercept = 0.05, color = "grey20", linewidth = 1.2, linetype = "dashed") +
  scale_colour_manual(name = TeX("Correlation ($\\rho$)"), values = rho_colours) +
  labs(x = TeX("$\\frac{\\widehat{\\theta} - \\theta}{\\theta}$"), y = "Type I error", title = "C") +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 1))

ggsave("figures/Figure1.pdf", (fig1A + fig1B) / fig1C + plot_layout(heights = c(1.2, 1)) &
  theme(plot.title = element_text(face = "bold", size = 12), text = element_text(size = 13)),
width = 185, height = 130, units = "mm")

#------------------------------------------------------------------------------#
# Figure 2: NB data thinning in a homogeneous population of granulocytes

app <- readRDS(file.path(res_dir, "app06_granulocytes.rds"))
genes <- app$genes

fig2A <- ggplot(genes, aes(CorOriginal, CorThinning)) +
  geom_abline(slope = 1, intercept = 0, col = "red", linewidth = 1.2, alpha = .7) +
  geom_point(size = 1, alpha = .5) +
  labs(x = TeX("Cor$\\left(X_1, X_j\\right)$"), y = TeX("Cor$\\left(X^{(1)}_1, X^{(2)}_j\\right)$")) +
  lims(x = c(-1, 1), y = c(-1, 1))

fig2B <- ggplot(tidyr::pivot_longer(genes, c(PvalMultivariate, PvalUnivariate))) +
  geom_abline(slope = 1, intercept = 0, col = "red", linewidth = 1.2, alpha = .7) +
  stat_qq(aes(sample = value, colour = name), distribution = qunif, size = 1.5) +
  scale_colour_manual(name = "Method", values = c("#8C33FF", "#2CA02C"),
    labels = c("Multivariate \n Clustering", "Univariate \n Clustering")) +
  labs(x = "Theoretical Quantiles", y = "Empirical Quantiles") +
  lims(x = c(0, 1), y = c(0, 1))

# C: gene with the smallest p-value after multivariate clustering (PGK1)
g <- which.min(genes$PvalMultivariate)
set.seed(20260122 + g)  # univariate clustering, as in R/06_app_granulocytes.R
cl_univ <- kmeans(log2(app$Xtrain[, g] + 1), centers = 2, nstart = 100)$cluster
df <- data.frame(expr = log2(c(app$Xtrain[, g], app$Xtest[, g]) + 1),
  X = factor(rep(c("Train", "Test"), each = nrow(app$Xtrain)), levels = c("Train", "Test")),
  Multivariate = factor(unname(app$cluster)), Univariate = factor(unname(cl_univ))) %>%
  tidyr::pivot_longer(c(Multivariate, Univariate), names_to = "Method", values_to = "Cluster") %>%
  mutate(Method = paste(Method, "Clustering"))
pvals <- df %>%
  group_by(Method, X) %>%
  summarise(label = paste0("Wilcoxon~italic(p)~\"=\"~\"", scales::pvalue(wilcox_pvalue(expr, Cluster)), "\""),
    .groups = "drop")

# split violins: half violins (left for cluster 1, right for cluster 2)
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, draw_group = function(self, data, panel_params, coord, ...) {
  left <- data$group[1] %% 2 == 1
  data$x <- if (left) data$x - data$violinwidth * (data$x - data$xmin) else data$x + data$violinwidth * (data$xmax - data$x)
  data <- data[order(if (left) data$y else -data$y), ]
  data <- rbind(data[1, ], data, data[nrow(data), ], data[1, ])
  data[c(1, nrow(data) - 1, nrow(data)), "x"] <- round(data[1, "x"])
  GeomPolygon$draw_panel(data, panel_params, coord)
})

fig2C <- ggplot(df, aes(X, expr, fill = Cluster)) +
  layer(stat = "ydensity", geom = GeomSplitViolin, position = "identity", params = list(alpha = .5, trim = FALSE)) +
  geom_boxplot(width = .2, show.legend = FALSE) +
  geom_text(data = pvals, aes(X, max(df$expr) * 1.25, label = label), inherit.aes = FALSE, parse = TRUE, size = 4) +
  scale_fill_manual(name = "Cluster", values = rev(colorspace::lighten(c("#294122", "#EB3D00"), 0.25)),
    labels = expression(C[1], C[2])) +
  scale_x_discrete(labels = c(Train = expression(X^{
    (1)
  }), Test = expression(X^{
    (2)
  }))) +
  facet_grid(~Method) +
  labs(x = "", y = paste(genes$GeneSymbol[g], "log2(counts + 1)"))

ggsave("figures/Figure2.pdf", (fig2A + fig2B) / fig2C + plot_annotation(tag_levels = "A") +
  plot_layout(heights = c(2, 4)) & theme_classic() &
  theme(text = element_text(size = 16), plot.tag = element_text(face = "bold", size = 16)),
width = 250, height = 200, units = "mm")

#------------------------------------------------------------------------------#
# Figure 3 of the main manuscript: conditional vs marginal NB data thinning on
# the mixture of granulocytes and CD24 neutrophils with gene-wise clustering
# (from the results of R/07), and numbers for Section 4.4

res <- readRDS(file.path(res_dir, "app07_granulocytes_neutrophils.rds"))$results

# Numbers for the text: genes with p < 0.05 and ARI quartiles of H1 genes
settings <- unique(res[, c("Thinning", "Clustering")])
numbers <- do.call(rbind, lapply(seq_len(nrow(settings)), FUN = function(i) {
  df <- res[res$Thinning == settings$Thinning[i] & res$Clustering == settings$Clustering[i], ]
  ari_h1 <- quantile(df$ARI[df$Hyp == "H1"])
  return(data.frame(settings[i, ],
    H0_significant = sum(df$pvalue[df$Hyp == "H0"] < 0.05), H0_total = sum(df$Hyp == "H0"),
    H1_significant = sum(df$pvalue[df$Hyp == "H1"] < 0.05), H1_total = sum(df$Hyp == "H1"),
    H1_ARI_min = ari_h1[[1]], H1_ARI_q1 = ari_h1[[2]], H1_ARI_median = ari_h1[[3]],
    H1_ARI_q3 = ari_h1[[4]], H1_ARI_max = ari_h1[[5]]))
}))
write.table(numbers, "results/Figure3_numbers.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(numbers, row.names = FALSE)

# Plotted data: gene-wise clustering only, p-values floored at the smallest
# positive double before the -log10 transform
uni <- res[res$Clustering == "Univariate", ]
message(sum(uni$pvalue < .Machine$double.xmin), " p-value(s) floored at .Machine$double.xmin")
uni$y <- -log10(pmax(uni$pvalue, .Machine$double.xmin))
panel_titles <- c(Conditional = "Conditional thinning\n(oracle: uses the true cell types)",
  Marginal = "Marginal thinning\n(feasible in practice)")
uni$Thinning <- factor(panel_titles[uni$Thinning], levels = panel_titles)
uni$Hyp <- factor(uni$Hyp, levels = c("H1", "H0"))

# In-panel counts
counts <- numbers[numbers$Clustering == "Univariate", ]
annotations <- data.frame(
  Thinning = factor(panel_titles[counts$Thinning], levels = panel_titles),
  de = sprintf("differentially expressed:\n%d/%d detected", counts$H1_significant, counts$H1_total),
  null = sprintf("no difference:\n%d/%d significant", counts$H0_significant, counts$H0_total),
  de_y = c(7, 3.3),  # right of the H1 genes in each panel
  null_y = c(4, 180)  # above the H0 genes, away from the labels
)

# Top 5 genes with the smallest p-values in each panel
top5 <- do.call(rbind, lapply(split(uni, uni$Thinning), FUN = function(df) {
  return(head(df[order(df$pvalue), ], 5))
}))

figure3 <- ggplot(uni, aes(ARI, y, colour = ARI, shape = Hyp)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey30") +
  geom_point(size = 2, alpha = 0.5) +
  # Differentially expressed genes (regular font)
  geom_text(
    data = annotations,
    aes(x = 1.05, y = de_y, label = "differentially expressed:"),
    inherit.aes = FALSE,
    hjust = 1, vjust = 0,
    fontface = "bold",
    size = 2.8,
    colour = "grey25"
  ) +
  geom_text(
    data = annotations,
    aes(x = 1.05, y = de_y, label = sprintf("%d/%d detected",
      counts$H1_significant,
      counts$H1_total)
    ),
    inherit.aes = FALSE,
    hjust = 1, vjust = 0,
    size = 2.8,
    colour = "grey25",
    nudge_y = -0.3
  ) +
  # Bold first line
  geom_text(
    data = annotations,
    aes(x = 0.04, y = null_y, label = "no difference:"),
    inherit.aes = FALSE,
    hjust = 0, vjust = 0.2,
    fontface = "bold",
    size = 2.8,
    colour = "grey25"
  ) +

  # Regular second line
  geom_text(
    data = annotations,
    aes(
      x = 0.04,
      y = null_y,
      label = sprintf("%d/%d significant",
        counts$H0_significant,
        counts$H0_total)
    ),
    inherit.aes = FALSE,
    hjust = 0, vjust = 0,
    size = 2.8,
    colour = "grey25",
    nudge_y = -0.35
  ) +
  ggrepel::geom_label_repel(data = top5, aes(label = GeneSymbol), size = 2.5, fill = alpha("white", 0.8),
    min.segment.length = 0, seed = 1, show.legend = FALSE) +
  facet_wrap(~Thinning) +
  scale_shape_manual(values = c(H1 = 16, H0 = 4),
    labels = c(H1 = sprintf("True differentially expressed genes (n=%d)", counts$H1_total[1]),
      H0 = sprintf("Genes without true difference (n=%d)", counts$H0_total[1]))) +
  # scale_colour_viridis_c(option = "plasma", direction = -1, limits = c(min(res$ARI), 1),
  #  guide = guide_colourbar(barwidth = unit(3, "cm"), barheight = unit(0.3, "cm"))) +
  scale_colour_gradientn(
    colours = MetBrewer::met.brewer("Greek", direction = -1),
    limits = c(min(uni$ARI), 1), breaks = seq(0, 1, by = 0.25),
    labels = c("0", "0.25", "0.5", "0.75", "1"),
    oob = scales::squish) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.1, base = 10),
    breaks = c(0, 1.3, 10, 100), labels = expression(1, "5%", 10^-10, 10^-100),
    minor_breaks = outer(1:9, 10^(-1:2))) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = expression("ARI: gene-wise univariate clusters on" ~ X^{
    (1)
  } ~ "against reference true cell types"),
  y = expression(atop("p-value on" ~ X^{
    (2)
  },
  "(pseudo-log scale)")), shape = NULL) +
  guides(shape = guide_legend(override.aes = list(size = 2.5, alpha = 0.8, colour = "grey40"),
    position = "bottom", direction = "vertical"),
  colour = guide_colourbar(theme = theme_sub_legend(
    key.width = unit(3.5, "cm"),
    title.position = "top",
    margin = margin_auto(NA, 0)  # remove left/right margins
  ), position = "bottom")) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom", # strip.background = element_blank(),
    strip.text = element_text(hjust = 0, size = 12),
    legend.title = element_text(hjust = 0.5, size = 9.5),
    legend.text = element_text(size = 7.75),
    panel.grid.minor.y = element_line(colour = "grey90", linewidth = 0.25),
    axis.text = element_text(size = 7.5)
  )

figure3
ggsave("figures/Figure3.pdf", figure3, width = 8, height = 3.5)
performances <- res %>%
  group_by(Thinning, Clustering) %>%
  summarise(TypeI = mean(pvalue[Hyp == "H0"] < 0.05), Power = mean(pvalue[Hyp == "H1"] < 0.05), .groups = "drop")



#------------------------------------------------------------------------------#
# Numbers reported in the text (Sections 2.2.1 and 4.4, Appendix A7)

app7 <- readRDS(file.path(res_dir, "app07_granulocytes_neutrophils.rds"))
numbers <- c(granulocytes_TSP14 = app$n_cells[["TSP14"]],
  granulocytes_kept_homogeneous = app$n_cells[["kept"]],
  cd24_neutrophils_TSP14 = app7$n_cells[["cd24 neutrophil"]],
  granulocytes_kept_mixture = app7$n_cells[["kept.granulocyte"]],
  cd24_neutrophils_kept_mixture = app7$n_cells[["kept.cd24 neutrophil"]],
  typeI_multivariate_clustering = mean(genes$PvalMultivariate < 0.05),
  typeI_univariate_clustering = mean(genes$PvalUnivariate < 0.05))
write.csv(data.frame(Quantity = names(numbers), Value = numbers), "tables/Numbers_in_text.csv", row.names = FALSE)
