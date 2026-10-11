# PCA analysis for M hapla RILs

library(tidyverse)

# Setting working directory
setwd("~/Desktop/PhD/04_projects/04_mhapla_infection_intensity")
source("00_R/plot_utils.R")


# Loading infection intensity and gene expression data
infection = read.csv("03_infection_intensity/infectionIntensity.csv")
expr = read.csv("04_feature_counts/hapla_log2normCPM.csv")


# Preprocess the data
geneids = expr[,1]
# preparing the expression data for PC analysis
expr = expr[,2:ncol(expr)]

# Removing outliers
expr = expr[,!colnames(expr) %in% c("EA03",
                                    "EA32",
                                    "FB17",
                                    "EA02")]

# Transpose the expression matrix
dt = t(expr)
dt = as.data.frame(dt)

# Scale the data
dt_scaled = scale(dt)

# Run Principal Component Analysis
pca_result <- prcomp(dt_scaled, center = TRUE, scale. = TRUE)
sum = summary(pca_result)

# Extract variance explained
var_exp = sum$importance
var_exp = t(var_exp)
var_exp = as.data.frame(var_exp)
write.csv(var_exp,
          "05_PCA_hapla/pc_variance_explained.csv")

var_exp$PC = 1:93

ggplot(var_exp,
       aes(x = PC,
           y = `Proportion of Variance`)) +
  geom_point(fill = "#533943") +
  theme_tyler() +
  labs(x = "")
ggsave("05_PCA_hapla/proportion_variance_explained.pdf",
       height = 3, width = 4)

ggplot(var_exp,
       aes(x = PC,
           y = `Cumulative Proportion`)) +
  geom_point(fill = "#533943") +
  theme_tyler() 
ggsave("05_PCA_hapla/cumulative_variance_explained.pdf",
       height = 3, width = 4)


# Extract PCA scores
pca_scores <- as.data.frame(pca_result$x)

# Add sample names from rownames
pca_scores$sample <- rownames(pca_scores)

# Add infection intensity data
pca_scores <- pca_scores %>%
  left_join(infection, by = "sample")

# Percent variance explained
pc_var <- summary(pca_result)$importance["Proportion of Variance", ] * 100

ggplot(pca_scores, aes(x = PC1, y = PC2, color = infection_intensity)) +
  geom_point(size = 3, alpha = 0.6) +
  scale_color_gradientn(colors = wes_palette("Zissou1",
                                          100,
                                          "continuous")) +
  labs(
    x = paste0("PC1 (", round(pc_var["PC1"], 1), "%)"),
    y = paste0("PC2 (", round(pc_var["PC2"], 1), "%)"),
    color = "Infection\nintensity"
  ) +
  theme_tyler()

ggplot(pca_scores, aes(x = PC1, y = PC2, color = infection_intensity)) +
  geom_point(size = 3) +
  scale_color_gradient(
    low = "#F2EFE8",   # creamy off-white / warm grey
    high = "#593745"   # muted plum purple
  ) +
  labs(
    x = paste0("PC1 (", round(pc_var["PC1"], 1), "%)"),
    y = paste0("PC2 (", round(pc_var["PC2"], 1), "%)"),
    color = "Infection\nintensity"
  ) +
  theme_tyler()

ggsave("05_PCA_hapla/PCA_hapla.png",
       dpi = 600,
       height = 4, width = 6)

ggsave("05_PCA_hapla/PCA_hapla.png",
       dpi = 600,
       height = 4, width = 5)

# Expression of female marker gene MhA1_Contig298.frz3.gene8
colnames(dt) = geneids
female_marker = dt[,colnames(dt) %in% "MhA1_Contig298.frz3.gene8_mRNA"]
female_marker = cbind(female_marker, row.names(dt))
colnames(female_marker) = c("female_marker",
                            "sample")
female_marker = as.data.frame(female_marker)

pca_scores <- pca_scores %>%
  left_join(female_marker, by = "sample")
pca_scores$female_marker = as.numeric(pca_scores$female_marker)

ggplot(pca_scores, aes(x = PC1, y = PC2, color = female_marker)) +
  geom_point(size = 3) +
  scale_color_gradient(
    low = "#F1EFEA",   # creamy off-white / warm grey
    high = "#8A4F3D"   # muted plum purple
  ) +
  labs(
    x = paste0("PC1 (", round(pc_var["PC1"], 1), "%)"),
    y = paste0("PC2 (", round(pc_var["PC2"], 1), "%)"),
    color = " "
  ) +
  theme_tyler()

ggsave("05_PCA_hapla/female_marker_PCA.png",
       dpi = 600,
       height = 2, width = 4)

# Expression of J4 gene MhA1_Contig2096.frz3.gene8
j4 = dt[,colnames(dt) %in% "MhA1_Contig2096.frz3.gene8_mRNA"]
j4 = cbind(j4, row.names(dt))
colnames(j4) = c("j4_marker",
                            "sample")
j4 = as.data.frame(j4)

pca_scores <- pca_scores %>%
  left_join(j4, by = "sample")
pca_scores$j4_marker = as.numeric(pca_scores$j4_marker)

ggplot(pca_scores, aes(x = PC1, y = PC2, color = j4_marker)) +
  geom_point(size = 3) +
  scale_color_gradient(
    low = "#F1EFEA",   # creamy off-white / warm grey
    high = "#3E5F73"   # muted plum purple
  ) +
  labs(
    x = paste0("PC1 (", round(pc_var["PC1"], 1), "%)"),
    y = paste0("PC2 (", round(pc_var["PC2"], 1), "%)"),
    color = ""
  ) +
  theme_tyler() +
  theme(axis.title.x = element_blank())


ggsave("05_PCA_hapla/j4_marker_PCA.png",
       dpi = 600,
       height = 2, width = 4)


# Expression of DPY18 (MhA1_Contig2702.frz3.gene1)
colnames(dt) = geneids
dpy18 = dt[,colnames(dt) %in% "MhA1_Contig2702.frz3.gene1_mRNA"]
dpy18 = cbind(dpy18, row.names(dt))
colnames(dpy18) = c("dpy18",
                            "sample")
dpy18 = as.data.frame(dpy18)

pca_scores <- pca_scores %>%
  left_join(dpy18, by = "sample")
pca_scores$dpy18 = as.numeric(pca_scores$dpy18)

ggplot(pca_scores, aes(x = PC1, y = PC2, color = dpy18)) +
  geom_point(size = 3) +
  scale_color_gradient(
    low = "#F1EFEA",   # creamy off-white / warm grey
    high = "#762A83",   # muted plum purple
    limits = c(5,7),
    oob = scales::squish
  ) +
  labs(
    x = paste0("PC1 (", round(pc_var["PC1"], 1), "%)"),
    y = paste0("PC2 (", round(pc_var["PC2"], 1), "%)"),
    color = " "
  ) +
  theme_tyler()

ggsave("05_PCA_hapla/DPY18_PCA.png",
       dpi = 600,
       height = 4, width = 5)


# Expression of ALG2 (-	MhA1_Contig684.frz3.gene16 )
alg2 = dt[,colnames(dt) %in% "MhA1_Contig684.frz3.gene16_mRNA"]
alg2 = cbind(alg2, row.names(dt))
colnames(alg2) = c("alg2",
                    "sample")
alg2 = as.data.frame(alg2)

pca_scores <- pca_scores %>%
  left_join(alg2, by = "sample")
pca_scores$alg2 = as.numeric(pca_scores$alg2)

ggplot(pca_scores, aes(x = PC1, y = PC2, color = alg2)) +
  geom_point(size = 3) +
  scale_color_gradient(
    low = "#F1EFEA",   # creamy off-white / warm grey
    high = "#4F8AA3",   # muted plum purple
   # limits = c(5,7),
    oob = scales::squish
  ) +
  labs(
    x = paste0("PC1 (", round(pc_var["PC1"], 1), "%)"),
    y = paste0("PC2 (", round(pc_var["PC2"], 1), "%)"),
    color = " "
  ) +
  theme_tyler()

ggsave("05_PCA_hapla/ALG2_PCA.png",
       dpi = 600,
       height = 4, width = 5)


# Expression of ERP3 (	MhA1_Contig593.frz3.gene3  )
erp3 = dt[,colnames(dt) %in% "MhA1_Contig593.frz3.gene3_mRNA"]
erp3 = cbind(erp3, row.names(dt))
colnames(erp3) = c("erp3",
                   "sample")
erp3 = as.data.frame(erp3)

pca_scores <- pca_scores %>%
  left_join(erp3, by = "sample")
pca_scores$erp3 = as.numeric(pca_scores$erp3)

ggplot(pca_scores, aes(x = PC1, y = PC2, color = erp3)) +
  geom_point(size = 3) +
  scale_color_gradient(
    low = "#F1EFEA",   # creamy off-white / warm grey
    high = "#C27A43",   # muted plum purple
    # limits = c(5,7),
    oob = scales::squish
  ) +
  labs(
    x = paste0("PC1 (", round(pc_var["PC1"], 1), "%)"),
    y = paste0("PC2 (", round(pc_var["PC2"], 1), "%)"),
    color = " "
  ) +
  theme_tyler()

ggsave("05_PCA_hapla/ERP3_PCA.png",
       dpi = 600,
       height = 4, width = 5)

## PLOTTING MODULE EIGENGENES

library(tidyverse)
library(patchwork)
library(scales)

MEs = read.csv("coXnetwork/moduleEigen.csv")

# ------------------------------------------------------------
# Convert module eigengenes to long format
# ------------------------------------------------------------

eigengene_long <- MEs %>%
  as.data.frame() %>%
  pivot_longer(
    cols = -c(sample, X),
    names_to = "module",
    values_to = "eigengene"
  )

module_pca <- pca_scores %>%
  select(sample, PC1, PC2, infection_intensity) %>%
  left_join(eigengene_long, by = "sample")

module_pca = left_join(module_pca, parasite_modules,
                       by = "module")

module_pca = module_pca[,c(1:4, 6, 7, 12)]

modules_to_plot <- c(
  "parasite7",
  "parasite6",
  "parasite8",
  "parasite2",
  "parasite1",
  "parasite3",
  "parasite9"
)

module_pca_plot <- module_pca %>%
  filter(name %in% modules_to_plot) %>%
  mutate(
    name = factor(name, levels = modules_to_plot)
  )

eigengene_limit <- max(
  abs(module_pca_plot$eigengene),
  na.rm = TRUE
)

ggplot(
  module_pca_plot,
  aes(x = PC1, y = PC2, color = eigengene)
) +
  geom_point(
    size = 2.2,
    alpha = 0.9
  ) +
  facet_wrap(
    ~ name,
    ncol = 3
  ) +
  scale_color_gradientn(
    colors = interaction_palette,
    limits = c(-eigengene_limit, eigengene_limit),
    values = scales::rescale(seq(-eigengene_limit,
                                 eigengene_limit,
                                 length.out = length(interaction_palette))),
    oob = scales::squish,
    name = "Module\neigengene"
  ) +
  labs(
    x = paste0("PC1 (", round(pc_var["PC1"], 1), "%)"),
    y = paste0("PC2 (", round(pc_var["PC2"], 1), "%)")
  ) +
  theme_tyler() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 10
    ),
    panel.spacing = unit(0.7, "lines"),
    legend.title = element_text(
      face = "bold"
    ),
    legend.position = "right"
  )

ggsave("09_coexpression_hapla/sig_modules_eigen_PCA.pdf",
       width = 10, height = 6,
       bg = "transparent", dpi = 600)
