#### Principal Component Analysis with Medicago Gene Expression ####
## Tyler Inkseep, Groen Lab, UC Riverside ##
# 1 July 2025 #

library(tidyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(princurve)
library(wesanderson)

### SETWD ###
setwd("~/Desktop/PhD/04_projects/04_mhapla_infection_intensity")

source("00_R/plot_utils.R")

### IMPORT DATASETS ###
# Import normalized counts
counts = read.csv("04_feature_counts/parsed_mtrunc_log2normCPM.csv")
counts$geneid = gsub("rna-", "", counts$X)
geneids = counts$geneid

# Transpose so samples are rows
geneids = counts$geneid
expr_t <- t(counts[,2:84])
expr_t <- expr_t[!rownames(expr_t) %in% c("EA02",
                                          "EA32",
                                          "FB17"),] # 80 genotypes

pca_res <- prcomp(expr_t, scale. = TRUE, center = TRUE)

# Getting loadings for genes into each PC
loadings_df <- as.data.frame(pca_res$rotation)
loadings_df$geneid = geneids

loadings_df_clean <- loadings_df %>%
  distinct(geneid, .keep_all = TRUE)

# Grabbing the top 5% of genes into each PC & pivocitng the loading df into tidy format
top_pct <- 0.05  # Top 5%
loading_long <- loadings_df_clean %>%
  pivot_longer(-geneid, names_to = "PC", values_to = "Loading") %>%
  group_by(PC) %>%
  mutate(
    absLoading = abs(Loading)
  ) %>%
  arrange(desc(absLoading), .by_group = TRUE) %>%
  mutate(
    rank = row_number(),
    n_top = ceiling(n() * top_pct),
    TopGene = rank <= n_top
  ) %>%
  ungroup()


# Getting variance explained
var_explained <- pca_res$sdev^2 / sum(pca_res$sdev^2)

# Plotting the PC
pca_scores <- as.data.frame(pca_res$x) %>%
  tibble::rownames_to_column("sample")

# Loading the infection intensity score
infection_intensity = read.csv("03_infection_intensity/infectionIntensity.csv")

pca_scores <- pca_scores %>%
  left_join(infection_intensity, by = "sample")

plot1 = ggplot(pca_scores, aes(x = PC1, y = PC2,
                                color = infection_intensity)) +
  geom_point(size = 3) +
  scale_color_gradient(
    low = "#F2EFE8",   # creamy off-white / warm grey
    high = "#3A6547"   # muted plum purple
  ) +
  labs(x = "PC1 (16.48% variance)",
       y = "PC2 (14.03% variance)",
       color = "Infection\nIntensity") +
  theme_tyler()
plot1
ggsave("06_PCA_medicago/pc_infection_intensity_medicago.pdf",
       plot1,
       height = 4, width = 5,
       bg = "transparent",
       dpi = 300)

plot1 = ggplot(pca_scores, aes(x = PC3, y = PC4,
                               color = infection_intensity)) +
  geom_point(size = 3) +
  scale_color_gradient(
    low = "#F2EFE8",   # creamy off-white / warm grey
    high = "#3A6547"   # muted plum purple
  ) +
  labs(x = "PC3 (8.74% variance)",
       y = "PC4 (6.99% variance)",
       color = "Infection\nIntensity") +
  theme_tyler()
plot1
ggsave("06_PCA_medicago/pc_infection_intensity_medicago_pc3_pc4.pdf",
       plot1,
       height = 6, width = 5,
       bg = "transparent",
       dpi = 300)

sum = summary(pca_res)

# Extract variance explained
var_exp = sum$importance
var_exp = t(var_exp)
var_exp = as.data.frame(var_exp)

var_exp$PC = 1:80

ggplot(var_exp,
       aes(x = PC,
           y = `Proportion of Variance`)) +
  geom_point(color = "#3A6547") +
  theme_tyler() +
  labs(x = "")
ggsave("06_PCA_medicago/proportion_variance_explained.pdf",
       height = 3, width = 4)

ggplot(var_exp,
       aes(x = PC,
           y = `Cumulative Proportion`)) +
  geom_point(color = "#3A6547") +
  theme_tyler() 
ggsave("06_PCA_medicago/cumulative_variance_explained.pdf",
       height = 3, width = 4)


# Top PCs
PC_cols <- grep("^PC", colnames(pca_scores), value = TRUE)

PC_model_results <- lapply(PC_cols, function(mod) {
  
  form <- as.formula(paste(mod, "~ infection_intensity"))
  fit <- lm(form, data = pca_scores)
  
  broom::tidy(fit) %>%
    filter(term == "infection_intensity") %>%
    mutate(
      Module = gsub("^ME", "", mod),
      Eigengene = mod
    )
  
}) %>%
  bind_rows() %>%
  mutate(
    p_adj = p.adjust(p.value, method = "BH"),
    neglog10_p = -log10(p.value),
    neglog10_padj = -log10(p_adj),
    direction = ifelse(estimate > 0, "Positive", "Negative")
  )

PC_model_results$Module = gsub("PC", "", PC_model_results$Module)
PC_model_results = PC_model_results[which(PC_model_results$Module %in% c(1:10)),]

plot_PC_infection <- PC_model_results %>%
  mutate(
    Module = reorder(Module, estimate)
  ) %>%
  ggplot(aes(x = Module, y = estimate, fill = neglog10_padj)) +
  geom_col(width = 0.75) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  scale_fill_gradient(
    low = "#D8D8D8",
    high = "#3A6547",
    name = expression(-log[10](p))
  ) +
  labs(
    x = "Principal Component",
    y = "β (Infection Intensity)"
  ) +
  coord_flip() +
  theme_tyler() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

plot_PC_infection

ggsave("06_PCA_medicago/PC_beta_1_10.png",
       plot_PC_infection,
       dpi = 600,
       height = 6, width = 4)

# Multivariate selection model on PC scores
fit <- lm(infection_intensity ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11, data = pca_scores)
sum = summary(fit)
model = sum$coefficients

# PCA scores
pca_scores <- as.data.frame(pca_res$x)
pca_scores$sample <- rownames(pca_scores)

# Make sure infection_intensity has sample IDs
# Example assumes infection_intensity has columns: sample and infection_intensity
dat <- merge(pca_scores, infection_intensity, by = "sample")

# Check that merge worked
dim(pca_scores)
dim(infection_intensity)
dim(dat)

# Model
fit <- lm(infection_intensity ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11,
          data = dat)

summary(fit)

library(dplyr)
library(ggplot2)
library(broom)

pc_cols <- paste0("PC", 1:11)

dat_scaled <- dat %>%
  mutate(
    infection_intensity_z = as.numeric(scale(infection_intensity)),
    across(all_of(pc_cols), ~ as.numeric(scale(.x)))
  )

fit_scaled <- lm(infection_intensity_z ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11,
                 data = dat_scaled)

coef_scaled <- tidy(fit_scaled) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    PC = gsub("PC", "", term),
    PC = factor(PC, levels = as.character(1:11)),
    neglog10p = -log10(p.value)
  )

ggplot(coef_scaled, aes(x = estimate, y = PC, fill = neglog10p)) +
  geom_col(width = 0.75) +
  geom_vline(xintercept = 0, linewidth = 0.7) +
  scale_fill_gradient(
    low = "grey80",
    high = "#214f32",
    name = expression(-log[10](p))
  ) +
  labs(
    x = "β, Infection Intensity",
    y = "Principal Component"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

ggsave("06_PCA_medicago/PC_beta_1_10.png",
       dpi = 600,
       height = 6, width = 4)


medicago_genes = read.csv("medicagoGeneAnnot.csv")

pc3 = loading_long[which(loading_long$PC == "PC3"),]
pc3 = pc3[which(pc3$TopGene == TRUE),]
pc3 = left_join(pc3, medicago_genes,
                by = c("geneid" = "NCBIGeneID"))
write.csv(pc3,
          "06_PCA_medicago/PC3_genes.csv")


pc4 = loading_long[which(loading_long$PC == "PC4"),]
pc4 = pc4[which(pc4$TopGene == TRUE),]
pc4 = left_join(pc4, medicago_genes,
                by = c("geneid" = "NCBIGeneID"))
write.csv(pc4,
          "06_PCA_medicago/PC4_genes.csv")

library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(clusterProfiler)

# ----------------------------
# 1. Read gene lookup table
# ----------------------------

lookup_raw <- read_excel("00_source_data/medicago_gene_lookup.xlsx",
                         sheet = 2)

lookup_long <- lookup_raw %>%
  select(
    Medtr_ID = Mt4.0v1_ID,
    Mtrun_ID = `MtrunA17r5.0-ANR_ID`
  ) %>%
  filter(!is.na(Medtr_ID), !is.na(Mtrun_ID)) %>%
  separate_rows(Mtrun_ID, sep = "\\|") %>%
  mutate(
    Medtr_ID = str_trim(Medtr_ID),
    Mtrun_ID = str_trim(Mtrun_ID)
  ) %>%
  distinct(Mtrun_ID, Medtr_ID)

head(lookup_long)

lines <- readLines("00_source_data/medicago_go.txt")

go_raw <- read_csv(
  "00_source_data/medicago_go.txt",
  show_col_types = FALSE
)

go_raw <- go_raw %>%
  rename(
    Medtr_ID = `Gene Name`,
    GO = `GO ID`,
    GO_description = `GO Description`
  )

go_clean <- go_raw %>%
  mutate(
    Medtr_ID = str_trim(Medtr_ID),
    GO = str_trim(GO),
    GO_description = str_trim(GO_description)
  ) %>%
  filter(!is.na(GO), GO != "", str_detect(GO, "^GO:")) %>%
  distinct(Medtr_ID, GO, GO_description)

term2gene <- go_clean %>%
  distinct(GO, Medtr_ID)

term2name <- go_clean %>%
  distinct(GO, GO_description)

library(dplyr)
library(stringr)
library(clusterProfiler)

# ----------------------------
# 1. Map PCA/loadings genes to Medtr IDs
# ----------------------------

loadings_mapped <- loadings_df_clean %>%
  left_join(medicago_genes, by = c("geneid" = "NCBIGeneID"))
  
loadings_mapped <- loadings_mapped %>%
  left_join(lookup_long, by = c("GENEID" = "Mtrun_ID"))


# Check mapping success
mapping_summary <- loadings_mapped %>%
  summarise(
    total_pca_genes = n(),
    mapped_to_medtr = sum(!is.na(Medtr_ID)),
    mapping_rate = mean(!is.na(Medtr_ID))
  )

mapping_summary

# all PCA genes that mapped to Medtr IDs
background_genes <- loadings_mapped %>%
  filter(!is.na(Medtr_ID)) %>%
  pull(Medtr_ID) %>%
  unique()

# restrict to genes that actually have GO annotations
background_genes <- intersect(background_genes, unique(term2gene$Medtr_ID))

length(background_genes) #10,248

n_top <- ceiling(nrow(loadings_mapped) * 0.05)  # top 0.5%

pc3_pos <- loadings_mapped %>%
  arrange(desc(PC3)) %>%
  slice_head(n = n_top) %>%
  filter(!is.na(Medtr_ID)) %>%
  pull(Medtr_ID) %>%
  unique() %>%
  intersect(background_genes)

pc3_neg <- loadings_mapped %>%
  arrange(PC3) %>%
  slice_head(n = n_top) %>%
  filter(!is.na(Medtr_ID)) %>%
  pull(Medtr_ID) %>%
  unique() %>%
  intersect(background_genes)

pc4_pos <- loadings_mapped %>%
  arrange(desc(PC4)) %>%
  slice_head(n = n_top) %>%
  filter(!is.na(Medtr_ID)) %>%
  pull(Medtr_ID) %>%
  unique() %>%
  intersect(background_genes)

pc4_neg <- loadings_mapped %>%
  arrange(PC4) %>%
  slice_head(n = n_top) %>%
  filter(!is.na(Medtr_ID)) %>%
  pull(Medtr_ID) %>%
  unique() %>%
  intersect(background_genes)

data.frame(
  gene_set = c("PC3 positive", "PC3 negative", "PC4 positive", "PC4 negative"),
  n_genes = c(length(pc3_pos), length(pc3_neg), length(pc4_pos), length(pc4_neg))
)

library(clusterProfiler)

ego_pc3_pos <- enricher(
  gene = pc3_pos,
  universe = background_genes,
  TERM2GENE = term2gene,
  TERM2NAME = term2name,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.10
)

ego_pc3_neg <- enricher(
  gene = pc3_neg,
  universe = background_genes,
  TERM2GENE = term2gene,
  TERM2NAME = term2name,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.10
)

ego_pc4_pos <- enricher(
  gene = pc4_pos,
  universe = background_genes,
  TERM2GENE = term2gene,
  TERM2NAME = term2name,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.10
)

ego_pc4_neg <- enricher(
  gene = pc4_neg,
  universe = background_genes,
  TERM2GENE = term2gene,
  TERM2NAME = term2name,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.10
)

pc3_pos_go <- as.data.frame(ego_pc3_pos)
pc3_neg_go <- as.data.frame(ego_pc3_neg)
pc4_pos_go <- as.data.frame(ego_pc4_pos)
pc4_neg_go <- as.data.frame(ego_pc4_neg)

head(pc3_pos_go, 20)
head(pc3_neg_go, 20)
head(pc4_pos_go, 20)
head(pc4_neg_go, 20)

# GSEA
library(dplyr)
library(clusterProfiler)

# Keep mapped genes with GO background
loadings_gsea <- loadings_mapped %>%
  filter(!is.na(Medtr_ID)) %>%
  filter(Medtr_ID %in% background_genes)

collapse_rank <- function(df, pc_col) {
  df %>%
    select(Medtr_ID, loading = all_of(pc_col)) %>%
    filter(!is.na(loading)) %>%
    group_by(Medtr_ID) %>%
    slice_max(order_by = abs(loading), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc(loading)) %>%
    { setNames(.$loading, .$Medtr_ID) }
}

rank_pc3 <- collapse_rank(loadings_gsea, "PC3")
rank_pc4 <- collapse_rank(loadings_gsea, "PC4")

gsea_pc3 <- GSEA(
  geneList = rank_pc3,
  TERM2GENE = term2gene,
  TERM2NAME = term2name,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 1,
  verbose = FALSE
)

gsea_pc4 <- GSEA(
  geneList = rank_pc4,
  TERM2GENE = term2gene,
  TERM2NAME = term2name,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 1,
  verbose = FALSE
)

gsea_pc3_df <- as.data.frame(gsea_pc3)
gsea_pc4_df <- as.data.frame(gsea_pc4)

head(gsea_pc3_df)
head(gsea_pc4_df)

library(ggplot2)
library(forcats)

gsea_all <- bind_rows(
  gsea_pc3_df %>% mutate(PC = "PC3"),
  gsea_pc4_df %>% mutate(PC = "PC4")
) %>%
  mutate(
    infection_side = case_when(
      NES < 0 ~ "Higher infection side",
      NES > 0 ~ "Lower infection side"
    ),
    neglog10_fdr = -log10(p.adjust)
  )

gsea_plot <- gsea_all %>%
  filter(p.adjust < 0.10) %>%
  group_by(PC) %>%
  slice_min(order_by = p.adjust, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(Description = fct_reorder(Description, NES))

ggplot(gsea_plot, aes(x = NES, y = Description)) +
  geom_vline(xintercept = 0, linewidth = 0.7) +
  geom_point(aes(size = setSize, color = neglog10_fdr)) +
  scale_color_gradient(
    low = "#F2EFE8",   # creamy off-white / warm grey
    high = "#3A6547"   # deep green
  ) +
  facet_wrap(~ PC, scales = "free_y") +
  labs(
    x = "Normalized enrichment score",
    y = NULL,
    size = "Gene set size",
    color = expression(-log[10]("FDR"))
  ) +
  theme_tyler()

pc3_go_plot = gsea_pc3_df[which(gsea_pc3_df$p.adjust <= 0.01),]
pc3_go_plot$neglog10_fdr = -log10(pc3_go_plot$p.adjust)
pc3_go_plot = pc3_go_plot %>%
  arrange(NES)

pc3_go_plot <- gsea_pc3_df %>%
  filter(p.adjust <= 0.01) %>%
  mutate(
    neglog10_fdr = -log10(p.adjust),
    Description = fct_reorder(Description, NES)
  )

ggplot(pc3_go_plot, aes(x = NES, y = Description)) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  geom_point(aes(size = setSize, color = neglog10_fdr)) +
  scale_color_gradient(
    low = "#F2EFE8",
    high = "#3A6547"
  ) +
  labs(
    x = "Normalized enrichment score",
    y = NULL,
    size = "Gene set size",
    color = expression(-log[10]("FDR"))
  ) +
  theme_tyler()

ggsave("06_PCA_medicago/pc3_gsea.pdf",
       height = 4, width = 8)

pc4_go_plot = gsea_pc3_df[which(gsea_pc4_df$p.adjust <= 0.01),]
pc4_go_plot$neglog10_fdr = -log10(pc4_go_plot$p.adjust)
pc4_go_plot = pc4_go_plot %>%
  arrange(NES)

pc4_go_plot <- gsea_pc4_df %>%
  filter(p.adjust <= 0.01) %>%
  mutate(
    neglog10_fdr = -log10(p.adjust),
    Description = fct_reorder(Description, NES)
  )

ggplot(pc4_go_plot, aes(x = NES, y = Description)) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  geom_point(aes(size = setSize, color = neglog10_fdr)) +
  scale_color_gradient(
    low = "#F2EFE8",
    high = "#3A6547"
  ) +
  labs(
    x = "Normalized enrichment score",
    y = NULL,
    size = "Gene set size",
    color = expression(-log[10]("FDR"))
  ) +
  theme_tyler()

ggsave("06_PCA_medicago/pc4_gsea.pdf",
       height = 5, width = 9)

## OTHER FIGURES
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

## Rebuild PC scores in case pca_scores was overwritten later
pca_scores_plot <- as.data.frame(pca_res$x) %>%
  tibble::rownames_to_column("sample")

## Standardize the expression-table gene IDs
counts_plot <- counts %>%
  mutate(
    expression_geneid = X %>%
      str_remove("^rna-") %>%
      str_remove("\\.[0-9]+$")
  )


plot_gene_pc34 <- function(
    locus_tag,
    gene_label,
    high_color,
    output_file
) {
  
  ## Find transcript(s) associated with the requested locus tag
  target_transcripts <- host_results_mapped %>%
    filter(LocusTag == locus_tag) %>%
    pull(transcript_id) %>%
    str_remove("^rna-") %>%
    str_remove("\\.[0-9]+$") %>%
    unique()
  
  if (length(target_transcripts) == 0) {
    stop("No transcript found for ", locus_tag)
  }
  
  ## Extract expression and average if multiple transcripts map
  gene_expression <- counts_plot %>%
    filter(expression_geneid %in% target_transcripts) %>%
    dplyr::select(expression_geneid, all_of(rownames(expr_t))) %>%
    pivot_longer(
      cols = -expression_geneid,
      names_to = "sample",
      values_to = "expression"
    ) %>%
    group_by(sample) %>%
    summarise(
      expression = mean(expression, na.rm = TRUE),
      .groups = "drop"
    )
  
  if (nrow(gene_expression) == 0) {
    stop(
      "Transcript found in host_results_mapped but not in counts: ",
      paste(target_transcripts, collapse = ", ")
    )
  }
  
  plot_data <- pca_scores_plot %>%
    left_join(gene_expression, by = "sample")
  
  p <- ggplot(
    plot_data,
    aes(
      x = PC3,
      y = PC4,
      color = expression
    )
  ) +
    geom_point(
      size = 3,
      alpha = 0.95
    ) +
    scale_color_gradient(
      low = "#F2EFE8",
      high = high_color
    ) +
    labs(
      x = paste0(
        "PC3 (",
        round(var_explained[3] * 100, 2),
        "% variance)"
      ),
      y = paste0(
        "PC4 (",
        round(var_explained[4] * 100, 2),
        "% variance)"
      ),
      color = paste0(
        gene_label,
        "\nexpression"
      )
    ) +
    theme_tyler()
  
  ggsave(
    output_file,
    plot = p,
    height = 4,
    width = 5,
    dpi = 600,
    bg = "transparent"
  )
  
  p
}

tctp_pc_plot <- plot_gene_pc34(
  locus_tag = "MtrunA17_Chr1g0192201",
  gene_label = "MtTCTP",
  high_color = "#8A4F3D",
  output_file = "06_PCA_medicago/PC3_PC4_TCTP_expression.pdf"
)

tctp_pc_plot


## PLOTTING MODULE EIGENGENES

library(tidyverse)
library(patchwork)
library(scales)

MEs = read.csv("coXnetwork/medicago/moduleEigen.csv")

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
  select(sample, PC3, PC4, infection_intensity) %>%
  left_join(eigengene_long, by = "sample")

module_pca = left_join(module_pca, host_modules,
                       by = "module")

module_pca = module_pca[,c(1:4, 6, 7, 12)]

modules_to_plot <- c(
  "host26",
  "host19",
  "host18",
  "host21",
  "host14",
  "host16",
  "host28",
  "host20"
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
  aes(x = PC3, y = PC4, color = eigengene)
) +
  geom_point(
    size = 2.2,
    alpha = 0.9
  ) +
  facet_wrap(
    ~ name,
    ncol = 2
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
    x = paste0("PC3 (", round(pc_var["PC3"], 1), "%)"),
    y = paste0("PC4 (", round(pc_var["PC4"], 1), "%)")
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
    legend.position = "bottom"
  )

ggsave("10_coepxression_medicago/sig_modules_eigen_PCA.pdf",
       width = 7, height = 6,
       bg = "transparent", dpi = 600)
