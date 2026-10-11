#### Principal Component Analysis with Medicago Gene Expression - scMarker Genes ####
## Tyler Inkseep, Groen Lab, UC Riverside ##
# 14 July 2025 #

library(tidyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(princurve)
library(wesanderson)
library(readxl)
library(purrr)

### SETWD ###
setwd("~/Desktop/groenlab/haplaRIL")

### IMPORT DATASETS ###
# Import normalized counts
counts = read.csv("featureCounts/parsed_mtrunc_log2normCPM.csv")
counts$geneid = gsub("rna-", "", counts$X)
geneids = counts$geneid

medicagoGeneAnnot <- read_csv("medicagoGeneAnnot.csv")

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
loadings_df = left_join(loadings_df, medicagoGeneAnnot,
                        by = c("geneid" = "NCBIGeneID"))
loadings_df = loadings_df[,c(1:80,83)]
loadings_df$GENEID = gsub("_", "", loadings_df$GENEID)

loadings_df_clean <- loadings_df %>%
  distinct(GENEID, .keep_all = TRUE)

# Grabbing the top 5% of genes into each PC & pivocitng the loading df into tidy format
top_pct <- 0.05  # Top 5%
loading_long <- loadings_df_clean %>%
  pivot_longer(-GENEID, names_to = "PC", values_to = "Loading") %>%
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

# Loading the cell type marker genes from Pereira et al.
# List all sheet names
excel_file <- "PereiraCellRepscRNA.xlsx"
sheets <- excel_sheets(excel_file)
sheets = sheets[7:30]

# Import all sheets into one tidy dataframe
marker_df <- map_dfr(sheets, ~ read_excel(excel_file, sheet = .x))
marker_df$cluster = as.character(marker_df$cluster)

# Add cluster annotations
cluster_annot = 1:24
cluster_annot = cbind(cluster_annot, c("Cortex-1", "Atrichoblast-1", "Cortex-2",
                            "Atrichoblast-2", "Stele-1", "Nodule-1",
                            "Pericycle", "Endodermis", "Nodule-2",
                            "Suberized Endodermis", "Trichoblast", "Unidentified-1",
                            "Phloem", "Dividing Stele", "Nodule-3",
                            "Infected Stele", "Xylem-1", "Stele-2",
                            "Xylem-2", "Unidentified-2", "Lateral Root-1",
                            "Lateral Root-2", "Unidentified-3", "Unidentified-4"))
colnames(cluster_annot) = c("cluster", "ClusterType")
cluster_annot = as.data.frame(cluster_annot)

marker_df = left_join(marker_df, cluster_annot,
                      by = "cluster")

# Merge module data with the loadings & pivot to longer
loadings_df <- loadings_df %>%
  left_join(marker_df, by = c("GENEID" = "gene_id"))

top_pct <- 0.01 

loading_long <- loadings_df %>%
  pivot_longer(cols = starts_with("PC"), names_to = "PC", values_to = "Loading") %>%
  group_by(PC) %>%
  mutate(absLoading = abs(Loading)) %>%
  arrange(desc(absLoading), .by_group = TRUE) %>%
  mutate(
    rank = row_number(),
    n_top = ceiling(n() * top_pct),
    TopGene = rank <= n_top
  ) %>%
  ungroup()

pcs <- unique(loading_long$PC)
clusters <- 1:24

results_list <- list()

for (pc in pcs) {
  for (mod in clusters) {
    df_pc <- loading_long[loading_long$PC == pc, ]
    
    mod_genes <- df_pc$cluster == mod & !is.na(df_pc$cluster)
    non_mod_genes <- is.na(df_pc$cluster) | df_pc$cluster != mod
    
    TopGene <- df_pc$TopGene
    
    TopClusterGenes <- sum(TopGene & mod_genes)
    TopNonClusterGenes <- sum(TopGene & non_mod_genes)
    OtherClusterGenes <- sum(!TopGene & mod_genes)
    OtherNonClusterGenes <- sum(!TopGene & non_mod_genes)
    
    # Initialize odds ratio
    odds <- NA
    
    # Check if valid contingency table
    if ((TopClusterGenes + OtherClusterGenes) > 0 & (TopNonClusterGenes + OtherNonClusterGenes) > 0) {
      test <- fisher.test(matrix(c(TopClusterGenes, TopNonClusterGenes,
                                   OtherClusterGenes, OtherNonClusterGenes), nrow = 2))
      pval <- test$p.value
      odds <- as.numeric(test$estimate)  # odds ratio
    } else {
      pval <- NA
    }
    
    # Save to list
    results_list[[length(results_list) + 1]] <- data.frame(
      PC = pc,
      Cluster = mod,
      TopClusterGenes = TopClusterGenes,
      TopNonClusterGenes = TopNonClusterGenes,
      OtherClusterGenes = OtherClusterGenes,
      OtherNonClusterGenes = OtherNonClusterGenes,
      p_value = pval,
      odds_ratio = odds
    )
  }
}

enrichment_results <- do.call(rbind, results_list)

# Adjust p-values
enrichment_results$p_adj <- p.adjust(enrichment_results$p_value,
                                     method = "BH")
enrichment_results_sub = enrichment_results[which(enrichment_results$PC %in% c("PC1", "PC2",
                                                                               "PC3", "PC4",
                                                                               "PC5", "PC6",
                                                                               "PC7", "PC8",
                                                                               "PC9", "PC10",
                                                                               "PC11", "PC12",
                                                                               "PC13", "PC14",
                                                                               "PC15", "PC16",
                                                                               "PC17", "PC18")),]
enrich_sig = enrichment_results_sub[which(enrichment_results_sub$p_adj <= 0.05),]
enrich_sig$logp = -log10(enrich_sig$p_adj)
enrich_sig$PC = factor(enrich_sig$PC, levels = c("PC1", "PC2",
                                                 "PC3", "PC4",
                                                 "PC5", "PC6",
                                                 "PC7", "PC8",
                                                 "PC9", "PC10",
                                                 "PC11", "PC12",
                                                 "PC13", "PC14",
                                                 "PC15", "PC16",
                                                 "PC17", "PC18"))

plot5 = ggplot(enrich_sig, aes(x = PC, y = Cluster)) +
  geom_point(aes(size = odds_ratio, color = logp)) +
  scale_color_gradientn(colors = wes_palette("Zissou1", 100, "continuous")) +
  scale_size_continuous(range = c(2, 8)) +
  labs(
    x = "Principal Component",
    y = "Cluster",
    color = "-log10(adj p-value)",
    size = "Odds ratio"
  ) +
  theme_linedraw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  )

plot5

ggsave("pcMedicago/20250714PCEnrichClusterMarkers.pdf",
       plot5,
       height = 8, width = 10,
       bg = "transparent",
       dpi = 300)

# I want to try GSEA on marker genes
library(dplyr)
library(fgsea)

# Create a list of your cell-type marker gene sets
marker_sets <- split(marker_df$gene_id, marker_df$cluster)

# Prepare an empty list to hold results
gsea_results_list <- list()

# Loop through each PC
for (pc in unique(loading_long$PC)) {
  
  # Subset for this PC
  df_pc <- loading_long %>% filter(PC == pc)
  
  # Build named vector of loadings
  gene_ranks <- df_pc$Loading
  names(gene_ranks) <- df_pc$GENEID
  
  # Sort decreasing — fgsea expects highest ranks first if "enriched at top"
  gene_ranks <- sort(gene_ranks, decreasing = TRUE)
  
  # Run fgsea
  fgseaRes <- fgsea(pathways = marker_sets,
                    stats = gene_ranks,
                    minSize = 20,    # Adjust depending on marker set sizes
                    maxSize = 500,
                    nperm = 10000)
  
  # Add PC info
  fgseaRes$PC <- pc
  
  # Save
  gsea_results_list[[pc]] <- fgseaRes
}

# Combine all results
gsea_all <- bind_rows(gsea_results_list)

# Adjust p-values globally if you want (across all PCs and clusters)
gsea_all$padj_global <- p.adjust(gsea_all$pval, method = "BH")
gsea_sub = gsea_all[which(gsea_all$PC %in% c("PC1", "PC2",
                                             "PC3", "PC4",
                                             "PC5", "PC6",
                                             "PC7", "PC8",
                                             "PC9", "PC10",
                                             "PC11", "PC12",
                                             "PC13", "PC14",
                                             "PC15", "PC16",
                                             "PC17", "PC18"))]
gsea_sub$PC = factor(gsea_sub$PC, levels = c("PC1", "PC2",
                                             "PC3", "PC4",
                                             "PC5", "PC6",
                                             "PC7", "PC8",
                                             "PC9", "PC10",
                                             "PC11", "PC12",
                                             "PC13", "PC14",
                                             "PC15", "PC16",
                                             "PC17", "PC18"))
gsea_sub = left_join(gsea_sub, cluster_annot,
                     by = c("pathway" = "cluster"))
gsea_sig = gsea_sub[which(gsea_sub$padj_global <= 0.05),]

# Look at it!
plot1 = ggplot(gsea_sig, aes(x = PC, y = ClusterType, fill = NES)) +
  geom_tile() +
  scale_fill_gradientn(colors=c("#00441b","#1b7837","#5aae61","#a6dba0",
                                 "#d9f0d3","#f7f7f7","#e7d4e8","#c2a5cf",
                                 "#9970ab","#762a83","#40004b"),limits = c(-4,4),
                        oob = scales::squish) +
  theme_linedraw() +
  labs(x = "",
       y = "",
       fill = "NES") +
  theme(
    axis.title = element_text(size = 14,
                              face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12,
                                face = "bold"),
    legend.text = element_text(size = 10)
  )

ggsave("pcMedicago/20250714PCEnrichmentSCMarkerGenes.pdf",
       plot1,
       dpi = 300,
       bg = "transparent",
       width = 12, height = 6)
