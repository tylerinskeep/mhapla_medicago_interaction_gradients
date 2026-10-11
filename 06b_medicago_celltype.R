# 9 June 2026
library(Seurat)

seu = readRDS("00_source_data/GSE210881_sNucRNA-seq_Medicago_integrated_root_atlas.rds")
DimPlot(seu)

# I only want cells from mock treatment 
seu_mock <- subset(
  seu,
  subset = ppod1ocol %in% c("medicago_mock1_spx",
                            "medicago_mock2_spx",
                            "medicago_mock3_spx")
)
rm(seu)
DimPlot(seu_mock)

cluster_to_celltype <- c(
  "0" = "Cortex",
  "1" = "Cortex",
  "2" = "Stele",
  "3" = "Cortex",
  "4" = "Epidermis",
  "5" = "Endodermis",
  "6" = "Pericycle",
  "7" = "Cortex",
  "8" = "Epidermis",
  "9" = "Cortex",
  "10" = "Cortex",
  "11" = "Stem Cell Niche",
  "12" = "Cortex",
  "13" = "Endodermis",
  "14" = "Endodermis",
  "15" = "Phloem",
  "16" = "Unassigned",
  "17" = "Epidermis",
  "18" = "Epidermis",
  "19" = "Xylem",
  "20" = "Epidermis",
  "21" = "Endodermis",
  "22" = "Xylem",
  "23" = "Xylem",
  "24" = "Endodermis"
)

seu_mock$cell_type <- unname(cluster_to_celltype[as.character(seu_mock$seurat_clusters)])

DimPlot(seu_mock, 
        group.by = "cell_type")

library(Seurat)
library(dplyr)
library(pheatmap)

# Make sure this exists
table(seu_mock$cell_type)

# Set identity to cell_type
Idents(seu_mock) <- "cell_type"

# Average normalized expression by cell type
avg_expr <- AverageExpression(
  seu_mock,
  assays = "RNA",
  slot = "data",
  group.by = "cell_type"
)$RNA

# Identifying cell type specific genes
celltype_markers <- FindAllMarkers(
  seu_mock,
  assay = "RNA",
  slot = "data",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

top_markers <- celltype_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE) %>%
  ungroup()

DoHeatmap(
  seu_mock,
  features = unique(top_markers$gene),
  group.by = "cell_type",
  size = 3
) +
  NoLegend()

library(pheatmap)

avg_expr <- AverageExpression(
  seu_mock,
  assays = "RNA",
  slot = "data",
  group.by = "cell_type"
)$RNA

genes_to_plot <- unique(celltype_markers$gene)
genes_to_plot <- intersect(genes_to_plot, rownames(avg_expr))

mat <- avg_expr[genes_to_plot, , drop = FALSE]

# Z-score each gene across cell types
mat_z <- t(scale(t(mat)))

# Optional: order genes by marker cell type
gene_anno <- top_markers %>%
  select(gene, cluster) %>%
  distinct()

gene_anno <- as.data.frame(gene_anno)
rownames(gene_anno) <- gene_anno$gene
gene_anno$gene <- NULL
colnames(gene_anno) <- "marker_cell_type"

gene_anno <- gene_anno[rownames(mat_z), , drop = FALSE]

pheatmap(
  mat_z,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  annotation_row = gene_anno,
  fontsize_row = 6,
  fontsize_col = 10,
  main = "Cell type marker expression"
)

features <- read.delim(
  "00_source_data/pereira_snrna/features.tsv",
  header = FALSE,
  stringsAsFactors = FALSE
)

head(features)

feature_lookup <- features %>%
  transmute(
    gene_id = V1,
    feature_name = V2
  )

sum(rownames(seu_mock) %in% feature_lookup$feature_name)
mean(rownames(seu_mock) %in% feature_lookup$feature_name)

DefaultAssay(seu_mock) <- "RNA"

nrow(seu_mock)
length(rownames(seu_mock))
head(rownames(seu_mock), 50)

library(Seurat)
library(dplyr)
library(ggplot2)
pca_ids_clean <- gsub("MtrunA17_Chr", "MtrunA17Chr", loadings_mapped$GENEID)

sum(pca_ids_clean %in% rownames(seu_mock))
mean(pca_ids_clean %in% rownames(seu_mock))

library(dplyr)
library(Seurat)

DefaultAssay(seu_mock) <- "RNA"

loadings_mapped$GENEID_seurat <- gsub("MtrunA17_Chr", "MtrunA17Chr", loadings_mapped$GENEID)

loadings_overlap <- loadings_mapped %>%
  filter(GENEID_seurat %in% rownames(seu_mock))

nrow(loadings_overlap)


n_top <- 100

pc_gene_sets <- list(
  PC3_positive = loadings_mapped %>%
    arrange(desc(PC3)) %>%
    slice_head(n = n_top) %>%
    pull(GENEID),
  
  PC3_negative = loadings_mapped %>%
    arrange(PC3) %>%
    slice_head(n = n_top) %>%
    pull(GENEID),
  
  PC4_positive = loadings_mapped %>%
    arrange(desc(PC4)) %>%
    slice_head(n = n_top) %>%
    pull(GENEID),
  
  PC4_negative = loadings_mapped %>%
    arrange(PC4) %>%
    slice_head(n = n_top) %>%
    pull(GENEID)
)

pc_gene_sets <- lapply(pc_gene_sets, function(x) {
  gsub("MtrunA17_Chr", "MtrunA17Chr", x)
})

sapply(pc_gene_sets, function(x) sum(x %in% rownames(seu_mock)))

seu_mock <- AddModuleScore(
  object = seu_mock,
  features = pc_gene_sets,
  assay = "RNA",
  name = "PC_loading"
)

score_cols <- grep("^PC_loading", colnames(seu_mock@meta.data), value = TRUE)
score_cols <- tail(score_cols, length(pc_gene_sets))

for (i in seq_along(score_cols)) {
  colnames(seu_mock@meta.data)[
    colnames(seu_mock@meta.data) == score_cols[i]
  ] <- names(pc_gene_sets)[i]
}

VlnPlot(
  seu_mock,
  features = c("PC3_positive", "PC3_negative", "PC4_positive", "PC4_negative"),
  group.by = "cell_type",
  pt.size = 0,
  ncol = 2
)

# How many total PCA genes can be found in Seurat?
pca_ids <- loadings_mapped$GENEID
pca_ids_clean <- gsub("MtrunA17_Chr", "MtrunA17Chr", pca_ids)

sum(pca_ids_clean %in% rownames(seu_mock))
mean(pca_ids_clean %in% rownames(seu_mock))



VlnPlot(
  seu_mock,
  features = c("PC3_positive", "PC3_negative", "PC4_positive", "PC4_negative"),
  group.by = "cell_type",
  pt.size = 0,
  ncol = 2
)



# So many common names among the locus IDs... Trying the Pereira Data instead...
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)

data_dir <- "00_source_data/pereira_snrna/"

counts <- Read10X(
  data.dir = "00_source_data/pereira_snrna/",
  gene.column = 1
)

seu_pereira <- CreateSeuratObject(
  counts = counts,
  project = "Pereira_Medicago_root_0h",
  min.cells = 3,
  min.features = 200
)

seu_pereira <- NormalizeData(seu_pereira)
seu_pereira <- FindVariableFeatures(seu_pereira, nfeatures = 3000)
seu_pereira <- ScaleData(seu_pereira)
seu_pereira <- RunPCA(seu_pereira)
ElbowPlot(seu_pereira)
seu_pereira <- FindNeighbors(seu_pereira, dims = 1:20)
seu_pereira <- FindClusters(seu_pereira, resolution = 0.5)
seu_pereira <- RunUMAP(seu_pereira, dims = 1:30)

DimPlot(seu_pereira, label = TRUE)


library(dplyr)
library(tidyr)
library(ggplot2)
meta <- seu_mock@meta.data

# Keep last occurrence of each column name
meta_clean <- meta[, !duplicated(colnames(meta), fromLast = TRUE)]

seu_mock@meta.data <- meta_clean

library(dplyr)
library(tidyr)
library(ggplot2)

module_long <- seu_mock@meta.data %>%
  select(cell_type, PC3_positive, PC3_negative, PC4_positive, PC4_negative) %>%
  pivot_longer(
    cols = c(PC3_positive, PC3_negative, PC4_positive, PC4_negative),
    names_to = "module",
    values_to = "score"
  ) %>%
  separate(module, into = c("PC", "loading_side"), sep = "_") %>%
  mutate(
    loading_side = factor(loading_side, levels = c("negative", "positive")),
    PC = factor(PC, levels = c("PC3", "PC4"))
  )


pc3_pos = module_long[which(module_long$PC == "PC3" & module_long$loading_side == "positive"),]
pc3_pos = pc3_pos[!pc3_pos$cell_type %in% c("Unassigned"),]
ggplot(pc3_pos,
       aes(x = cell_type, 
           y = score, 
           fill = cell_type)) +
  geom_violin(
    scale = "width",
    trim = TRUE,
    position = position_dodge(width = 0.85),
    linewidth = 0.25
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    position = position_dodge(width = 0.85),
    linewidth = 0.25
  ) +
  scale_fill_manual(
    values = wes_palette("Zissou1", 8, "continuous")
  ) +
  labs(
    x = NULL,
    y = "upper loadings",
    fill = NULL
  ) +
  theme_tyler() +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
ggsave("06_PCA_medicago/pc3_upper_100_modulescore_cervantes-perez.pdf",
       height = 3, width = 6)

#### Cell-type marker enrichment across PC loading gene sets ####

library(dplyr)
library(tidyr)
library(ggplot2)

outdir <- "06_PCA_medicago"
### -----------------------------------------
### 1. Make sure you have cell-type markers
### -----------------------------------------

DefaultAssay(seu_mock) <- "RNA"
Idents(seu_mock) <- "cell_type"

# Run this only if you have not already calculated markers
celltype_markers <- FindAllMarkers(
  seu_mock,
  assay = "RNA",
  slot = "data",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.csv(
  celltype_markers,
  file.path(outdir, "mock_cell_type_markers.csv"),
  row.names = FALSE
)

### -----------------------------------------
### 2. Define cell-type-specific marker genes
### -----------------------------------------

marker_fdr_cutoff <- 0.05
marker_logfc_cutoff <- 0.5
top_n_markers <- 300

celltype_marker_sets <- celltype_markers %>%
  filter(p_val_adj < marker_fdr_cutoff) %>%
  filter(avg_log2FC > marker_logfc_cutoff) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = top_n_markers) %>%
  ungroup() %>%
  select(cell_type = cluster, gene) %>%
  distinct()

marker_sets_list <- split(celltype_marker_sets$gene, celltype_marker_sets$cell_type)

sapply(marker_sets_list, length)

### -----------------------------------------
### 3. Define PC loading gene sets
### -----------------------------------------

n_top_loading <- 300

loadings_mapped <- loadings_mapped %>%
  mutate(
    GENEID_seurat = gsub("MtrunA17_Chr", "MtrunA17Chr", GENEID)
  )

# Use only genes present in the Seurat RNA assay
loadings_overlap <- loadings_mapped %>%
  filter(GENEID_seurat %in% rownames(seu_mock))

pc_loading_sets <- list(
  PC1_positive = loadings_overlap %>%
    arrange(desc(PC1)) %>%
    slice_head(n = n_top_loading) %>%
    pull(GENEID_seurat),
  
  PC1_negative = loadings_overlap %>%
    arrange(PC1) %>%
    slice_head(n = n_top_loading) %>%
    pull(GENEID_seurat),
  
  PC2_positive = loadings_overlap %>%
    arrange(desc(PC2)) %>%
    slice_head(n = n_top_loading) %>%
    pull(GENEID_seurat),
  
  PC2_negative = loadings_overlap %>%
    arrange(PC2) %>%
    slice_head(n = n_top_loading) %>%
    pull(GENEID_seurat),
  
  PC3_positive = loadings_overlap %>%
    arrange(desc(PC3)) %>%
    slice_head(n = n_top_loading) %>%
    pull(GENEID_seurat),
  
  PC3_negative = loadings_overlap %>%
    arrange(PC3) %>%
    slice_head(n = n_top_loading) %>%
    pull(GENEID_seurat),
  
  PC4_positive = loadings_overlap %>%
    arrange(desc(PC4)) %>%
    slice_head(n = n_top_loading) %>%
    pull(GENEID_seurat),
  
  PC4_negative = loadings_overlap %>%
    arrange(PC4) %>%
    slice_head(n = n_top_loading) %>%
    pull(GENEID_seurat)
)

pcs_to_test <- paste0("PC", 1:11)

pc_loading_sets <- list()

for (pc in pcs_to_test) {
  pc_loading_sets[[paste0(pc, "_positive")]] <- loadings_overlap %>%
    arrange(desc(.data[[pc]])) %>%
    slice_head(n = n_top_loading) %>%
    pull(GENEID_seurat) %>%
    unique()
  
  pc_loading_sets[[paste0(pc, "_negative")]] <- loadings_overlap %>%
    arrange(.data[[pc]]) %>%
    slice_head(n = n_top_loading) %>%
    pull(GENEID_seurat) %>%
    unique()
}

sapply(pc_loading_sets, length)

### -----------------------------------------
### 4. Fisher enrichment test
### -----------------------------------------

# The background should be genes available to both:
# 1. PCA loadings
# 2. Seurat RNA assay
# 3. marker calling universe
background_genes <- unique(loadings_overlap$GENEID_seurat)

# Optional stricter background: genes in Seurat assay and PCA overlap
background_genes <- intersect(background_genes, rownames(seu_mock))

fisher_enrich <- function(query_genes, marker_genes, background_genes) {
  query_genes <- intersect(unique(query_genes), background_genes)
  marker_genes <- intersect(unique(marker_genes), background_genes)
  
  a <- length(intersect(query_genes, marker_genes))
  b <- length(setdiff(query_genes, marker_genes))
  c <- length(setdiff(marker_genes, query_genes))
  d <- length(setdiff(background_genes, union(query_genes, marker_genes)))
  
  mat <- matrix(c(a, b, c, d), nrow = 2)
  
  test <- fisher.test(mat, alternative = "two.sided")
  
  data.frame(
    overlap = a,
    query_size = length(query_genes),
    marker_size = length(marker_genes),
    background_size = length(background_genes),
    odds_ratio = unname(test$estimate),
    log2_odds_ratio = log2(unname(test$estimate)),
    p_value = test$p.value
  )
}

enrich_results <- bind_rows(
  lapply(names(pc_loading_sets), function(pc_set_name) {
    bind_rows(
      lapply(names(marker_sets_list), function(celltype_name) {
        res <- fisher_enrich(
          query_genes = pc_loading_sets[[pc_set_name]],
          marker_genes = marker_sets_list[[celltype_name]],
          background_genes = background_genes
        )
        
        res$pc_set <- pc_set_name
        res$cell_type <- celltype_name
        res
      })
    )
  })
)

enrich_results <- enrich_results %>%
  separate(pc_set, into = c("PC", "loading_side"), sep = "_", remove = FALSE) %>%
  group_by(PC, loading_side) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    neglog10_fdr = -log10(p_adj),
    enriched = p_adj < 0.05 & odds_ratio > 1
  ) %>%
  ungroup()

write.csv(
  enrich_results,
  file.path(outdir, "PC_loading_cell_type_marker_enrichment_fisher.csv"),
  row.names = FALSE
)

### -----------------------------------------
### 5. Plot enrichment as log2 odds ratios
### -----------------------------------------

plot_df <- enrich_results %>%
  filter(p_adj < 0.10 | abs(log2_odds_ratio) > 0.5) %>%
  mutate(
    label = paste(PC, loading_side, sep = " "),
    cell_type = factor(cell_type)
  )

ggplot(plot_df, aes(x = log2_odds_ratio, y = cell_type, fill = neglog10_fdr)) +
  geom_vline(xintercept = 0, linewidth = 0.7) +
  geom_col(width = 0.75) +
  facet_grid(PC ~ loading_side, scales = "free_y") +
  scale_fill_gradient(
    low = "#D8D8D8",
    high = "#214F32",
    name = expression(-log[10]("FDR"))
  ) +
  labs(
    x = expression(log[2]("odds ratio")),
    y = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "right"
  )

library(dplyr)
library(ggplot2)
library(stringr)

# choose PCs to show
pcs_show <- paste0("PC", 1:11)

heat_df <- enrich_results %>%
  filter(PC %in% pcs_show) %>%
  mutate(
    module = paste0(PC, "_", loading_side),
    module = factor(
      module,
      levels = as.vector(rbind(
        paste0(pcs_show, "_positive"),
        paste0(pcs_show, "_negative")
      ))
    ),
    sig_label = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

heat_df$PC = gsub("PC", "", heat_df$PC)
heat_df$PC = factor(heat_df$PC,
                    levels = 1:11)
heat_df = heat_df[!heat_df$cell_type %in% c("Unassigned"),]

p_celltype_enrich_heatmap <- ggplot(
  heat_df,
  aes(x = cell_type, y = PC, fill = log2_odds_ratio)
) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(
    aes(label = sig_label),
    color = "black",
    size = 4,
    fontface = "bold"
  ) +
  scale_fill_gradientn(
    colors = c(
      "#40004b", "#762a83", "#9970ab", "#c2a5cf",
      "#e7d4e8", "#f7f7f7", "#d9f0d3", "#a6dba0",
      "#5aae61", "#1b7837", "#00441b"
    ),
    limits = c(-4, 4),
    name = expression(log[2]("odds ratio"))
    
  ) +
  labs(
    x = NULL,
    y = "Principal component"
  ) +
  facet_grid(~loading_side) +
  theme_tyler() +
  theme(
    axis.text.x = element_text(angle = 45,
                               vjust = 0.6)
  )

p_celltype_enrich_heatmap

ggsave("06_PCA_medicago/PC_enrich_cell_type_markers.pdf",
       p_celltype_enrich_heatmap,
       height = 6, width = 6)
