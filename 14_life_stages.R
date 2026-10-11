#### Predicting Life Stage-Specific Expression for M hapla ####
## Tyler Inkseep, Groen Lab, UC Riverside ##
# 30 June 2025 #

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
orthology = read.delim("00_source_data/hapla_incognita_orthologs.tsv")
ortholog_long <- orthology %>%
  separate_rows(m_hapla, sep = ",\\s*") %>%   # Split Mh genes into rows
  separate_rows(m_incognita, sep = ",\\s*")     # Split Minc genes into rows

lifeStage = read.csv("00_source_data/m_incognita_lifeStageTau.csv")

expression = read.delim("14_life_stages/choi_incognita_molts_tpm.tsv")

# Plotting the expression of life stage-specific genes
lss = lifeStage[which(lifeStage$tau >= 0.7),]

lss_expr = expression[which(expression$gene_id %in% lss$GENEID),]

heatmap_mat <- lss_expr[,2:6 ]

# Z-score by gene, so each row shows relative expression pattern
heatmap_z <- t(scale(t(heatmap_mat)))
heatmap_z[is.na(heatmap_z)] <- 0

# One continuous color gradient
heat_cols <- colorRampPalette(c("#00441b","#1b7837","#5aae61","#a6dba0",
                                "#d9f0d3","#f7f7f7","#e7d4e8","#c2a5cf",
                                "#9970ab","#762a83","#40004b"))(101)

plot1 = pheatmap(
  heatmap_z,
  color = heat_cols,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  fontsize_col = 9,
  border_color = NA,
  main = "Life-Stage Expression M. incognita"
)

ggsave("14_life_stages/m_incognita_life_stage_expression_heatmap.pdf",
       plot1,
       height = 4, width = 4)

stage_filter_counts <- tibble(
  
  step = factor(
    
    c("Total\nM. hapla genes",
      
      "Expressed\nin dataset",
      
      "Transferred\nstage labels"),
    
    levels = c("Total\nM. hapla genes",
               
               "Expressed\nin dataset",
               
               "Transferred\nstage labels")
    
  ),
  
  n = c(10754, 8686, 6149)
  
)

ggplot(stage_filter_counts, aes(x = step, y = n)) +
  geom_col(width = 0.65, fill = "#533943") +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    x = NULL,
    y = "Number of genes"
  ) +
  theme_tyler() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1))

ggsave("14_life_stages/mhapla_life_stage_annotation_stats.pdf",
       height = 4, width = 4)

# ============================================================
# Add sample/run IDs from cell barcodes

# Merge the datasets
dt = ortholog_long %>%
  left_join(lifeStage %>%
              select(GENEID, tau, preferred_stage),
            by = c("m_incognita" = "GENEID"))
dt = dt %>%
  mutate(preferred_stage = ifelse(is.na(tau),
                                  NA,
                                  preferred_stage))
write.csv(dt,
          "14_life_stages/hapla_incognita_lifeStageTau.csv")

# I first wanna summarize the variablity in stage specificity across orthologs
mh_variability <- dt %>%
  group_by(m_hapla) %>%
  summarize(
    NumMincParalogs = n_distinct(m_incognita),
    NumStageCategories = n_distinct(preferred_stage[!is.na(preferred_stage)]),
    StageCategories = paste(sort(unique(preferred_stage[!is.na(preferred_stage)])), collapse = ";"),
    ConsensusStatus = case_when(
      NumStageCategories == 1 ~ "Conserved",
      NumStageCategories > 1 ~ "Divergent",
      TRUE ~ "Unknown"
    )
  )

write.csv(mh_variability,
          "14_life_stages/hapla_ortholog_lifeStageVariability.csv")

plot1 = ggplot(mh_variability,
               aes(NumMincParalogs)) +
  geom_histogram()

plot1

total_mh_genes <- nrow(mh_variability) #10,754 genes in M hapla
num_multi_mh_genes <- sum(mh_variability$NumMincParalogs > 1,
                          na.rm = TRUE) # 7,524 of those genes have multiple orthologs in M. incognita
7542/10754 # that's 70.13% of all genes in the genome

multi_mh_genes <- mh_variability %>%
  filter(NumMincParalogs > 1)

num_divergent_in_multi <- sum(multi_mh_genes$ConsensusStatus == "Divergent",
                              na.rm = TRUE) # of genes w/ multiple orthologs, 1,280 have divergent expression profiles
1280/7542 # 16.97% of such genes exhibit paralogs with divergent expression profiles

# That's 11.9% of all genes in the genome - that's not terrible!

# Now, I want to know where these genes are coming from - are there life stages that
# seem to be more variable?
dt <- dt %>%
  left_join(mh_variability, by = "m_hapla")

# Which genes have orthologs with divergent expression patterns
divergent_orthologs = dt %>%
  filter(ConsensusStatus == "Divergent",
         !is.na(preferred_stage))

conserved_orthologs = dt %>%
  filter(ConsensusStatus == "Conserved",
         !is.na(preferred_stage))

dt2 = bind_rows(divergent_orthologs %>% mutate(Group = "Divergent"),
                conserved_orthologs %>% mutate(Group = "Conserved")
)

dt2$preferred_stage = gsub("adult.female",
                           "Adult female",
                           dt2$preferred_stage)
dt2$preferred_stage = gsub("J2.juveniles",
                           "J2",
                           dt2$preferred_stage)
dt2$preferred_stage = gsub("Constitutive",
                           "Broad",
                           dt2$preferred_stage)
dt2$preferred_stage = factor(dt2$preferred_stage,
                             levels = c("Broad",
                                        "Egg",
                                        "J2",
                                        "J3",
                                        "J4",
                                        "Adult female"))

plot2 = ggplot(dt2, aes(x = preferred_stage,
                fill = Group)) +
  geom_bar(position = "dodge") +
  labs(x = "Life Stage", y = "Number of Paralogs",
       fill = "Expression \nConservation \nAmong Paralogs") +
  scale_fill_manual(values = wes_palette("Zissou1", 2, "continuous")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14,
                                  face = "bold"),
        legend.title = element_text(size = 14,
                                    face = "bold"),
        legend.text = element_text(size = 12))
plot2

ggsave("14_life_stages/incognita_genes_expression_among_paralogs.pdf",
       plot2,
       height = 4, width = 4,
       bg = "transparent",
       dpi = 600)

# Test for enrichment/depletion
# Contingency table: rows = Group, columns = Category
table_minc <- table(dt2$Group, dt2$preferred_stage)

# Print table
print(table_minc)

#.        Constitutive   Egg    J2    J3    J4      Adult female
#Conserved        17099   644  1152   191   301          414
#Divergent         2543   673   703   283   305          333

# Test for enrichment/depletion
chisq_test <- chisq.test(table_minc)
print(chisq_test)
# X-squared = 2913.1, df = 5, p-value < 2.2e-16
# Using the standardized residuals to infer enrichment; plotting
std_residuals <- chisq_test$stdres
residuals_df <- as.data.frame(as.table(std_residuals))
colnames(residuals_df) <- c("Group", "Category", "StdResidual")

plot3 = ggplot(residuals_df, aes(x = Category, y = Group, fill = StdResidual)) +
  geom_tile() +
  scale_fill_gradientn(colors=c("#00441b","#1b7837","#5aae61","#a6dba0",
                                "#d9f0d3","#f7f7f7","#e7d4e8","#c2a5cf",
                                "#9970ab","#762a83","#40004b"),limits = c(-50,50),
                       oob = scales::squish) +
  labs(title = "Chi-Square Standardized Residuals",
       x = "Life Stage", y = "Expression Conservation Among Paralogs") +
  theme_classic() +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14,
                                    face = "bold"),
          legend.title = element_text(size = 14,
                                      face = "bold"),
          legend.text = element_text(size = 12))

plot3

ggsave("lifeStage/20250630ExpressionConservationStdResidual.pdf",
       plot3,
       height = 6, width = 8,
       bg = "transparent",
       dpi = 300)

# What impact does gene family size have on expression conservation?
plot4 = ggplot(mh_variability %>% filter(ConsensusStatus %in% c("Conserved",
                                                                "Divergent")),
       aes(x = ConsensusStatus,
           y = NumMincParalogs,
           fill = ConsensusStatus)) +
  geom_boxplot() +
  labs(x = "Expression Conservation Among M. incognita Paralogs",
       y = "Number of M. incognita orthologs per M. hapla gene") +
  scale_fill_manual(values = wes_palette("Zissou1", 2, "continuous")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14,
                                  face = "bold"),
        legend.title = element_blank(),
        legend.text = element_blank())

plot4

ggsave("lifeStage/20250630ExpressionConservationByOrthogroupSize.pdf",
       plot4,
       height = 6, width = 6,
       bg = "transparent",
       dpi = 300)

# Filter out unknown cases
filtered_summary <- mh_variability %>%
  filter(ConsensusStatus %in% c("Conserved", "Divergent"))

# Wilcoxon rank-sum test (non-parametric)
wilcox.test(NumMincParalogs ~ ConsensusStatus, data = filtered_summary)
#Wilcoxon rank sum test with continuity correction
#data:  NumMincParalogs by ConsensusStatus
#W = 2652694, p-value < 2.2e-16
#alternative hypothesis: true location shift is not equal to 0


# I want to measure the similarity of tau across orthogroups as well - Can I transfer these?
tau_similarity <- dt %>%
  group_by(m_hapla) %>%
  filter(n() > 1, !is.na(tau)) %>%
  summarize(
    NumMincParalogs = n(),
    Tau_SD = sd(tau),
    Tau_Range = max(tau) - min(tau),
    Tau_Mean = mean(tau),
    Tau_Min = min(tau),
    Tau_Max = max(tau),
    NumStageCategories = n_distinct(preferred_stage[!is.na(preferred_stage)]),
    ExpressionDivergence = case_when(
      NumStageCategories == 1 ~ "Conserved",
      NumStageCategories > 1 ~ "Divergent",
      TRUE ~ "Unknown"
    )
  )

plot5 = ggplot(tau_similarity, aes(x = ExpressionDivergence,
                           y = Tau_SD,
                           fill = ExpressionDivergence)) +
  geom_boxplot() +
  labs(x = "Expression Conservation Among M. Incognita Paralogs",
       y = "Stdev (Orthogroup LSS)") +
  scale_fill_manual(values = wes_palette("Zissou1", 2, "continuous")) +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14,
                                  face = "bold"),
        legend.title = element_blank(),
        legend.text = element_blank())

plot5

ggsave("lifeStage/20250630ExpressionConservationByOrthogroupTauSTDEV.pdf",
       plot5,
       height = 6, width = 6,
       bg = "transparent",
       dpi = 300)

# Filter out unknown cases
filtered_summary <- tau_similarity %>%
  filter(ExpressionDivergence %in% c("Conserved", "Divergent"))

# Wilcoxon rank-sum test (non-parametric)
wilcox.test(Tau_SD ~ ExpressionDivergence, data = filtered_summary)
#Wilcoxon rank sum test with continuity correction
#data:  NumMincParalogs by ExpressionDivergence
# W = 1837696, p-value < 2.2e-16

plot6 = ggplot(tau_similarity, aes(x = ExpressionDivergence,
                                   y = Tau_Mean,
                                   fill = ExpressionDivergence)) +
  geom_boxplot() +
  labs(x = "Expression Conservation Among M. Incognita Paralogs",
       y = "Mean Orthogroup LSS") +
  scale_fill_manual(values = wes_palette("Zissou1", 2, "continuous")) +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14,
                                  face = "bold"),
        legend.title = element_blank(),
        legend.text = element_blank())

plot6

ggsave("lifeStage/20250630ExpressionConservationByOrthogroupTauMEAN.pdf",
       plot6,
       height = 6, width = 6,
       bg = "transparent",
       dpi = 300)

# Wilcoxon rank-sum test (non-parametric)
wilcox.test(Tau_Mean ~ ExpressionDivergence, data = filtered_summary)
#Wilcoxon rank sum test with continuity correction
#data:  Tau_Mean by ExpressionDivergence
# W = 1286984, p-value < 2.2e-16


### IMPORTING INTERACTION MODEL ###
model = read.csv("selectionModel/hapla_selectionModel.csv")
model$geneid = gsub("_mRNA", "",
                    model$geneid)

# Merging the data with the interaction model
dt = left_join(tau_similarity,
               model,
               by = c("m_hapla" = "geneid"))
dt_conserved = dt[which(dt$ExpressionDivergence == "Conserved"),]               

# Removing transcripts that weren't detected
dt_conserved = dt_conserved %>%
  filter(!is.na(S)) # 3,533 genes retained (40% of all the genes)

# Is this small sample similar in distribution to the entirety?
wilcox.test(model$S, dt_conserved$S)
#Wilcoxon rank sum test with continuity correction
#data:  model$S and dt_conserved$S
#W = 15312560, p-value = 0.8675

wilcox.test(model$C, dt_conserved$C)
#Wilcoxon rank sum test with continuity correction
#data:  model$C and dt_conserved$C
#W = 15413876, p-value = 0.6845

# Binning the life stage specificity scores
# Step 1: Bin the tau values
dt_conserved$tau_bins <- cut(dt_conserved$Tau_Mean, breaks=seq(0, 1, by=0.1), 
                   labels=paste0(seq(0, 0.9, by=0.1),
                                 "-", seq(0.1, 1, by=0.1)),
                   include.lowest=TRUE)

dt_conserved = dt_conserved[complete.cases(dt_conserved$tau_bins),]

palette <- wes_palette("Zissou1", 10,
                       type = "continuous")

plot7 = ggplot(dt_conserved, aes(x=tau_bins, 
                       y=S,
                       fill = tau_bins)) +
  geom_boxplot(alpha = 0.4,
               show.legend = FALSE) +
  scale_fill_manual(values = palette) +
  theme_linedraw() +
  labs(x = "Life Stage Specificity",
       y = "S") +
  theme(
    axis.text = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),       # Size of axis titles
    axis.title.x = element_text(size = 14, face = "bold"), # Custom style for x-axis title
    axis.title.y = element_text(size = 14, face = "bold")
  ) 
plot7

ggsave("lifeStage/20250630lifeStageTauS.pdf",
       plot7,
       dpi = 300,
       width = 8, height = 6,
       bg = "transparent")

kruskal = kruskal.test(S ~ tau_bins,
                       data = dt_conserved)
print(kruskal)
#Kruskal-Wallis rank sum test
#data:  S by tau_bins
#Kruskal-Wallis chi-squared = 12.174, df = 9, p-value = 0.2037

dt_conserved$absS = abs(dt_conserved$S)

kruskal = kruskal.test(absS ~ tau_bins,
                       data = dt_conserved)
print(kruskal)

#Kruskal-Wallis rank sum test
#data:  absS by tau_bins
#Kruskal-Wallis chi-squared = 241.49, df = 9, p-value < 2.2e-16

# Dunn's Post-Hoc Test
test = FSA::dunnTest(absS ~ tau_bins,
                     data = dt_conserved,
                     method = "bonferroni")
test
# 0-0.1 = AB
# 0.9-1 = B
# ALL ELSE = A

plot8 = ggplot(dt_conserved, aes(x=tau_bins, 
                       y=absS,
                       fill = tau_bins)) +
  geom_boxplot(alpha = 0.4,
               show.legend = FALSE) +
  scale_fill_manual(values = palette) +
  theme_linedraw() +
  labs(x = "Life Stage Specificity (in M. incognita)",
       y = "|S|") +
  theme(
    axis.text = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),       # Size of axis titles
    axis.title.x = element_text(size = 14, face = "bold"), # Custom style for x-axis title
    axis.title.y = element_text(size = 14, face = "bold")
  ) 
plot8

ggsave("lifeStage/20250630lifeStageabsS.pdf",
       plot8,
       dpi = 300,
       width = 8, height = 6,
       bg = "transparent")

kruskal = kruskal.test(C ~ tau_bins,
                       data = dt_conserved)
print(kruskal)
#Kruskal-Wallis rank sum test
#data:  C by tau_bins
#Kruskal-Wallis chi-squared = 145.64, df = 9, p-value < 2.2e-16
# Dunn's Post-Hoc Test
test = FSA::dunnTest(C ~ tau_bins,
                     data = dt_conserved,
                     method = "bonferroni")
test

dt_summary <- dt_conserved %>%
  group_by(tau_bins) %>%
  summarize(
    mean_C = mean(C, na.rm = TRUE),
    se_C = sd(C, na.rm = TRUE) / sqrt(n())
  )

plot9 = ggplot(dt_summary, aes(x = tau_bins, y = mean_C, fill = tau_bins)) +
  geom_col(alpha = 0.7, show.legend = FALSE) +
  geom_errorbar(aes(ymin = mean_C - se_C, ymax = mean_C + se_C), width = 0.2) +
  scale_fill_manual(values = palette) +
  theme_linedraw() +
  labs(x = "Life Stage Specificity (in M. incognita)",
       y = "C") +
  theme(
    axis.text = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold")
  )

plot9

ggsave("lifeStage/20250630lifeStageC.pdf",
       plot9,
       dpi = 300,
       width = 8, height = 6,
       bg = "transparent")


# Now, I want to see if the association is driven by any given life stage
dt_conserved = left_join(dt_conserved,
                         conserved_orthologs,
                         by = "m_hapla")
dt_conserved$StageCategories = gsub("adult.female",
                                    "Adult female",
                                    dt_conserved$StageCategories)
dt_conserved$StageCategories = gsub("J2.juveniles",
                                    "J2",
                                    dt_conserved$StageCategories)
dt_conserved$StageCategories = gsub("Constitutive",
                                    "Broad",
                                    dt_conserved$StageCategories)
dt_conserved$StageCategories = factor(dt_conserved$StageCategories,
                                      levels = c("Broad",
                                                 "Egg",
                                                 "J2",
                                                 "J3",
                                                 "J4",
                                                 "Adult female"))

kruskal = kruskal.test(S ~ StageCategories,
                       data = dt_conserved)
print(kruskal)
#Kruskal-Wallis rank sum test
#data:  S by StageCategories
#Kruskal-Wallis chi-squared = 489.08, df = 5, p-value < 2.2e-16
# Dunn's Post-Hoc Test
test = FSA::dunnTest(S ~ StageCategories,
                     data = dt_conserved,
                     method = "bonferroni")
test


plot10 = ggplot(dt_conserved,
                aes(x = S,
                    y = StageCategories,
                    fill = StageCategories)) +
  geom_boxplot(alpha = 0.4,
               show.legend = FALSE) +
  scale_fill_manual(values = wes_palette("Zissou1", 6, "continuous")) +
  theme_linedraw() +
  labs(x = "S",
       y = "Stage-Specific Genes (M. incognita)") +
  theme(
    axis.text = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),       # Size of axis titles
    axis.title.x = element_text(size = 14, face = "bold"), # Custom style for x-axis title
    axis.title.y = element_text(size = 14, face = "bold")
  ) 
plot10

ggsave("lifeStage/20250630stageSpecificS.pdf",
       plot10,
       dpi = 300,
       width = 7, height = 6,
       bg = "transparent")


## NOW... IM TRYING TO SEE IF I CAN CLUSTER SAMPLES ON THE BASIS OF STAGE-SPECIFIC
# MARKER GENE EXPRESSION
# Get the stage-specific genes
# Filter to stage-specific only
stage_markers_df <- dt_conserved %>%
  filter(ExpressionDivergence != "Broad")  # Keep only stage-specific

# Vector of marker gene IDs
stage_marker_genes <- stage_markers_df$m_hapla # 643 stage-specific marker genes

# Import normalized counts
counts = read.csv("featureCounts/parsed_hapla_log2normCPM.csv")
counts$geneid = gsub("_mRNA", "", counts$X)

# Subsetting the count matrix to keep only the stage marker genes
stage_counts <- counts[counts$geneid %in% stage_markers, ] # grabbed 218 genes
geneids = stage_counts$geneid

# Scale expression per gene
expr_scaled <- t(scale(t(stage_counts[,2:84])))
row.names(expr_scaled) = geneids

# Annotating stage-specific genes
gene_ann <- stage_markers_df %>%
  filter(m_hapla %in% rownames(expr_scaled)) %>%
  distinct(m_hapla, .keep_all = TRUE) %>%   # <-- Remove duplicates, keep first
  select(m_hapla, StageCategories) %>%
  tibble::column_to_rownames("m_hapla")


plot11 = pheatmap(expr_scaled,
         annotation_row = gene_ann,
         # annotation_col = sample_metadata, # optional if you have
         cluster_cols = TRUE,        # Cluster samples
         cluster_rows = TRUE,        # Cluster genes
         show_rownames = FALSE,      # Hide gene names if too many
         show_colnames = FALSE,
         scale = "none",             # Already scaled
         fontsize_row = 6)
plot11

ggsave("lifeStage/20250630clusterByStageMarkers.pdf",
       plot11,
       dpi = 300,
       width = 8, height = 8,
       bg = "transparent")

# LET'S JUST TRY PRINCIPLA COMPONENT ANALYSIS
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

# Getting loadings for genes into each PC
loadings_df <- as.data.frame(pca_result$rotation)
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


# Merging the loadings with stage-specific markers
stage_markers_df_clean <- conserved_orthologs[!duplicated(conserved_orthologs$m_hapla),]

loading_long$geneid = gsub("_mRNA", "", loading_long$geneid)

loading_long <- loading_long %>%
  left_join(stage_markers_df_clean, by = c("geneid" = "m_hapla")) 

# Get unique PCs and marker stages
pcs <- unique(loading_long$PC)
stages <- unique(loading_long$StageCategories)
stages <- stages[!is.na(stages)]  # remove NA

# Initialize list to hold results
results_list <- list()

# Loop through PCs and stage categories
for (pc in pcs) {
  for (stage in stages) {
    # Subset to this PC
    df_pc <- loading_long[loading_long$PC == pc, ]
    
    # Marker genes for this stage
    marker_genes <- df_pc$StageCategories == stage & !is.na(df_pc$StageCategories)
    
    # Non-marker genes
    non_marker_genes <- is.na(df_pc$StageCategories)
    
    # Top genes
    top_genes <- df_pc$TopGene
    
    # Build 2x2 counts
    TopMarkerGenes <- sum(top_genes & marker_genes)
    TopNonMarkerGenes <- sum(top_genes & non_marker_genes)
    OtherMarkerGenes <- sum(!top_genes & marker_genes)
    OtherNonMarkerGenes <- sum(!top_genes & non_marker_genes)
    
    # Check if table is valid
    if ((TopMarkerGenes + OtherMarkerGenes) > 0 & (TopNonMarkerGenes + OtherNonMarkerGenes) > 0) {
      pval <- fisher.test(matrix(c(TopMarkerGenes, TopNonMarkerGenes,
                                   OtherMarkerGenes, OtherNonMarkerGenes), nrow = 2))$p.value
    } else {
      pval <- NA
    }
    
    # Save to list
    results_list[[length(results_list) + 1]] <- data.frame(
      PC = pc,
      StageCategory = stage,
      TopMarkerGenes = TopMarkerGenes,
      TopNonMarkerGenes = TopNonMarkerGenes,
      OtherMarkerGenes = OtherMarkerGenes,
      OtherNonMarkerGenes = OtherNonMarkerGenes,
      p_value = pval
    )
  }
}

# Combine all results
enrichment_results <- do.call(rbind, results_list)

# Adjust p-values
enrichment_results$p_adj <- p.adjust(enrichment_results$p_value, method = "BH")

# Check results
print(enrichment_results)

# PC2 is enriched for adult female marker genes!
# Getting variance explained
var_explained <- pca_res$sdev^2 / sum(pca_res$sdev^2)

# Plotting the PC
pca_scores <- as.data.frame(pca_res$x) %>%
  tibble::rownames_to_column("sample")

# Loading the infection intensity score
infection_intensity = read.csv("infectionIntensity/infectionIntensity.csv")

pca_scores <- pca_scores %>%
  left_join(infection_intensity, by = "sample")

ggplot(pca_scores, aes(x = PC1, y = PC2,
                       color = infection_intensity)) +
  geom_point(size = 3) +
  scale_color_gradientn(colors = wes_palette("Zissou1", 100,
                                          type = "continuous")) +
  labs(title = "PCA of samples colored by infection intensity",
       x = "PC1",
       y = "PC2",
       color = "Infection\nIntensity") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )
View(pca_scores)
# EA02, EA32, FB17 - outliers on PC1!
# I'm going to re-try this after removing these

expr_t <- expr_t[!rownames(expr_t) %in% c("EA02",
                                          "EA32",
                                          "FB17"),] # 80 genotypes

colnames(expr_t) = counts[,1]

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


# Merging the loadings with stage-specific markers
stage_markers_df_clean <- stage_markers_df %>%
  distinct(m_hapla, .keep_all = TRUE)

loading_long <- loading_long %>%
  left_join(stage_markers_df_clean, by = c("geneid" = "m_hapla")) %>%
  mutate(IsMarker = !is.na(StageCategories))

# Get unique PCs and marker stages
pcs <- unique(loading_long$PC)
stages <- unique(loading_long$StageCategories)
stages <- stages[!is.na(stages)]  # remove NA

# Initialize list to hold results
results_list <- list()

# Loop through PCs and stage categories
for (pc in pcs) {
  for (stage in stages) {
    # Subset to this PC
    df_pc <- loading_long[loading_long$PC == pc, ]
    
    # Marker genes for this stage
    marker_genes <- df_pc$StageCategories == stage & !is.na(df_pc$StageCategories)
    
    # Non-marker genes
    non_marker_genes <- is.na(df_pc$StageCategories)
    
    # Top genes
    top_genes <- df_pc$TopGene
    
    # Build 2x2 counts
    TopMarkerGenes <- sum(top_genes & marker_genes)
    TopNonMarkerGenes <- sum(top_genes & non_marker_genes)
    OtherMarkerGenes <- sum(!top_genes & marker_genes)
    OtherNonMarkerGenes <- sum(!top_genes & non_marker_genes)
    
    # Check if table is valid
    if ((TopMarkerGenes + OtherMarkerGenes) > 0 & (TopNonMarkerGenes + OtherNonMarkerGenes) > 0) {
      pval <- fisher.test(matrix(c(TopMarkerGenes, TopNonMarkerGenes,
                                   OtherMarkerGenes, OtherNonMarkerGenes), nrow = 2))$p.value
    } else {
      pval <- NA
    }
    
    # Save to list
    results_list[[length(results_list) + 1]] <- data.frame(
      PC = pc,
      StageCategory = stage,
      TopMarkerGenes = TopMarkerGenes,
      TopNonMarkerGenes = TopNonMarkerGenes,
      OtherMarkerGenes = OtherMarkerGenes,
      OtherNonMarkerGenes = OtherNonMarkerGenes,
      p_value = pval
    )
  }
}

# Combine all results
enrichment_results <- do.call(rbind, results_list)

# Adjust p-values
enrichment_results$p_adj <- p.adjust(enrichment_results$p_value, method = "BH")

# Check results
print(enrichment_results)

# PC1 is enriched for adult female marker genes! p-adj = 1.127894e-26
# Getting variance explained
var_explained <- pca_res$sdev^2 / sum(pca_res$sdev^2)

# Plotting the PC
pca_scores <- as.data.frame(pca_res$x) %>%
  tibble::rownames_to_column("sample")

# Loading the infection intensity score
infection_intensity = read.csv("infectionIntensity/infectionIntensity.csv")

pca_scores <- pca_scores %>%
  left_join(infection_intensity, by = "sample")

plot12 = ggplot(pca_scores, aes(x = PC1, y = PC2,
                       color = infection_intensity)) +
  geom_point(size = 3) +
  scale_color_gradientn(colors = wes_palette("Zissou1", 100,
                                             type = "continuous")) +
  labs(x = "PC1 (21.75% variance)",
       y = "PC2 (9.8% variance)",
       color = "Infection\nIntensity") +
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 14,
                              face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12,
                                face = "bold"),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot12
ggsave("infectionIntensity/20250630PCinfectionIntensity.pdf",
       plot12,
       height = 6, width = 8,
       bg = "transparent",
       dpi = 300)
ggsave("infectionIntensity/20250630PCinfectionIntensity.png",
       plot12,
       height = 6, width = 8,
       bg = "transparent",
       dpi = 300)


# Plotting the expression of some stage marker genes

female_gene = as.data.frame(expr_t[,which(colnames(expr_t) == "MhA1_Contig447.frz3.gene2_mRNA")])
female_gene$sample = row.names(female_gene)
colnames(female_gene) = c("MhA1_Contig447.frz3.gene2", "sample")

pca_scores <- pca_scores %>%
  left_join(female_gene, by = "sample")

plot13 = ggplot(pca_scores, aes(x = PC1, y = PC2,
                                color = MhA1_Contig447.frz3.gene2)) +
  geom_point(size = 3) +
  scale_color_gradientn(colors = wes_palette("Zissou1", 100,
                                             type = "continuous")) +
  labs(x = "PC1 (21.75% variance)",
       y = "PC2 (9.8% variance)",
       color = "MhA1_Contig447.frz3.gene2") +
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 14,
                              face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12,
                                face = "bold"),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot13
  

female_gene = as.data.frame(expr_t[,which(colnames(expr_t) == "MhA1_Contig298.frz3.gene8_mRNA")])
female_gene$sample = row.names(female_gene)
colnames(female_gene) = c("MhA1_Contig298.frz3.gene8", "sample")

pca_scores <- pca_scores %>%
  left_join(female_gene, by = "sample")

plot14 = ggplot(pca_scores, aes(x = PC1, y = PC2,
                                color = MhA1_Contig298.frz3.gene8)) +
  geom_point(size = 3) +
  scale_color_gradientn(colors = wes_palette("Zissou1", 100,
                                             type = "continuous")) +
  labs(x = "PC1 (21.75% variance)",
       y = "PC2 (9.8% variance)",
       color = "MhA1_Contig298.frz3.gene8\n(Female Marker)") +
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 14,
                              face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12,
                                face = "bold"),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot14

ggsave("infectionIntensity/20250917_FemaleMarkerPCA.pdf",
       plot14,
       height = 6, width = 8,
       bg = "transparent",
       dpi = 300)

MhA1_Contig2096.frz3.gene8_mRNA

j4_gene = as.data.frame(expr_t[,which(colnames(expr_t) == "MhA1_Contig2096.frz3.gene8_mRNA")])
j4_gene$sample = row.names(j4_gene)
colnames(j4_gene) = c("MhA1_Contig2096.frz3.gene8", "sample")

pca_scores <- pca_scores %>%
  left_join(j4_gene, by = "sample")

plot15 = ggplot(pca_scores, aes(x = PC1, y = PC2,
                                color = MhA1_Contig2096.frz3.gene8)) +
  geom_point(size = 3) +
  scale_color_gradientn(colors = wes_palette("Zissou1", 100,
                                             type = "continuous")) +
  labs(x = "PC1 (21.75% variance)",
       y = "PC2 (9.8% variance)",
       color = "MhA1_Contig2096.frz3.gene8\n(J4 Marker)") +
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 14,
                              face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12,
                                face = "bold"),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot15

ggsave("infectionIntensity/20250917_J4MarkerPCA.pdf",
       plot15,
       height = 6, width = 8,
       bg = "transparent",
       dpi = 300)

ph4a3 = as.data.frame(expr_t[,which(colnames(expr_t) == "MhA1_Contig2702.frz3.gene1_mRNA")])
ph4a3$sample = row.names(ph4a3)
colnames(ph4a3) = c("Mh-P4HA3", "sample")

pca_scores <- pca_scores %>%
  left_join(ph4a3, by = "sample")

plot15 = ggplot(pca_scores, aes(x = PC1, y = PC2,
                                color = `Mh-P4HA3`)) +
  geom_point(size = 3) +
  scale_color_gradientn(colors = wes_palette("Zissou1", 100,
                                             type = "continuous")) +
  labs(x = "PC1 (21.75% variance)",
       y = "PC2 (9.8% variance)",
       color = "Mhap-P4HA3") +
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 14,
                              face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12,
                                face = "bold"),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot15

ggsave("infectionIntensity/20250917_MhP4HA3PCA.pdf",
       plot15,
       height = 6, width = 8,
       bg = "transparent",
       dpi = 300)

erp3 = as.data.frame(expr_t[,which(colnames(expr_t) == "MhA1_Contig593.frz3.gene3_mRNA")])
erp3$sample = row.names(erp3)
colnames(erp3) = c("Mh-ERP3", "sample")

pca_scores <- pca_scores %>%
  left_join(erp3, by = "sample")

plot15 = ggplot(pca_scores, aes(x = PC1, y = PC2,
                                color = `Mh-ERP3`)) +
  geom_point(size = 3) +
  scale_color_gradientn(colors = wes_palette("Zissou1", 100,
                                             type = "continuous")) +
  labs(x = "PC1 (21.75% variance)",
       y = "PC2 (9.8% variance)",
       color = "Mhap-ERP3") +
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 14,
                              face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12,
                                face = "bold"),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot15

ggsave("infectionIntensity/20250917_MhERP3_PCA.pdf",
       plot15,
       height = 6, width = 8,
       bg = "transparent",
       dpi = 300)

alg2 = as.data.frame(expr_t[,which(colnames(expr_t) == "MhA1_Contig684.frz3.gene16_mRNA")])
alg2$sample = row.names(alg2)
colnames(alg2) = c("Mh-ALG2", "sample")

pca_scores <- pca_scores %>%
  left_join(alg2, by = "sample")

plot15 = ggplot(pca_scores, aes(x = PC1, y = PC2,
                                color = `Mh-ALG2`)) +
  geom_point(size = 3) +
  scale_color_gradientn(colors = wes_palette("Zissou1", 100,
                                             type = "continuous")) +
  labs(x = "PC1 (21.75% variance)",
       y = "PC2 (9.8% variance)",
       color = "Mhap-ALG2") +
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 14,
                              face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12,
                                face = "bold"),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot15

ggsave("infectionIntensity/20250917_MhALG2_PCA.pdf",
       plot15,
       height = 6, width = 8,
       bg = "transparent",
       dpi = 300)

# Collagen genes - rhythmic expression in C elegans molting!
col125 = as.data.frame(expr_t[,which(colnames(expr_t) == "MhA1_Contig1140.frz3.gene4_mRNA")])
col125$sample = row.names(col125)
colnames(col125) = c("Mhap-col-125", "sample")

pca_scores <- pca_scores %>%
  left_join(col125, by = "sample")

plot15 = ggplot(pca_scores, aes(x = PC1, y = PC2,
                                color = `Mhap-col-125`)) +
  geom_point(size = 3) +
  scale_color_gradientn(colors = wes_palette("Zissou1", 100,
                                             type = "continuous")) +
  labs(x = "PC1 (21.75% variance)",
       y = "PC2 (9.8% variance)",
       color = "Mhap-col-125\nCuticle collagen N-terminal\ndomain-containing protein") +
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 14,
                              face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12,
                                face = "bold"),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot15

ggsave("infectionIntensity/20250925_col-125PCA.pdf",
       plot15,
       height = 6, width = 8,
       bg = "transparent",
       dpi = 300)


sqt3 = as.data.frame(expr_t[,which(colnames(expr_t) == "MhA1_Contig1140.frz3.gene4_mRNA")])
sqt3$sample = row.names(sqt3)
colnames(sqt3) = c("Mhap-sqt-3", "sample")

pca_scores <- pca_scores %>%
  left_join(sqt3, by = "sample")

plot15 = ggplot(pca_scores, aes(x = PC1, y = PC2,
                                color = `Mhap-sqt-3`)) +
  geom_point(size = 3) +
  scale_color_gradientn(colors = wes_palette("Zissou1", 100,
                                             type = "continuous")) +
  labs(x = "PC1 (21.75% variance)",
       y = "PC2 (9.8% variance)",
       color = "Mhap-sqt-3\nCuticle collagen 1") +
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 14,
                              face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12,
                                face = "bold"),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot15

ggsave("infectionIntensity/20250925_sqt-3PCA.pdf",
       plot15,
       height = 6, width = 8,
       bg = "transparent",
       dpi = 300)
MhA1_Contig309.frz3.gene32

cut4 = as.data.frame(expr_t[,which(colnames(expr_t) == "MhA1_Contig309.frz3.gene32_mRNA")])
cut4$sample = row.names(cut4)
colnames(cut4) = c("Mhap-cut-4", "sample")

pca_scores <- pca_scores %>%
  left_join(cut4, by = "sample")

plot15 = ggplot(pca_scores, aes(x = PC1, y = PC2,
                                color = `Mhap-cut-4`)) +
  geom_point(size = 3) +
  scale_color_gradientn(colors = wes_palette("Zissou1", 100,
                                             type = "continuous")) +
  labs(x = "PC1 (21.75% variance)",
       y = "PC2 (9.8% variance)",
       color = "Mhap-cut-4\nCuticulin 4") +
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 14,
                              face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12,
                                face = "bold"),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot15

ggsave("infectionIntensity/20250925_cut-4PCA.pdf",
       plot15,
       height = 6, width = 8,
       bg = "transparent",
       dpi = 300)


cutl12 = as.data.frame(expr_t[,which(colnames(expr_t) == "MhA1_Contig1584.frz3.fgene1_mRNA")])
cutl12$sample = row.names(cutl12)
colnames(cutl12) = c("Mhap-cutl-12", "sample")

pca_scores <- pca_scores %>%
  left_join(cutl12, by = "sample")

plot15 = ggplot(pca_scores, aes(x = PC1, y = PC2,
                                color = `Mhap-cutl-12`)) +
  geom_point(size = 3) +
  scale_color_gradientn(colors = wes_palette("Zissou1", 100,
                                             type = "continuous")) +
  labs(x = "PC1 (21.75% variance)",
       y = "PC2 (9.8% variance)",
       color = "Mhap-cutl-12\nCuticulin-like 12") +
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 14,
                              face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12,
                                face = "bold"),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot15

ggsave("infectionIntensity/20250925_cutl-12PCA.pdf",
       plot15,
       height = 6, width = 8,
       bg = "transparent",
       dpi = 300)

twk18 = as.data.frame(expr_t[,which(colnames(expr_t) == "MhA1_Contig25.frz3.gene31_mRNA")])
twk18$sample = row.names(twk18)
colnames(twk18) = c("Mhap-twk-18", "sample")

pca_scores <- pca_scores %>%
  left_join(twk18, by = "sample")

plot15 = ggplot(pca_scores, aes(x = PC1, y = PC2,
                                color = `Mhap-twk-18`)) +
  geom_point(size = 3) +
  scale_color_gradientn(colors = wes_palette("Zissou1", 100,
                                             type = "continuous")) +
  labs(x = "PC1 (21.75% variance)",
       y = "PC2 (9.8% variance)",
       color = "Mhap-twk-18\nTWiK family \npotassium channel 18") +
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 14,
                              face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12,
                                face = "bold"),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot15

ggsave("infectionIntensity/20250925_twk-18PCA.pdf",
       plot15,
       height = 6, width = 8,
       bg = "transparent",
       dpi = 300)

twk18b = as.data.frame(expr_t[,which(colnames(expr_t) == "MhA1_Contig25.frz3.gene33_mRNA")])
twk18b$sample = row.names(twk18b)
colnames(twk18b) = c("Mhap-twk-18b", "sample")

pca_scores <- pca_scores %>%
  left_join(twk18b, by = "sample")

plot15 = ggplot(pca_scores, aes(x = PC1, y = PC2,
                                color = `Mhap-twk-18b`)) +
  geom_point(size = 3) +
  scale_color_gradientn(colors = wes_palette("Zissou1", 100,
                                             type = "continuous")) +
  labs(x = "PC1 (21.75% variance)",
       y = "PC2 (9.8% variance)",
       color = "Mhap-twk-18b\nTWiK family \npotassium channel 18") +
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 14,
                              face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12,
                                face = "bold"),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot15

ggsave("infectionIntensity/20250925_twk-18bPCA.pdf",
       plot15,
       height = 6, width = 8,
       bg = "transparent",
       dpi = 300)

unc9 = as.data.frame(expr_t[,which(colnames(expr_t) == "MhA1_Contig651.frz3.fgene1_mRNA")])
unc9$sample = row.names(unc9)
colnames(unc9) = c("Mhap-unc-9", "sample")

pca_scores <- pca_scores %>%
  left_join(unc9, by = "sample")

plot15 = ggplot(pca_scores, aes(x = PC1, y = PC2,
                                color = `Mhap-unc-9`)) +
  geom_point(size = 3) +
  scale_color_gradientn(colors = wes_palette("Zissou1", 100,
                                             type = "continuous")) +
  labs(x = "PC1 (21.75% variance)",
       y = "PC2 (9.8% variance)",
       color = "Mhap-unc-9\nInnexin") +
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 14,
                              face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12,
                                face = "bold"),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot15

ggsave("infectionIntensity/20250925_unc-9PCA.pdf",
       plot15,
       height = 6, width = 8,
       bg = "transparent",
       dpi = 300)

MhA1_Contig301.frz3.gene22

nas1 = as.data.frame(expr_t[,which(colnames(expr_t) == "MhA1_Contig1377.frz3.gene10_mRNA")])
nas1$sample = row.names(nas1)
colnames(nas1) = c("Mhap-nas-1", "sample")

pca_scores <- pca_scores %>%
  left_join(nas1, by = "sample")

plot15 = ggplot(pca_scores, aes(x = PC1, y = PC2,
                                color = `Mhap-nas-1`)) +
  geom_point(size = 3) +
  scale_color_gradientn(colors = wes_palette("Zissou1", 100,
                                             type = "continuous")) +
  labs(x = "PC1 (21.75% variance)",
       y = "PC2 (9.8% variance)",
       color = "Mhap-nas-1\nZinc metalloproteinase") +
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 14,
                              face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12,
                                face = "bold"),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot15

ggsave("infectionIntensity/20250925_nas-1PCA.pdf",
       plot15,
       height = 6, width = 8,
       bg = "transparent",
       dpi = 300)

slc17.3 = as.data.frame(expr_t[,which(colnames(expr_t) == "MhA1_Contig298.frz3.gene5_mRNA")])
slc17.3$sample = row.names(slc17.3)
colnames(slc17.3) = c("Mhap-slc-17.3", "sample")

pca_scores <- pca_scores %>%
  left_join(slc17.3, by = "sample")

plot15 = ggplot(pca_scores, aes(x = PC1, y = PC2,
                                color = `Mhap-slc-17.3`)) +
  geom_point(size = 3) +
  scale_color_gradientn(colors = wes_palette("Zissou1", 100,
                                             type = "continuous")) +
  labs(x = "PC1 (21.75% variance)",
       y = "PC2 (9.8% variance)",
       color = "Mhap-slc-17.3\nUncharacterized transpoerter") +
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 14,
                              face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12,
                                face = "bold"),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot15

ggsave("infectionIntensity/20250925_slc-17.3PCA.pdf",
       plot15,
       height = 6, width = 8,
       bg = "transparent",
       dpi = 300)

vsra1 = as.data.frame(expr_t[,which(colnames(expr_t) == "MhA1_Contig208.frz3.fgene3_mRNA")])
vsra1$sample = row.names(vsra1)
colnames(vsra1) = c("Mhap-vsra-1", "sample")

pca_scores <- pca_scores %>%
  left_join(vsra1, by = "sample")

plot15 = ggplot(pca_scores, aes(x = PC1, y = PC2,
                                color = `Mhap-vsra-1`)) +
  geom_point(size = 3) +
  scale_color_gradientn(colors = wes_palette("Zissou1", 100,
                                             type = "continuous")) +
  labs(x = "PC1 (21.75% variance)",
       y = "PC2 (9.8% variance)",
       color = "Mhap-vsra-1\nVersatile small RNA/nArgonaute") +
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 14,
                              face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12,
                                face = "bold"),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot15

ggsave("infectionIntensity/20250925_vsra-1PCA.pdf",
       plot15,
       height = 6, width = 8,
       bg = "transparent",
       dpi = 300)


nrde3 = as.data.frame(expr_t[,which(colnames(expr_t) == "MhA1_Contig1305.frz3.gene6_mRNA")])
nrde3$sample = row.names(nrde3)
colnames(nrde3) = c("Mhap-nrde-3", "sample")

pca_scores <- pca_scores %>%
  left_join(nrde3, by = "sample")

plot15 = ggplot(pca_scores, aes(x = PC1, y = PC2,
                                color = `Mhap-nrde-3`)) +
  geom_point(size = 3) +
  scale_color_gradientn(colors = wes_palette("Zissou1", 100,
                                             type = "continuous")) +
  labs(x = "PC1 (21.75% variance)",
       y = "PC2 (9.8% variance)",
       color = "Mhap-nrde-3\nNuclear RNAi Defective 3") +
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 14,
                              face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12,
                                face = "bold"),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot15

ggsave("infectionIntensity/20250925_nrde-3PCA.pdf",
       plot15,
       height = 6, width = 8,
       bg = "transparent",
       dpi = 300)

slc17.4 = as.data.frame(expr_t[,which(colnames(expr_t) == "MhA1_Contig298.frz3.gene5_mRNA")])
slc17.4$sample = row.names(slc17.4)
colnames(slc17.4) = c("Mhap-slc-17.4", "sample")

pca_scores <- pca_scores %>%
  left_join(slc17.4, by = "sample")

plot15 = ggplot(pca_scores, aes(x = PC1, y = PC2,
                                color = `Mhap-slc-17.4`)) +
  geom_point(size = 3) +
  scale_color_gradientn(colors = wes_palette("Zissou1", 100,
                                             type = "continuous")) +
  labs(x = "PC1 (21.75% variance)",
       y = "PC2 (9.8% variance)",
       color = "Mhap-slc-17.4\nPredicted transporter") +
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 14,
                              face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12,
                                face = "bold"),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot15

ggsave("infectionIntensity/20250925_slc-17.4PCA.pdf",
       plot15,
       height = 6, width = 8,
       bg = "transparent",
       dpi = 300)

arf5 = as.data.frame(expr_t[,which(colnames(expr_t) == "MhA1_Contig256.frz3.gene3_mRNA")])
arf5$sample = row.names(arf5)
colnames(arf5) = c("Mhap-arf-5", "sample")

pca_scores <- pca_scores %>%
  left_join(arf5, by = "sample")

plot15 = ggplot(pca_scores, aes(x = PC1, y = PC2,
                                color = `Mhap-arf-5`)) +
  geom_point(size = 3) +
  scale_color_gradientn(colors = wes_palette("Zissou1", 100,
                                             type = "continuous")) +
  labs(x = "PC1 (21.75% variance)",
       y = "PC2 (9.8% variance)",
       color = "Mhap-arf-5\nADP ribosylation factor") +
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 14,
                              face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12,
                                face = "bold"),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot15

ggsave("infectionIntensity/20250925_arf-5PCA.pdf",
       plot15,
       height = 6, width = 8,
       bg = "transparent",
       dpi = 300)

gly6 = as.data.frame(expr_t[,which(colnames(expr_t) == "MhA1_Contig116.frz3.fgene6_mRNA")])
gly6$sample = row.names(gly6)
colnames(gly6) = c("Mhap-gly-6", "sample")

pca_scores <- pca_scores %>%
  left_join(gly6, by = "sample")

plot15 = ggplot(pca_scores, aes(x = PC1, y = PC2,
                                color = `Mhap-gly-6`)) +
  geom_point(size = 3) +
  scale_color_gradientn(colors = wes_palette("Zissou1", 100,
                                             type = "continuous")) +
  labs(x = "PC1 (21.75% variance)",
       y = "PC2 (9.8% variance)",
       color = "Mhap-gly-6\nGalNAc Transferase") +
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 14,
                              face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12,
                                face = "bold"),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot15

ggsave("infectionIntensity/20250925_gly-6PCA.pdf",
       plot15,
       height = 6, width = 8,
       bg = "transparent",
       dpi = 300)

gly6 = as.data.frame(expr_t[,which(colnames(expr_t) == "MhA1_Contig116.frz3.fgene6_mRNA")])
gly6$sample = row.names(gly6)
colnames(gly6) = c("Mhap-gly-6", "sample")

pca_scores <- pca_scores %>%
  left_join(gly6, by = "sample")

plot15 = ggplot(pca_scores, aes(x = PC1, y = PC2,
                                color = `Mhap-gly-6`)) +
  geom_point(size = 3) +
  scale_color_gradientn(colors = wes_palette("Zissou1", 100,
                                             type = "continuous")) +
  labs(x = "PC1 (21.75% variance)",
       y = "PC2 (9.8% variance)",
       color = "Mhap-gly-6\nGalNAc Transferase") +
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 14,
                              face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12,
                                face = "bold"),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot15

ggsave("infectionIntensity/20250925_gly-6PCA.pdf",
       plot15,
       height = 6, width = 8,
       bg = "transparent",
       dpi = 300)

gly11 = as.data.frame(expr_t[,which(colnames(expr_t) == "MhA1_Contig706.frz3.gene4_mRNA")])
gly11$sample = row.names(gly11)
colnames(gly11) = c("Mhap-gly-11", "sample")

pca_scores <- pca_scores %>%
  left_join(gly11, by = "sample")

plot15 = ggplot(pca_scores, aes(x = PC1, y = PC2,
                                color = `Mhap-gly-11`)) +
  geom_point(size = 3) +
  scale_color_gradientn(colors = wes_palette("Zissou1", 100,
                                             type = "continuous")) +
  labs(x = "PC1 (21.75% variance)",
       y = "PC2 (9.8% variance)",
       color = "Mhap-gly-11\nGalNAc Transferase") +
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 14,
                              face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12,
                                face = "bold"),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot15

ggsave("infectionIntensity/20250925_gly-11PCA.pdf",
       plot15,
       height = 6, width = 8,
       bg = "transparent",
       dpi = 300)

gly3 = as.data.frame(expr_t[,which(colnames(expr_t) == "MhA1_Contig1588.frz3.gene15_mRNA")])
gly3$sample = row.names(gly3)
colnames(gly3) = c("Mhap-gly-3", "sample")

pca_scores <- pca_scores %>%
  left_join(gly3, by = "sample")

plot15 = ggplot(pca_scores, aes(x = PC1, y = PC2,
                                color = `Mhap-gly-3`)) +
  geom_point(size = 3) +
  scale_color_gradientn(colors = wes_palette("Zissou1", 100,
                                             type = "continuous")) +
  labs(x = "PC1 (21.75% variance)",
       y = "PC2 (9.8% variance)",
       color = "Mhap-gly-3\nGalNAc Transferase") +
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 14,
                              face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12,
                                face = "bold"),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot15

ggsave("infectionIntensity/20250925_gly-3PCA.pdf",
       plot15,
       height = 6, width = 8,
       bg = "transparent",
       dpi = 300)

mbk2 = as.data.frame(expr_t[,which(colnames(expr_t) == "MhA1_Contig1678.frz3.gene9_mRNA")])
mbk2$sample = row.names(mbk2)
colnames(mbk2) = c("Mhap-mbk-2", "sample")

pca_scores <- pca_scores %>%
  left_join(mbk2, by = "sample")

plot15 = ggplot(pca_scores, aes(x = PC1, y = PC2,
                                color = `Mhap-mbk-2`)) +
  geom_point(size = 3) +
  scale_color_gradientn(colors = wes_palette("Zissou1", 100,
                                             type = "continuous")) +
  labs(x = "PC1 (21.75% variance)",
       y = "PC2 (9.8% variance)",
       color = "Mhap-mbk-2\nDual Ser/Thr &\nTyr kinase") +
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 14,
                              face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12,
                                face = "bold"),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot15

ggsave("infectionIntensity/20250925_mbk-2PCA.pdf",
       plot15,
       height = 6, width = 8,
       bg = "transparent",
       dpi = 300)

cdc48.2 = as.data.frame(expr_t[,which(colnames(expr_t) == "MhA1_Contig1238.frz3.fgene1_mRNA")])
cdc48.2$sample = row.names(cdc48.2)
colnames(cdc48.2) = c("Mhap-cdc-48.2", "sample")

pca_scores <- pca_scores %>%
  left_join(cdc48.2, by = "sample")

plot15 = ggplot(pca_scores, aes(x = PC1, y = PC2,
                                color = `Mhap-cdc-48.2`)) +
  geom_point(size = 3) +
  scale_color_gradientn(colors = wes_palette("Zissou1", 100,
                                             type = "continuous")) +
  labs(x = "PC1 (21.75% variance)",
       y = "PC2 (9.8% variance)",
       color = "Mhap-cdc-48.2\nTransitional ER ATPase") +
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 14,
                              face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12,
                                face = "bold"),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot15

ggsave("infectionIntensity/20250925_cdc-18.2PCA.pdf",
       plot15,
       height = 6, width = 8,
       bg = "transparent",
       dpi = 300)

rpt3 = as.data.frame(expr_t[,which(colnames(expr_t) == "MhA1_Contig302.frz3.gene4_mRNA")])
rpt3$sample = row.names(rpt3)
colnames(rpt3) = c("Mhap-rpt-3", "sample")

pca_scores <- pca_scores %>%
  left_join(rpt3, by = "sample")

plot15 = ggplot(pca_scores, aes(x = PC1, y = PC2,
                                color = `Mhap-rpt-3`)) +
  geom_point(size = 3) +
  scale_color_gradientn(colors = wes_palette("Zissou1", 100,
                                             type = "continuous")) +
  labs(x = "PC1 (21.75% variance)",
       y = "PC2 (9.8% variance)",
       color = "Mhap-rpt-3\n26S proteasome regulatory\nsubunit 6B") +
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 14,
                              face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12,
                                face = "bold"),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot15

ggsave("infectionIntensity/20250925_rpt-3PCA.pdf",
       plot15,
       height = 6, width = 8,
       bg = "transparent",
       dpi = 300)

rpt5 = as.data.frame(expr_t[,which(colnames(expr_t) == "MhA1_Contig596.frz3.gene13_mRNA")])
rpt5$sample = row.names(rpt5)
colnames(rpt5) = c("Mhap-rpt-5", "sample")

pca_scores <- pca_scores %>%
  left_join(rpt5, by = "sample")

plot15 = ggplot(pca_scores, aes(x = PC1, y = PC2,
                                color = `Mhap-rpt-5`)) +
  geom_point(size = 3) +
  scale_color_gradientn(colors = wes_palette("Zissou1", 100,
                                             type = "continuous")) +
  labs(x = "PC1 (21.75% variance)",
       y = "PC2 (9.8% variance)",
       color = "Mhap-rpt-5\n26S proteasome regulatory\nsubunit 6A") +
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 14,
                              face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12,
                                face = "bold"),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot15

ggsave("infectionIntensity/20250925_rpt-5PCA.pdf",
       plot15,
       height = 6, width = 8,
       bg = "transparent",
       dpi = 300)

pgp2 = as.data.frame(expr_t[,which(colnames(expr_t) == "MhA1_Contig76.frz3.fgene1_mRNA")])
pgp2$sample = row.names(pgp2)
colnames(pgp2) = c("Mhap-pgp-2", "sample")

pca_scores <- pca_scores %>%
  left_join(pgp2, by = "sample")

plot15 = ggplot(pca_scores, aes(x = PC1, y = PC2,
                                color = `Mhap-pgp-2`)) +
  geom_point(size = 3) +
  scale_color_gradientn(colors = wes_palette("Zissou1", 100,
                                             type = "continuous")) +
  labs(x = "PC1 (21.75% variance)",
       y = "PC2 (9.8% variance)",
       color = "Mhap-pgp-2\nATPase-coupled transporter") +
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 14,
                              face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12,
                                face = "bold"),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot15

ggsave("infectionIntensity/20250925_pgp-2PCA.pdf",
       plot15,
       height = 6, width = 8,
       bg = "transparent",
       dpi = 300)

znf782 = as.data.frame(expr_t[,which(colnames(expr_t) == "MhA1_Contig193.frz3.gene15_mRNA")])
znf782$sample = row.names(znf782)
colnames(znf782) = c("Mhap-znf-782", "sample")

pca_scores <- pca_scores %>%
  left_join(znf782, by = "sample")

plot15 = ggplot(pca_scores, aes(x = PC1, y = PC2,
                                color = `Mhap-znf-782`)) +
  geom_point(size = 3) +
  scale_color_gradientn(colors = wes_palette("Zissou1", 100,
                                             type = "continuous")) +
  labs(x = "PC1 (21.75% variance)",
       y = "PC2 (9.8% variance)",
       color = "Mhap-znf-782\nZinc finger TF") +
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 14,
                              face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12,
                                face = "bold"),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot15

ggsave("infectionIntensity/20250925_znf-782PCA.pdf",
       plot15,
       height = 6, width = 8,
       bg = "transparent",
       dpi = 300)

mccc1 = as.data.frame(expr_t[,which(colnames(expr_t) == "MhA1_Contig354.frz3.gene72_mRNA")])
mccc1$sample = row.names(mccc1)
colnames(mccc1) = c("Mhap-mccc-1", "sample")

pca_scores <- pca_scores %>%
  left_join(mccc1, by = "sample")

plot15 = ggplot(pca_scores, aes(x = PC1, y = PC2,
                                color = `Mhap-mccc-1`)) +
  geom_point(size = 3) +
  scale_color_gradientn(colors = wes_palette("Zissou1", 100,
                                             type = "continuous")) +
  labs(x = "PC1 (21.75% variance)",
       y = "PC2 (9.8% variance)",
       color = "Mhap-mccc-1\nMethylcrotonyl-CoA carboxylase\nsubunit alpha") +
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 14,
                              face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12,
                                face = "bold"),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot15

ggsave("infectionIntensity/20250925_mccc-1PCA.pdf",
       plot15,
       height = 6, width = 8,
       bg = "transparent",
       dpi = 300)

MhA1_Contig1126.frz3.gene3

oac46 = as.data.frame(expr_t[,which(colnames(expr_t) == "MhA1_Contig512.frz3.gene22_mRNA")])
oac46$sample = row.names(oac46)
colnames(oac46) = c("Mhap-oac-46", "sample")

pca_scores <- pca_scores %>%
  left_join(oac46, by = "sample")

plot15 = ggplot(pca_scores, aes(x = PC1, y = PC2,
                                color = `Mhap-oac-46`)) +
  geom_point(size = 3) +
  scale_color_gradientn(colors = wes_palette("Zissou1", 100,
                                             type = "continuous")) +
  labs(x = "PC1 (21.75% variance)",
       y = "PC2 (9.8% variance)",
       color = "Mhap-oac-46\nO-Acyltransferase") +
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 14,
                              face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12,
                                face = "bold"),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot15

ggsave("infectionIntensity/20250925_oac-46PCA.pdf",
       plot15,
       height = 6, width = 8,
       bg = "transparent",
       dpi = 300)

ttll5 = as.data.frame(expr_t[,which(colnames(expr_t) == "MhA1_Contig309.frz3.gene23_mRNA")])
ttll5$sample = row.names(ttll5)
colnames(ttll5) = c("Mhap-ttll-5", "sample")

pca_scores <- pca_scores %>%
  left_join(ttll5, by = "sample")

plot15 = ggplot(pca_scores, aes(x = PC1, y = PC2,
                                color = `Mhap-ttll-5`)) +
  geom_point(size = 3) +
  scale_color_gradientn(colors = wes_palette("Zissou1", 100,
                                             type = "continuous")) +
  labs(x = "PC1 (21.75% variance)",
       y = "PC2 (9.8% variance)",
       color = "Mhap-ttll-5\nTubulin polyglutamylase") +
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 14,
                              face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12,
                                face = "bold"),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot15

ggsave("infectionIntensity/20250925_ttll-5PCA.pdf",
       plot15,
       height = 6, width = 8,
       bg = "transparent",
       dpi = 300)

dpy9 = as.data.frame(expr_t[,which(colnames(expr_t) == "MhA1_Contig206.frz3.gene2_mRNA")])
ttll5$sample = row.names(ttll5)
colnames(ttll5) = c("Mhap-ttll-5", "sample")

pca_scores <- pca_scores %>%
  left_join(ttll5, by = "sample")

plot15 = ggplot(pca_scores, aes(x = PC1, y = PC2,
                                color = `Mhap-ttll-5`)) +
  geom_point(size = 3) +
  scale_color_gradientn(colors = wes_palette("Zissou1", 100,
                                             type = "continuous")) +
  labs(x = "PC1 (21.75% variance)",
       y = "PC2 (9.8% variance)",
       color = "Mhap-ttll-5\nTubulin polyglutamylase") +
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 14,
                              face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12,
                                face = "bold"),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot15

ggsave("infectionIntensity/20250925_ttll-5PCA.pdf",
       plot15,
       height = 6, width = 8,
       bg = "transparent",
       dpi = 300)

eigens = read.csv("coXnetwork/moduleEigen.csv")

pca_scores <- pca_scores %>%
  left_join(eigens, by = "sample")

plot15 = ggplot(pca_scores, aes(x = PC1, y = PC2,
                                color = MEred)) +
  geom_point(size = 3) +
  scale_color_gradientn(colors = wes_palette("Zissou1", 100,
                                             type = "continuous")) +
  labs(x = "PC1 (21.75% variance)",
       y = "PC2 (9.8% variance)",
       color = "MEred") +
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 14,
                              face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12,
                                face = "bold"),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot15

ggsave("coXnetwork/20250925_MEredPCA.pdf",
       plot15,
       height = 6, width = 8,
       bg = "transparent",
       dpi = 300)

plot15 = ggplot(pca_scores, aes(x = PC1, y = PC2,
                                color = MEbrown)) +
  geom_point(size = 3) +
  scale_color_gradientn(colors = wes_palette("Zissou1", 100,
                                             type = "continuous")) +
  labs(x = "PC1 (21.75% variance)",
       y = "PC2 (9.8% variance)",
       color = "MEbrown") +
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 14,
                              face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12,
                                face = "bold"),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot15

## WHOAH.....

ggsave("coXnetwork/20250925_MEbrownPCA.pdf",
       plot15,
       height = 6, width = 8,
       bg = "transparent",
       dpi = 300)

plot15 = ggplot(pca_scores, aes(x = PC1, y = PC2,
                                color = MEblack)) +
  geom_point(size = 3) +
  scale_color_gradientn(colors = wes_palette("Zissou1", 100,
                                             type = "continuous")) +
  labs(x = "PC1 (21.75% variance)",
       y = "PC2 (9.8% variance)",
       color = "MEblack") +
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 14,
                              face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12,
                                face = "bold"),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot15

## Chop...

ggsave("coXnetwork/20250925_MEblackPCA.pdf",
       plot15,
       height = 6, width = 8,
       bg = "transparent",
       dpi = 300)

plot15 = ggplot(pca_scores, aes(x = PC1, y = PC2,
                                color = MEyellow)) +
  geom_point(size = 3) +
  scale_color_gradientn(colors = wes_palette("Zissou1", 100,
                                             type = "continuous")) +
  labs(x = "PC1 (21.75% variance)",
       y = "PC2 (9.8% variance)",
       color = "MEyellow") +
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 14,
                              face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12,
                                face = "bold"),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot15

## Kind of a ki...

ggsave("coXnetwork/20250925_MEyellowPCA.pdf",
       plot15,
       height = 6, width = 8,
       bg = "transparent",
       dpi = 300)

# This looks like a developmental coninuum! Can we order these sampels by pseudotime?
coords <- as.matrix(pca_scores[, c("PC1", "PC2")])
fit <- principal_curve(coords)

pca_scores$Pseudotime <- fit$lambda

plot13 = ggplot(pca_scores, aes(x = PC1, y = PC2, color = Pseudotime)) +
  geom_point(size = 3) +
  scale_color_viridis_c() +
  theme_linedraw() +
  labs(title = "Samples ordered by principal curve pseudotime",
       x = "PC1 (21.75% variance)",
       y = "PC2 (9.8% variance)",
       color = "Pseudotime") +
  theme(
    axis.title = element_text(size = 14,
                              face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12,
                                face = "bold"),
    legend.text = element_text(size = 10)
  )
plot13

ggsave("infectionIntensity/20250701PCPseudotime.pdf",
       plot13,
       height = 6, width = 8,
       bg = "transparent",
       dpi = 300)

# I want to compare infection intensity and pseudotime scores
plot14 = ggplot(pca_scores,
       aes(x = Pseudotime,
           y = infection_intensity)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "sienna3", linetype = "dashed") +
  labs(x = "Pseudotime Score",
       y = "Infection Intensity") +
  theme_linedraw() +
  theme(
    axis.title = element_text(size = 14,
                              face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12,
                                face = "bold"),
    legend.text = element_text(size = 10)
  )
plot14

ggsave("infectionIntensity/20250701PseudotimeVsInfectionIntensity.pdf",
       plot14,
       height = 6, width = 6,
       bg = "transparent",
       dpi = 300)

cor.test(pca_scores$Pseudotime, pca_scores$infection_intensity,
         method = "spearman")
#Spearman's rank correlation rho
#data:  pca_scores$Pseudotime and pca_scores$infection_intensity
#S = 17932, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho 
#0.7898265 

# FIGURE - ODDS RATIO
pc_enrich <- enrichment_results %>%
  mutate(
    PC_num = as.numeric(gsub("PC", "", PC)),
    StageCategory = case_when(
      StageCategory == "adult.female" ~ "Adult female",
      StageCategory == "J2.juveniles" ~ "J2",
      StageCategory == "Constitutive" ~ "Broad",
      TRUE ~ StageCategory
    ),
    StageCategory = factor(
      StageCategory,
      levels = c("Broad", "Egg", "J2", "J3", "J4", "Adult female")
    ),
    odds_ratio = (TopMarkerGenes * OtherNonMarkerGenes) /
      (TopNonMarkerGenes * OtherMarkerGenes),
    log2_or = log2(odds_ratio),
    neglog10_padj = -log10(p_adj)
  )


pc_top10 <- pc_enrich %>%
  filter(PC %in% c("PC1", "PC2", "PC3", "PC4", "PC5",
                   "PC6", "PC7", "PC8", "PC9", "PC10")) 
pc_top10$PC = factor(pc_top10$PC,
                     levels = c("PC1", "PC2", "PC3", "PC4", "PC5",
                                "PC6", "PC7", "PC8", "PC9", "PC10")) 

plot_heatmap <- ggplot(pc_top10, aes(x = StageCategory, y = PC, fill = log2_or)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(
    aes(label = ifelse(p_adj < 0.05, "*", "")),
    size = 5,
    fontface = "bold"
  ) +
  scale_fill_gradientn(
    colors = colorRampPalette(c(
      "#00441b", "#1b7837", "#5aae61", "#a6dba0",
      "#d9f0d3", "#f7f7f7", "#e7d4e8", "#c2a5cf",
      "#9970ab", "#762a83", "#40004b"
    ))(101),
    limits = c(-3, 3),
    oob = scales::squish,
    name = "log2\n(odds ratio)"
  ) +
  labs(
    x = "Transferred life-stage label",
    y = "Principal component"
  ) +
  theme_tyler() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

plot_heatmap

ggsave("14_life_stages/pc_life_stage_enrichment.pdf",
       plot_heatmap,
       height = 6, width = 5)
