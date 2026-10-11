#### Univariate Selection on Gene Expression in M hapla RIL Population 

library(ggplot2)

### SETWD ###
setwd("~/Desktop/PhD/04_projects/04_mhapla_infection_intensity")
source("00_R/Plot_utils.R")


### IMPORT DATASETS ###
infection = read.csv("03_infection_intensity/infectionIntensity.csv")
expr = read.csv("04_feature_counts/parsed_hapla_log2normCPM.csv")
expr <- expr[,!colnames(expr) %in% c("EA02",
                                     "EA32",
                                     "FB17")] # 80 genotypes

### PREPROCESS THE DATA ###
geneids = expr[,1]
# Transposing and scaling the expression data
expr_scaled = scale(t(expr[,2:81]))
expr_scaled = as.data.frame(expr_scaled)
# Adding sample name for dataset merging
expr_scaled$sample = row.names(expr_scaled)

# Normalizing infection intensity
infection = infection[,c(1,9)]
infection$norm_infection = infection$infection_intensity / mean(infection$infection_intensity)

# Merging the datasets
dt_merge = merge(expr_scaled, infection,
                 by = "sample")


### RUNNING UNIVARIATE SELECTION MODELS ###
# Initialize variables to store results
model = data.frame(S = numeric(0),
                   S.error = numeric(0),
                   S.sig = numeric(0),
                   S.quad = numeric(0),
                   S.quad.err = numeric(0),
                   S.quad.sig = numeric(0),
                   C = numeric(0),
                   C.error = numeric(0),
                   C.sig = numeric(0))


# Loop through expression_scaled and run the linear model
for(i in dt_merge[,2:(ncol(dt_merge)-1)]){ # omitting the last column, with fitness trait
  m <- glm(norm_infection ~ i, data=dt_merge,
           family = gaussian) # Running the linear model first
  sum <- summary(m) # Extracting summary of the model
  S = sum$coefficients[2,1]
  S.err = sum$coefficients[2,2]
  S.sig = sum$coefficients[2,4]
  
  m <- glm(norm_infection ~ i + I(i^2), data=dt_merge,
           family = gaussian) # Running the quadratic model second
  sum <- summary(m) 
  S.quad = sum$coefficients[2,1]
  S.quad.err = sum$coefficients[2,2]
  S.quad.sig = sum$coefficients[2,4]
  C = sum$coefficients[3,1]
  C.err = sum$coefficients[3,2]
  C.sig = sum$coefficients[3,4]
  
  model = rbind(model,
                data.frame(S = S,
                           S.err = S.err,
                           S.sig = S.sig,
                           S.quad = S.quad,
                           S.quad.err = S.quad.err,
                           S.quad.sig = S.quad.sig,
                           C = C,
                           C.err = C.err,
                           C.sig = C.sig))
}

# Accidentally ran the last column which was the non-normalized infection intensity
model = model[1:8685,]

# Adding back in gene IDs
model$geneid = geneids

# Multiple testing corrections
model$S.padj = p.adjust(model$S.sig,
                        method = "bonferroni")
model$C.padj = p.adjust(model$C.sig,
                        method = "bonferroni")

# Saving to file
write.csv(model,
          "07_univariate_model_hapla/hapla_selectionModel.csv")
model = read.csv("07_univariate_model_hapla/hapla_selectionModel.csv")

### EXPLORING THE DISTRIBUTION OF SELECTION DIFFERENTIALS ###
# Range of S-values
range(model$S) # -0.4720740 to 0.4703802
median(model$S) # S(median) = 0.01377137 overall
median(abs(model$S)) # S(median) = 0.1787922 overall
pos = model[which(model$S > 0),] #4,153 transcripts with positive selection
median(pos$S) # S(median) = 0.1844354 for positive selection
neg = model[which(model$S < 0),] #4,532 transcripts with negative selection
median(neg$S) # S(median) = -0.174506 for negative selection
ks = ks.test(pos$S,
             abs(neg$S))
print(ks) # D = 0.046151, p-value = 0.0001958 NO BIAS TOWARD NEGATIVE SELECTION
# Asymptotic two-sided Kolmogorov-Smirnov test
plot1 <- ggplot(model, aes(x = S)) +
  geom_histogram(
    bins = 60,
    fill = "#533943",
    color = "white",
  ) +
  geom_vline(
    xintercept = 0,
    color = "gray25",
    linetype = "dashed",
    linewidth = 0.8
  ) +
  labs(
    y = "Number parasite transcripts",
    x = expression(bold("Interaction differential (" * S[p] * ")"))
  ) +
  theme_tyler()+
  theme(
    plot.margin = margin(t = 10, r = 25, b = 10, l = 10)
  )
plot1
ggsave("07_univariate_model_hapla/histogramS.pdf",
       plot1,
       dpi = 300,
       width = 4, height = 4,
       bg = "transparent")



# Range of C-values
range(model$C) # -0.8772070 to 0.9792685
median(model$C) # C(median)  -0.02061678 overall
div = model[which(model$C > 0),] # 3,409 transcripts with diversifying selection
median(div$C) # C(median) = 0.04736512 for diversifying selection
stab = model[which(model$C < 0),] #5,276 transcripts with stabilizing selection
median(stab$C) # C(median) = -0.05536337 for stabilizing selection
ut = wilcox.test(div$C,
                 abs(stab$C))
print(ut) # W = 8079645, p = 1.201e-15 TREND FOR STABILIZING SELECTION

plot2 = ggplot(model,
               aes(x = C))
plot2 = plot2 + geom_histogram(binwidth = 0.01) +
  theme_linedraw() +
  geom_vline(xintercept = 0,
             color = "sienna3", linetype = "dashed",
             size = 0.8) +  
  labs(y = "Number of transcripts") +
  theme(
    axis.text = element_text(size = 12),        # Size of axis tick labels
    axis.title = element_text(size = 14),       # Size of axis titles
    axis.title.x = element_text(size = 14, face = "bold"), # Custom style for x-axis title
    axis.title.y = element_text(size = 14, face = "bold"), # Custom style for y-axis title), # Custom style for y-axis title
  ) 
plot2
ggsave("selectionModel/histogramC.pdf",
       plot2,
       dpi = 300,
       width = 8, height = 6,
       bg = "transparent")

## HOW MANY ARE SIGNIFICANT? ##
sig_lin = model[which(model$S.sig <= 0.01),] #4,423 significant before multiple test correction
sig_lin_adj = model[which(model$S.padj <= 0.01),] # 1,673 significant after multiple test correction
write.csv(sig_lin_adj,
          "selectionModel/sig_lin_padj0.01.csv")
# Let's rank the strongest positive and negative linear selection
sig_lin_rank <- sig_lin_adj %>%
  arrange(desc(S))
write.csv(sig_lin_rank,
          "selectionModel/sig_lin_padj0.01_ranked.csv")

sig_quad = model[which(model$C.sig <= 0.01),] # 626 transcripts
sig_quad_adj = model[which(model$C.padj <= 0.01),] # 3 transcripts
sig_div = sig_quad_adj[which(sig_quad_adj$C > 0),] # 1 under diversifying selection
write.csv(sig_div,
          "selectionModel/sig_div_padj0.01.csv")
sig_stab = sig_quad_adj[which(sig_quad_adj$C < 0),] # 2 under stabilizing selection
write.csv(sig_stab,
          "selectionModel/sig_stab_padj0.01.csv")

### MAKING VOLCANO PLOTS ###
model <- model %>%
  mutate(
    log10_p = -log10(S.sig),
    sig = case_when(
      S.padj < 0.01 & S > 0 ~ "Positive selection",
      S.padj < 0.01 & S < 0 ~ "Negative selection",
      TRUE ~ "Not significant"
    )
  )

# Volcano plot
plot2 = ggplot(model,
               aes(x = S, y = log10_p,
                   color = sig)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c(
    "Positive selection" = "tomato1",
    "Negative selection" = "steelblue1",
    "Not significant" = "grey70"
  )) +
  labs(
    x = "Selection Differential (S)",
    y = "-log10(adj. p-value)",
    color = "Selection"
  ) +
  theme_linedraw(base_size = 14) +
  theme(legend.text = element_text(size = 12),
        legend.position = "right",
        legend.title = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14,
                                  face = "bold"))
plot2
ggsave("selectionModel/volcanoS.pdf",
       plot2,
       dpi = 300,
       bg = "transparent",
       width = 8, height = 6)

# Now for C
model = model[which(model$C.padj > 0.01),]
model$sig = "Not Significant"
sig_div$sig = "Diversifying Selection"
model = rbind(model, sig_div)
sig_stab$sig = "Stabilizing Selection"
model = rbind(model, sig_stab)
model$log10_p = -log10(model$C.sig)

# Volcano plot
plot3 = ggplot(model,
               aes(x = C, y = log10_p,
                   color = sig)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c(
    "Diversifying Selection" = "purple1",
    "Stabilizing Selection" = "springgreen3",
    "Not significant" = "grey70"
  )) +
  labs(
    x = "Selection Differential (C)",
    y = "-log10(adj. p-value)",
    color = "Selection"
  ) +
  theme_linedraw(base_size = 14) +
  theme(legend.text = element_text(size = 12),
        legend.position = "right",
        legend.title = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14,
                                  face = "bold"))
plot3
ggsave("selectionModel/volcanoC.png",
       plot3,
       dpi = 300,
       bg = "transparent",
       width = 8, height = 6)

# Genes enriched in dorsal gland
dorsalGenes = read.csv("00_source_data/rocha_DGenrich_MHapla_Orthogroups.csv")


model$geneid = gsub("_mRNA", "",
                    model$geneid)
model = model %>%
  left_join(dorsalGenes,
            by = c("geneid" = "m_hapla"))

model$DG_enriched <- ifelse(model$geneid %in% dorsalGenes$m_hapla,
                            "Dorsal gland", "Female head")

dorsal = model[which(model$DG_enriched == "Dorsal gland"),]
nonenrich = model[which(model$DG_enriched == "Female head"),]

wilcox.test(dorsal$S, nonenrich$S)
#Wilcoxon rank sum test with continuity correction
#data:  dorsal$S and nonenrich$S
#W = 169077, p-value < 2.2e-16


wilcox.test(abs(dorsal$S), abs(nonenrich$S))
#Wilcoxon rank sum test with continuity correction
#data:  dorsal$S and nonenrich$S
#W = 662409, p-value < 2.2e-16

plot3 = ggplot(model,
               aes(x =DG_enriched,  y = abs(S),
                   fill = DG_enriched)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = wes_palette("Zissou1", 2, "continuous")) +
  labs(
    y = "Interaction differential (|S|)",
    x = "Expression Enrichment",
    color = "Selection"
  ) +
  theme_tyler() +
  theme(axis.text.x = element_text(angle = 30,
                                   vjust = .5),
        axis.title.x = element_blank(),
        legend.position = "none")
plot3
ggsave("07_univariate_model_hapla/dorsalGland_absS.png",
       plot3,
       dpi = 300,
       bg = "transparent",
       width = 3, height = 4)

wilcox.test(dorsal$C, nonenrich$C)
#Wilcoxon rank sum test with continuity correction
#data:  dorsal$S and nonenrich$S
#W = 169077, p-= 0.0002912

summary_df <- model %>%
  group_by(DG_enriched) %>%
  summarise(
    mean_C = mean(C, na.rm = TRUE),
    se_C = sd(C, na.rm = TRUE) / sqrt(n())
  )

# Create plot
plot4 <- ggplot(summary_df,
                aes(x = DG_enriched, y = mean_C, fill = DG_enriched)) +
  geom_bar(stat = "identity", alpha = 0.7, width = 0.6) +
  geom_errorbar(aes(ymin = mean_C - se_C, ymax = mean_C + se_C), width = 0.2) +
  scale_fill_manual(values = wes_palette("Zissou1", 2, "continuous")) +
  labs(
    y = "Selection Differential (C)",
    x = "Dorsal Gland Enrichment"
  ) +
  theme_linedraw(base_size = 14) +
  theme(
    legend.text = element_text(size = 12),
    legend.position = "right",
    legend.title = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold")
  )

plot4
ggsave("lifeStage/dorsalGland_C.pdf",
       plot4,
       dpi = 300,
       bg = "transparent",
       width = 6, height = 6)

colnames(dt_merge) = c("Sample", geneids, "InfectionIntensity")

alg2 = dt_merge[,which(colnames(dt_merge) %in% c("InfectionIntensity", "MhA1_Contig684.frz3.gene16_mRNA"))]

plot5 = ggplot(alg2,
               aes(x = MhA1_Contig684.frz3.gene16_mRNA,
                   y = InfectionIntensity)) +
  geom_point() +
  geom_smooth(method = "lm",
              color = "sienna3") +
  theme_linedraw() +
  labs(
    y = "Infection Intensity",
    x = "MhALG2 (Scaled Expression)"
  ) +
  theme(
    legend.text = element_text(size = 12),
    legend.position = "right",
    legend.title = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  ) 
plot5
ggsave("selectionModel/MhALG2Expr.pdf",
       plot5,
       dpi = 300,
       width = 8, height = 6,
       bg = "transparent")



# Plotting ALG2 as a volcano plot
model <- model %>%
  mutate(color = ifelse(geneid == "MhA1_Contig684.frz3.gene16_mRNA", "highlight", "other"))


plot = ggplot(model, aes(x = S, y = -log10(S.sig), color = color)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("other" = "gray70", "highlight" = "sienna3")) +
  theme_linedraw() +
  labs(x = "S (Interaction Coeffiecient)", y = "-log10(p-value)") +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA)) +
  # optional: label the gene of interest
  geom_text(data = subset(model, geneid == "MhA1_Contig684.frz3.gene16_mRNA"),
            aes(label = "Mh-ALG2"), vjust = -1, color = "sienna3")
ggsave("selectionModel/MhALG2Volcano.pdf",
       plot,
       dpi = 300,
       width = 6, height = 6,
       bg = "transparent")
