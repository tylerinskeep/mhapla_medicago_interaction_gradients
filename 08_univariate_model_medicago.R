#### Univariate Selection on Gene Expression in Medicago

setwd("~/Desktop/PhD/04_projects/04_mhapla_infection_intensity")
source("00_R/plot_utils.R")


### IMPORT DATASETS ###
infection = read.csv("03_infection_intensity/infectionIntensity.csv")
expr = read.csv("04_feature_counts/mtrunc_log2normCPM.csv")


### PREPROCESS THE DATA ###
geneids = expr[,1]
# Transposing and scaling the expression data
expr_scaled = scale(t(expr[,2:96]))
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
model = model[1:18361,]

# Adding back in gene IDs
model$geneid = geneids

# Multiple testing corrections
model$S.padj = p.adjust(model$S.sig,
                        method = "bonferroni")
model$C.padj = p.adjust(model$C.sig,
                        method = "bonferroni")

# Saving to file
write.csv(model,
          "08_univariate_model_medicago/medicagoModel.csv")

model = read.csv( "08_univariate_model_medicago/medicagoModel.csv")

### EXPLORING THE DISTRIBUTION OF SELECTION DIFFERENTIALS ###
# Range of S-values
range(model$S) # -0.3999910 to 0.4556336
median(model$S) # S(median) = 0.00882929 overall
median(abs(model$S)) # S(median) = 0.10312 overall
pos = model[which(model$S > 0),] #9,578 transcripts with positive selection
median(pos$S) # S(median) = 0.1029774 for positive selection
neg = model[which(model$S < 0),] #8,783 transcripts with negative selection
median(neg$S) # S(median) = -0.1031434 for negative selection
ks = ks.test(pos$S,
             abs(neg$S))
print(ks) # D = 0.022281, p-value = 0.02116 BIAS FOR NEGATIVE
# Higher expression of many genes - maybe better for host?
# Asymptotic two-sided Kolmogorov-Smirnov test

plot1 <- ggplot(model, aes(x = S)) +
  geom_histogram(
    bins = 60,
    fill = "#3A6547",
    color = "white",
  ) +
  geom_vline(
    xintercept = 0,
    color = "gray25",
    linetype = "dashed",
    linewidth = 0.8
  ) +
  labs(
    y = "Number host transcripts",
    x = expression(bold("Interaction differential (" * S[h] * ")"))
  ) +
  theme_tyler() +
  theme(
    plot.margin = margin(t = 10, r = 25, b = 10, l = 10)
  )

plot1
ggsave("08_univariate_model_medicago/histogramS.png",
       plot1,
       dpi = 300,
       width = 4, height = 4,
       bg = "transparent")

# Range of C-values
range(model$C) # -0.2307136 to 0.2495333
median(model$C) # C(median)  -0.02356472 overall
div = model[which(model$C > 0),] # 6,086 transcripts with diversifying selection
median(div$C) # C(median) = 0.02793949 for diversifying selection
stab = model[which(model$C < 0),] #12,275 transcripts with stabilizing selection
median(stab$C) # C(median) = -0.04614407 for stabilizing selection
ut = wilcox.test(div$C,
                 abs(stab$C))
print(ut) # W = 26411551, p <2.2e-16 TREND FOR STABILIZING SELECTION

plot2 = ggplot(model,
               aes(x = C))
plot2 = plot2 + geom_histogram(binwidth = 0.001) +
  theme_linedraw() +
  geom_vline(xintercept = 0,
             color = "springgreen4", linetype = "dashed",
             size = 0.8) +  
  labs(y = "Number of transcripts") +
  theme(
    axis.text = element_text(size = 12),        # Size of axis tick labels
    axis.title = element_text(size = 14),       # Size of axis titles
    axis.title.x = element_text(size = 14, face = "bold"), # Custom style for x-axis title
    axis.title.y = element_text(size = 14, face = "bold"), # Custom style for y-axis title), # Custom style for y-axis title
  ) 
plot2
ggsave("selectionModel/medicago/histogramC.png",
       plot2,
       dpi = 300,
       width = 8, height = 6,
       bg = "transparent")

## HOW MANY ARE SIGNIFICANT? ##
sig_lin = model[which(model$S.sig <= 0.01),] #4,990 significant before multiple test correction
sig_lin_adj = model[which(model$S.padj <= 0.01),] # 357 significant after multiple test correction
write.csv(sig_lin_adj,
          "selectionModel/medicago/sig_lin_padj0.01.csv")
# Let's rank the strongest positive and negative linear selection
sig_lin_rank <- sig_lin_adj %>%
  arrange(desc(S))
write.csv(sig_lin_rank,
          "selectionModel/medicago/sig_lin_padj0.01_ranked.csv")

sig_quad = model[which(model$C.sig <= 0.01),] # 830 transcripts
sig_quad_adj = model[which(model$C.padj <= 0.01),] # 0 transcripts


### MAKING EXPRESSION COVARIANCE PLOTS ###
model = read.csv("selectionModel/medicago/medicagoModel.csv")
MtOPR2.L = "rna-XM_003610721.4"
MtCESA8.L = "rna-XM_013607142.2"
MtCESA4.L = "rna-XM_003594834.4"
MtIAA14.L = "rna-XM_003591758.4"
GDSL = "rna-XM_003591011.4"
# I think that these are all showing an enrichment for orthologs of AtVND7 target genes
# VND1-3 regulate xylem development in galls
MtPLT2.L = "rna-XM_003597400.4"
MtLBD16.L = "rna-XM_013594406.3"

infection = read.csv("infectionIntensity/infectionIntensity.csv")
expr = read.csv("featureCounts/mtrunc_log2normCPM.csv")

### PREPROCESS THE DATA ###
geneids = expr[,1]
# Transposing and scaling the expression data
expr_scaled = scale(t(expr[,2:96]))
expr_scaled = as.data.frame(expr_scaled)
# Adding sample name for dataset merging
expr_scaled$sample = row.names(expr_scaled)
colnames(expr_scaled) = geneids

# Normalizing infection intensity
infection = infection[,c(1,9)]
infection$norm_infection = infection$infection_intensity / mean(infection$infection_intensity)

# Merging the datasets
dt_merge = merge(expr_scaled, infection,
                 by = "sample")
dt_merge = as.data.frame(dt_merge)

opr2 = as.data.frame(expr_scaled[,which(colnames(expr_scaled) %in% MtOPR2.L)])
opr2$MtOPR2.L = opr2$`expr_scaled[, which(colnames(expr_scaled) %in% MtOPR2.L)]`
opr2$sample = row.names(expr_scaled)
opr2 = merge(opr2, infection,
             by = "sample")
opr2$sample
plot1 = ggplot(opr2,
               aes(x = MtOPR2,
                   y = infection_intensity))
plot1 = plot1 +
  geom_point() +
  geom_smooth(method = "lm",
              color = "#40ABAB") +
  theme_linedraw() +
  labs(x = "MtOPR2 Expression)",
       y = "Infection Intensity") +
  theme_tyler()
plot1

ggsave("08_univariate_model_medicago/opr2_expression.png",
       plot1,
       dpi = 600,
       height = 4, width = 4)

irx1 = as.data.frame(expr_scaled[,which(colnames(expr_scaled) %in% MtLBD16.L)])
irx1$MtIRX1 = irx1$`expr_scaled[, which(colnames(expr_scaled) %in% MtLBD16.L)]`
irx1$sample = row.names(expr_scaled)
irx1 = merge(irx1, infection,
             by = "sample")
plot2 = ggplot(irx1,
               aes(x = MtIRX1,
                   y = infection_intensity))
plot2 = plot2 +
  geom_point() +
  geom_smooth(method = "lm",
              color = "springgreen4") +
  theme_linedraw() +
  labs(x = "MtIRX1 (Scaled Expression)",
       y = "Infection Intensity") +
  theme(axis.title = element_text(size = 14,
                                  face = "bold"),
        axis.text = element_text(size = 12))
plot2

