#### COEXPRESSION NETWORK ACROSS THE M HAPLA RIL POPULATION ####
## Tyler Inskeep, Groen Lab, UC Riverside ##
# 11 April 2025 #

### SETWD ###
setwd("~/Desktop/PhD/04_projects/04_mhapla_infection_intensity")
source("00_R/plot_utils.R")


### IMPORT DATASETS ###
infection = read.csv("03_infection_intensity/infectionIntensity.csv")
expr = read.csv("04_feature_counts/hapla_log2normCPM.csv")
geneids = expr[,1]
dt = t(expr[,2:97])
colnames(dt) = geneids

### LOAD PACKAGES ###
library(WGCNA)
options(stringsAsFactors = FALSE)

### FILTER LOW VARIANCE GENES ###
gene_var <- apply(dt, 2, var)

# Keep the  most variable genes
top_var_genes <- names(sort(gene_var,
                            decreasing = TRUE))[1:round(length(gene_var)*0.25)]
dt_mod <- dt[,colnames(dt) %in% top_var_genes]
# 6478 genes retained at 75% variance
# 4318 genes retained at to 50%
# 2159 genes retained at 25% variance
geneids = colnames(dt_mod)

### CHOOSING A THRESHOLDING POWER ###
powers <- c(1:20)
sft <- pickSoftThreshold(dt_mod,
                         powerVector = powers,
                         networkType = "unsigned",
                         verbose = 5)

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = "Scale Independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     labels = powers, col = "blue")
abline(h = 0.8, col = "red")

# Selecting a soft thresholding power of 13 - R2 = 0.813 with good mean intramodular connectivity
# @ 25% vartiance picking soft threshold of 6
# Inferring the network
softPower <- 5

net <- blockwiseModules(dt_mod,
                        power = softPower,
                        TOMType = "signed",
                        minModuleSize = 20,
                        mergeCutHeight = 0.25,
                        numericLabels = TRUE,
                        verbose = 3)

moduleColors <- labels2colors(net$colors)
table(moduleColors)  # See how many genes per module

plotDendroAndColors(net$dendrograms[[1]], 
                    moduleColors[net$blockGenes[[1]]],
                    "Module Colors", 
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05)

# Exporting to inspect in cytoscape
# Calculate adjacency
adjacency <- adjacency(dt_mod,
                       power = softPower,
                       type = "signed")  
# Calculate TOM
TOM <- TOMsimilarity(adjacency,
                     TOMType = "signed")  # or "unsigned"
rownames(TOM) <- colnames(TOM) <- geneids  # Ensure proper naming

# Take only turquoise genes
names(moduleColors) = geneids
turquoise_genes <- names(moduleColors[moduleColors == "turquoise"])
blueGenes <- names(moduleColors[moduleColors == "blue"])
brownGenes <- names(moduleColors[moduleColors == "brown"])
# Subset the TOM matrix to only turquoise genes
TOM_turquoise <- TOM[turquoise_genes, turquoise_genes]
TOMblue = TOM[blueGenes, blueGenes]
TOMbrown = TOM[brownGenes, brownGenes]
threshold <- 0.1  # Adjust as needed
TOM_filtered <- TOMbrown
TOM_filtered[TOM_filtered < threshold] <- 0
# Convert to edge list
edge_list <- as.data.frame(as.table(TOM_filtered))
edge_list <- edge_list[edge_list$Freq > 0, ]  # Keep only non-zero edges
colnames(edge_list) <- c("source", "target", "weight")

# Export
write.table(edge_list,
            "coXnetwork/20250507brownModule_edges.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)


## EXTRACTING MODULE EIGENGENES ##
MEs_colored <- moduleEigengenes(dt_mod, colors = moduleColors)$eigengenes
MEs_colored <- orderMEs(MEs_colored)
MEs_colored$sample = row.names(MEs_colored)
write.csv(MEs_colored,
          "coXnetwork/moduleEigen.csv")

## INFERRING MODULE TRAIT RELATIONSHIPS ##
# Normalizing infection intensity
infection = infection[,c(1,9)]
infection$norm_infection = infection$infection_intensity / mean(infection$infection_intensity)

# Scaling the Eigenvectors for Lande-Arnold Regression
ME_scaled = as.data.frame(scale(MEs_colored[,1:9]))
ME_scaled$sample = MEs_colored$sample

# Merging the datasets
dt_merge = merge(ME_scaled, infection,
                 by = "sample")

# Modeling Univariate Association with Module Eigengenes
# Initialize variables to store results
model = data.frame(S = numeric(0),
                   S.error = numeric(0),
                   S.sig = numeric(0))


# Loop through expression_scaled and run the linear model
for(i in dt_merge[,2:10]){ 
  m <- glm(norm_infection ~ i, data=dt_merge,
           family = gaussian) # Running the linear model first
  sum <- summary(m) # Extracting summary of the model
  S = sum$coefficients[2,1]
  S.err = sum$coefficients[2,2]
  S.sig = sum$coefficients[2,4]
  model = rbind(model,
                data.frame(S = S,
                           S.err = S.err,
                           S.sig = S.sig))
}

model$module = colnames(dt_merge[,2:10])
model$padj = p.adjust(model$S.sig,
                      method = "bonferroni")
write.csv(model,
          "coXnetwork/module_infection_correlations.csv")

plot(MEs_colored$MEgreenyellow, infection$norm_infection,
     main = "Module vs. Infection Intensity",
     xlab = "MEgreenyellow", ylab = "Infection intensity")
abline(lm(infection$norm_infection ~ MEs_colored$MEgreenyellow), col = "red")


### Trying a multivariate model ###
model2 = glm(norm_infection ~ MEblue + MEred + MEblack +
               MEgreen + MEbrown + MEyellow +
               MEturquoise + MEpink +
               MEgrey ,
             data=dt_merge,
             family = gaussian) 
sum2 <- summary(model2) # 
coef2 = sum2$coefficients

## COMPUTING CONNECTIVITY SCORES ##
# Rebuild adjacency matrix
adjacency = adjacency(dt_mod,
                      power = softPower,
                      type = "signed")

# Build TOM from adjacency
TOM = TOMsimilarity(adjacency)

# Intramodular connectivity - relationships between transcripts within modules
IMconn <- intramodularConnectivity(TOM, colors = net$colors, scaleByMax = FALSE)
IMconn$geneid = colnames(dt_mod)
IMconn$module = net$colors
IMconn$moduleColor = moduleColors
IMconn$moduleLoyalty = IMconn$kWithin/(IMconn$kWithin + IMconn$kOut)
write.csv(IMconn,
          "coXnetwork/haplaConnectivity.csv")

hub_candidates <- IMconn %>%
  group_by(moduleColor) %>%
  arrange(desc(kWithin), kOut) %>%  # Highest kWithin, lowest kOut
  slice_head(n = 10) %>%             # Top 10 per module
  ungroup()
write.csv(hub_candidates,
          "coXnetwork/haplaHubHighWithinLowWithout.csv")

hub_candidates2 <- IMconn %>%
  group_by(moduleColor) %>%
  arrange(desc(moduleLoyalty)) %>%  # Highest module loyalty
  slice_head(n = 10) %>%             # Top 10 per module
  ungroup()
write.csv(hub_candidates2,
          "coXnetwork/haplaHubModuleLoyalty.csv")


# TRYING TO SPLIT UP SOME MODULES #
# Modules to be split
modules_to_split <- c("turquoise", "blue", "brown")

# Extract TOM matrices for each module
TOM_modules <- lapply(modules_to_split, function(module) {
  module_genes <- names(moduleColors[moduleColors == module])
  TOM[module_genes, module_genes]
})

names(TOM_modules) <- modules_to_split

# Initialize list to store new labels
submodule_labels <- list()

for (module in modules_to_split) {
  
  # Extract the TOM matrix for the module
  TOM_module <- TOM_modules[[module]]
  
  # Calculate dissimilarity matrix
  dissTOM <- 1 - TOM_module
  
  # Hierarchical clustering
  geneTree <- hclust(as.dist(dissTOM), method = "average")
  
  # Plot dendrogram to inspect structure
  pdf(file = paste0(module, "_dendrogram.pdf"), width = 12, height = 8)
  plot(geneTree, main = paste("Dendrogram of", module, "Module"), cex = 0.6)
  dev.off()
  
  # Apply dynamic tree cut
  sub_labels <- cutreeDynamic(
    dendro = geneTree,
    distM = dissTOM,
    deepSplit = 2,  # Adjust as needed
    minClusterSize = 30
  )
  
  # Ensure the labels have gene names
  names(sub_labels) <- colnames(TOM_module)
  
  # Store submodule labels
  submodule_labels[[module]] <- sub_labels
}

# Check structure and distribution of submodules
lapply(submodule_labels, table)
# Initialize new moduleColors vector as a copy of the original
moduleColors_refined <- moduleColors

# Loop through each module to update labels
for (module in names(submodule_labels)) {
  
  # Extract submodule labels
  labels <- submodule_labels[[module]]
  
  # Identify unique labels
  unique_labels <- unique(labels)
  
  # Create new, unique submodule labels
  refined_labels <- labels
  
  for (i in seq_along(unique_labels)) {
    sub_label <- unique_labels[i]
    new_label <- paste0(module, "_", sub_label)  # e.g., "turquoise_1", "turquoise_2"
    
    # Apply new labels
    refined_labels[labels == sub_label] <- new_label
  }
  
  # Update the moduleColors vector
  moduleColors_refined[names(refined_labels)] <- refined_labels
}

# Check the refined color distribution
table(moduleColors_refined)

# Check for missing genes
missing_genes <- setdiff(names(moduleColors), names(moduleColors_refined))
if (length(missing_genes) > 0) {
  cat("Missing genes in refined labels:\n")
  print(missing_genes)
} else {
  cat("No missing genes. All genes accounted for.\n")
}

# Subset expression data to relevant genes
dt_modRefined <- dt_mod[, names(moduleColors_refined)]

# Calculate eigengenes
MEs_refined <- moduleEigengenes(dt_modRefined,
                                colors = moduleColors_refined)$eigengenes

# Inspect the new MEs
head(MEs_refined)

write.csv(MEs_refined,
          "coXnetwork/moduleEigenRefined.csv")

## RECALCULATING MODULE TRAIT RELATIONSHIPS ##
# Normalizing infection intensity
infection = infection[,c(1,9)]
infection$norm_infection = infection$infection_intensity / mean(infection$infection_intensity)

# Scaling the Eigenvectors for Lande-Arnold Regression
ME_scaled = as.data.frame(scale(MEs_refined))
ME_scaled$sample = row.names(ME_scaled)

# Merging the datasets
dt_merge = merge(ME_scaled, infection,
                 by = "sample")

# Modeling Univariate Association with Module Eigengenes
# Initialize variables to store results
model = data.frame(S = numeric(0),
                   S.error = numeric(0),
                   S.sig = numeric(0))


# Loop through expression_scaled and run the linear model
for(i in dt_merge[,2:10]){ 
  m <- glm(norm_infection ~ i, data=dt_merge,
           family = gaussian) # Running the linear model first
  sum <- summary(m) # Extracting summary of the model
  S = sum$coefficients[2,1]
  S.err = sum$coefficients[2,2]
  S.sig = sum$coefficients[2,4]
  model = rbind(model,
                data.frame(S = S,
                           S.err = S.err,
                           S.sig = S.sig))
}

model$module = colnames(dt_merge[,2:10])
model$padj = p.adjust(model$S.sig,
                      method = "bonferroni")
write.csv(model,
          "coXnetwork/module_infection_correlationsRefined.csv")

plot(MEs_refined$MEbrown_1, infection$norm_infection,
     main = "Module vs. Infection Intensity",
     xlab = "MEbrown", ylab = "Infection intensity")
abline(lm(infection$norm_infection ~ MEs_refined$MEbrown_1), col = "red")


colnames(dt_merge)

### Trying a multivariate model ###
model2 = glm(norm_infection ~ MEblue_1 + MEblue_2 + MEbrown_0 +
               MEbrown_1 + MEgrey + MEturquoise_1 +
               MEturquoise_2 + MEyellow ,
             data=dt_merge,
             family = gaussian) 
sum2 <- summary(model2) # 
coef2 = sum2$coefficients

brown2Genes = names(moduleColors_refined[which(moduleColors_refined == "brown_1")])

# Create a data frame with gene names and their assigned modules
refined_modules_df <- data.frame(
  Gene = names(moduleColors_refined),
  Module = moduleColors_refined
)

# Inspect the first few rows
head(refined_modules_df)

write.csv(refined_modules_df,
          "coXnetwork/haplaCoXNetRefined.csv")

# Make a new dendrogram
dissTOM = 1 - TOM
# Recalculate the hierarchical clustering
geneTree_refined <- hclust(as.dist(dissTOM),
                           method = "average")

# Plot the dendrogram
plot(geneTree_refined, main = "Dendrogram of Refined Modules", cex = 0.6, hang = -1)
plotDendroAndColors(geneTree_refined, moduleColors_refined, 
                    "Module Colors", 
                    main = "Refined Dendrogram with Submodules", 
                    cex.dendroLabels = 0.5, 
                    hang = -1)
