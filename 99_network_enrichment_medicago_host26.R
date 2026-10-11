############################################################
## Cell-type marker enrichment of host module 26
## Formerly the brown module
############################################################

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(readxl)
library(ggplot2)
library(forcats)


## =========================================================
## 1. Read Pereira et al. marker genes
## =========================================================

excel_file <- "PereiraCellRepscRNA.xlsx"

sheets <- excel_sheets(excel_file)[7:30]

marker_df <- map_dfr(
  sheets,
  ~ read_excel(
    excel_file,
    sheet = .x
  )
) %>%
  mutate(
    cluster = as.character(cluster),
    
    ## Match the formatting used in your expression annotation.
    gene_id = str_replace_all(
      as.character(gene_id),
      "_",
      ""
    )
  ) %>%
  filter(
    !is.na(gene_id),
    gene_id != "",
    !is.na(cluster),
    cluster != ""
  ) %>%
  distinct(
    cluster,
    gene_id
  )


## Add cluster names.

cluster_annot <- tibble(
  cluster = as.character(1:24),
  
  ClusterType = c(
    "Cortex-1",
    "Atrichoblast-1",
    "Cortex-2",
    "Atrichoblast-2",
    "Stele-1",
    "Nodule-1",
    "Pericycle",
    "Endodermis",
    "Nodule-2",
    "Suberized Endodermis",
    "Trichoblast",
    "Unidentified-1",
    "Phloem",
    "Dividing Stele",
    "Nodule-3",
    "Infected Stele",
    "Xylem-1",
    "Stele-2",
    "Xylem-2",
    "Unidentified-2",
    "Lateral Root-1",
    "Lateral Root-2",
    "Unidentified-3",
    "Unidentified-4"
  )
)

marker_df <- marker_df %>%
  left_join(
    cluster_annot,
    by = "cluster"
  )

## =========================================================
## 2. Define the tested expression universe
## =========================================================

host_expression_universe <- loadings_df_clean %>%
  transmute(
    gene_id = str_replace_all(
      as.character(geneid),
      "_",
      ""
    )
  ) %>%
  filter(
    !is.na(gene_id),
    gene_id != ""
  ) %>%
  distinct() %>%
  pull(gene_id)


cat(
  "Genes in expression universe:",
  length(host_expression_universe),
  "\n"
)

## =========================================================
## 3. Map host 26 transcripts to Pereira-compatible gene IDs
## =========================================================

module26 <- read.csv(
  "coXnetwork/medicago/ModuleBrown.csv",
  stringsAsFactors = FALSE
)


module26_transcripts <- module26 %>%
  transmute(
    transcript_id = str_remove(
      as.character(x),
      "^rna-"
    )
  )
  ) %>%
  filter(
    !is.na(transcript_id),
    transcript_id != ""
  ) %>%
  distinct()


## Standardize the transcript IDs in the annotation table.
##
## Change NCBIGeneID below only if your transcript-accession
## column in medicagoGeneAnnot has a different name.


module26_mapped <- module26_transcripts %>%
  left_join(
    gene2accession,
    by = c("transcript_id" = "RNA_nucleotide_accession")
  )


module26_mapped = module_brown_mapped

module26_mapping_summary <- module26_mapped %>%
  summarise(
    ModuleTranscripts = n_distinct(transcript_id),
    
    TranscriptsWithGeneID = n_distinct(
      transcript_id[!is.na(gene_id)]
    ),
    
    UniqueGenes = n_distinct(
      gene_id[!is.na(gene_id)]
    )
  )

print(module26_mapping_summary)


module26_genes <- module26_mapped %>%
  filter(
    !is.na(),
    gene_id != ""
  ) %>%
  pull(gene_id) %>%
  unique()


module26_genes_in_universe <- intersect(
  module26_genes,
  host_expression_universe
)


cat(
  "Mapped host 26 genes:",
  length(module26_genes),
  "\n"
)

cat(
  "Host 26 genes in expression universe:",
  length(module26_genes_in_universe),
  "\n"
)

############################################################
## Map expression transcripts and host 26 to Pereira IDs
############################################################

library(dplyr)
library(stringr)


## =========================================================
## 1. Prepare transcript accession -> NCBI GeneID mapping
## =========================================================

transcript_to_ncbi <- gene2accession_clean %>%
  transmute(
    transcript_id = str_remove(
      as.character(transcript_id),
      "\\.[0-9]+$"
    ),
    GeneID = as.character(GeneID)
  ) %>%
  filter(
    !is.na(transcript_id),
    transcript_id != "",
    !is.na(GeneID),
    GeneID != ""
  ) %>%
  distinct(
    transcript_id,
    GeneID
  )


## =========================================================
## 2. Prepare NCBI GeneID -> Pereira-compatible locus ID
## =========================================================

ncbi_to_locus <- geneids %>%
  transmute(
    GeneID = as.character(geneid),
    
    ## Pereira IDs lack the underscore after A17.
    gene_id = str_replace_all(
      as.character(LocusTag),
      "_",
      ""
    )
  ) %>%
  filter(
    !is.na(GeneID),
    GeneID != "",
    !is.na(gene_id),
    gene_id != ""
  ) %>%
  distinct(
    GeneID,
    gene_id
  )

## =========================================================
## 3. Map the expressed-gene background
## =========================================================

expression_transcripts <- tibble(
  transcript_id = counts$geneid
) %>%
  mutate(
    transcript_id = str_remove(
      as.character(transcript_id),
      "^rna-"
    ),
  ) %>%
  filter(
    !is.na(transcript_id),
    transcript_id != ""
  ) %>%
  distinct()


gene2accession$GeneID = as.character(gene2accession$GeneID)

expression_background_mapped <- expression_transcripts %>%
  left_join(
    gene2accession,
    by = c("transcript_id" = "RNA_nucleotide_accession")
  ) %>%
  left_join(
    gene_info_mapping,
    by = c("GeneID" = "GeneSymbol")
  )


expression_mapping_summary <- expression_background_mapped %>%
  summarise(
    ExpressionTranscripts = n_distinct(transcript_id),
    
    TranscriptsWithNCBIGeneID = n_distinct(
      transcript_id[!is.na(GeneID)]
    ),
    
    TranscriptsWithLocusID = n_distinct(
      transcript_id[!is.na(gene_id)]
    ),
    
    UniqueLocusIDs = n_distinct(
      gene_id[!is.na(gene_id)]
    )
  )

print(expression_mapping_summary)


host_expression_universe <- expression_background_mapped %>%
  filter(
    !is.na(LocusTag),
    LocusTag != ""
  ) %>%
  pull(LocusTag) %>%
  unique()


cat(
  "Mapped genes in expression universe:",
  length(host_expression_universe),
  "\n"
)

module26_genes <- module26_mapped %>%
  
  filter(
    
    !is.na(LocusTag),
    
    LocusTag != ""
    
  ) %>%
  
  pull(LocusTag) %>%
  
  unique()

module26_genes_in_universe <- intersect(
  
  module26_genes,
  
  host_expression_universe
  
)

cat(
  
  "Mapped host 26 genes:",
  
  length(module26_genes),
  
  "\n"
  
)

cat(
  
  "Host 26 genes in expression universe:",
  
  length(module26_genes_in_universe),
  
  "\n"
  
)
## =========================================================
## 5. Check overlap with Pereira marker genes
## =========================================================

marker_genes <- marker_df %>%
  filter(
    !is.na(gene_id),
    gene_id != ""
  ) %>%
  pull(gene_id) %>%
  unique()


cat(
  "Marker genes:",
  length(marker_genes),
  "\n"
)

host_expression_universe <- gsub("_", "", host_expression_universe)

module26_genes <- gsub("_", "", module26_genes)

marker_df$gene_id <- gsub("_", "", marker_df$gene_id)
module26_genes_in_universe <- intersect(
  
  module26_genes,
  
  host_expression_universe
  
)

cat(
  "Marker genes represented in expression universe:",
  length(
    intersect(
      marker_genes,
      host_expression_universe
    )
  ),
  "\n"
)
# Marker genes represented in expression universe: 2959 

cat(
  "Host 26 genes that are cell-type markers:",
  length(
    intersect(
      module26_genes_in_universe,
      marker_genes
    )
  ),
  "\n"
)

# Host 26 genes that are cell-type markers: 117 

marker_df %>%
  filter(gene_id %in% host_expression_universe) %>%
  count(ClusterType, sort = TRUE)

marker_sets <- marker_df %>%
  filter(
    gene_id %in% host_expression_universe
  ) %>%
  group_by(
    ClusterType
  ) %>%
  summarise(
    MarkerGenes = list(unique(gene_id)),
    MarkerSetSize = n_distinct(gene_id),
    .groups = "drop"
  ) %>%
  filter(
    MarkerSetSize >= 5
  )

marker_enrichment <- marker_sets %>%
  rowwise() %>%
  mutate(
    
    ModuleMarkers = length(
      intersect(module26_genes_in_universe,
                MarkerGenes)
    ),
    
    ModuleNotMarkers =
      length(module26_genes_in_universe) -
      ModuleMarkers,
    
    BackgroundMarkers = length(
      intersect(
        setdiff(host_expression_universe,
                module26_genes_in_universe),
        MarkerGenes
      )
    ),
    
    BackgroundNotMarkers =
      length(
        setdiff(host_expression_universe,
                module26_genes_in_universe)
      ) -
      BackgroundMarkers,
    
    Fisher = list(
      fisher.test(
        matrix(
          c(
            ModuleMarkers,
            ModuleNotMarkers,
            BackgroundMarkers,
            BackgroundNotMarkers
          ),
          nrow = 2,
          byrow = TRUE
        ),
        alternative = "greater"
      )
    ),
    
    OddsRatio = Fisher$estimate,
    PValue = Fisher$p.value,
    
    FoldEnrichment =
      (ModuleMarkers /
         length(module26_genes_in_universe)) /
      (BackgroundMarkers /
         length(setdiff(host_expression_universe,
                        module26_genes_in_universe)))
    
  ) %>%
  ungroup() %>%
  mutate(
    FDR = p.adjust(PValue, method = "BH")
  ) %>%
  arrange(FDR)

marker_enrichment %>%
  filter(
    FDR <= 0.05,
    ModuleMarkers > 0
  )

marker_enrichment <- marker_enrichment %>%
  rowwise() %>%
  mutate(
    ModuleGenes = paste(
      sort(intersect(module26_genes_in_universe,
                     MarkerGenes)),
      collapse = ";"
    )
  ) %>%
  ungroup()

# Now plotting the stele-2 genes on the network
edges = read.csv("coXnetwork/medicago/brown_edges.csv")

library(dplyr)
library(stringr)

## Build a clean transcript -> locus ID lookup
transcript_locus_lookup <- module26_mapped %>%
  transmute(
    transcript_id = str_remove(
      as.character(transcript_id),
      "\\.[0-9]+$"
    ),
    gene_id = str_remove_all(
      as.character(LocusTag),
      "_"
    )
  ) %>%
  filter(
    !is.na(transcript_id),
    transcript_id != "",
    !is.na(gene_id),
    gene_id != ""
  ) %>%
  distinct(
    transcript_id,
    .keep_all = TRUE
  )

edges_locus <- edges %>%
  mutate(
    source_transcript = str_remove(
      str_remove(
        as.character(source),
        "^rna-"
      ),
      "\\.[0-9]+$"
    ),
    
    target_transcript = str_remove(
      str_remove(
        as.character(target),
        "^rna-"
      ),
      "\\.[0-9]+$"
    )
  ) %>%
  
  left_join(
    transcript_locus_lookup %>%
      rename(
        source_transcript = transcript_id,
        source_gene = gene_id
      ),
    by = "source_transcript"
  ) %>%
  
  left_join(
    transcript_locus_lookup %>%
      rename(
        target_transcript = transcript_id,
        target_gene = gene_id
      ),
    by = "target_transcript"
  )

edges_locus_clean <- edges_locus %>%
  filter(
    !is.na(source_gene),
    !is.na(target_gene)
  ) %>%
  transmute(
    source = source_gene,
    target = target_gene,
    weight = weight
  ) %>%
  filter(
    source != target
  ) %>%
  distinct()

library(igraph)

graph26 <- graph_from_data_frame(
  edges_locus_clean,
  directed = FALSE
)

stele_genes <- marker_df %>%
  filter(
    ClusterType %in% c(
      "Stele-1",
      "Stele-2",
      "Dividing Stele"
    )
  ) %>%
  pull(gene_id) %>%
  str_remove_all("_") %>%
  unique()

V(graph26)$IsStele <-
  V(graph26)$name %in% stele_genes

library(ggraph)
library(ggplot2)

V(graph26)$Degree <- degree(graph26)

cat("Nodes:", vcount(graph26), "\n")
cat("Edges:", ecount(graph26), "\n")
cat("Stele-2 nodes:", sum(V(graph26)$IsStele), "\n")

edge_cutoff <- quantile(
  E(graph26)$weight,
  probs = 0.99,
  na.rm = TRUE
)

edge_cutoff

graph26_plot <- delete_edges(
  graph26,
  E(graph26)[weight < edge_cutoff]
)

cat("Plotting nodes:", vcount(graph26_plot), "\n")
cat("Plotting edges:", ecount(graph26_plot), "\n")
cat("Stele-2 nodes:", sum(V(graph26_plot)$IsStele), "\n")

V(graph26_plot)$PlotDegree <- degree(graph26_plot)

library(ggraph)
library(ggplot2)

set.seed(2026)

layout26 <- create_layout(
  graph26_plot,
  layout = "fr",
  weights = E(graph26_plot)$weight
)

full_degree <- degree(graph26)

V(graph26_plot)$FullDegree <-
  full_degree[V(graph26_plot)$name]

network26_plot <- ggraph(layout26) +
  
  geom_edge_link(
    aes(alpha = weight),
    color = "grey65",
    linewidth = 0.3,
    show.legend = FALSE
  ) +
  
  scale_edge_alpha_continuous(
    range = c(0.08, 0.35)
  ) +
  
  geom_node_point(
    aes(filter = !IsStele),
    color = "grey78",
    alpha = 0.35,
    size = 1
  ) +
  
  geom_node_point(
    aes(
      filter = IsStele,
      size = FullDegree
    ),
    color = "#5F8B70",
    alpha = .7
  ) +
  
  scale_size_continuous(
    name = "Total Network Degree",
    range = c(2.5, 6)
  ) +
  
  coord_equal() +
  theme_void() +
  
  theme(
    legend.position = "bottom"
  )

network26_plot

ggsave(
  filename = paste0(
    "10_coepxression_medicago/",
    "moduleBrown_stele_enrichment.pdf"
  ),
  plot = network26_plot,
  width = 5,
  height = 5,
  bg = "transparent"
)


node_df <- tibble(
  gene = V(graph26)$name,
  Degree = degree(graph26),
  Stele2 = V(graph26)$IsStele
)

node_df %>%
  group_by(Stele2) %>%
  summarise(
    n = n(),
    mean_degree = mean(Degree),
    median_degree = median(Degree),
    sd_degree = sd(Degree)
  )

wilcox.test(
  Degree ~ Stele2,
  data = node_df
)

stele_marker_genes <- marker_df %>%
  filter(
    ClusterType == "Stele-2"
  ) %>%
  mutate(
    gene_id = str_remove_all(gene_id, "_")
  ) %>%
  distinct(gene_id)

module26_stele_genes <- marker_df %>%
  filter(
    ClusterType == "Stele-2"
  ) %>%
  mutate(
    gene_id = str_remove_all(gene_id, "_")
  ) %>%
  filter(
    gene_id %in% module26_genes_in_universe
  ) %>%
  distinct(gene_id)

module26_stele_genes
