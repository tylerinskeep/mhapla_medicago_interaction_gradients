############################################################
## Plant directional GO analysis
## Positive Sh genes versus negative Sh genes
## Host: Medicago truncatula
############################################################

## =========================================================
## 1. Packages
## =========================================================

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)
library(forcats)
library(biomaRt)


## =========================================================
## 2. User settings
## =========================================================

## Your host interaction-differential results.
##
## This object must contain:
##   1. a host gene-ID column
##   2. Sh, or another column containing the signed host
##      interaction differential
##   3. an adjusted P-value column
##
## Change these three names to match your table.

host_results <- read.csv("08_univariate_model_medicago/medicagoModel.csv")
host_results$geneid = gsub("rna-", "", host_results$geneid)
host_results$geneid = stringr::str_sub(host_results$geneid, 1, 12)

# Gene annotations
gene_info_names <- c(
  "tax_id",
  "Symbol",
  "geneid",
  "LocusTag",
  "Synonyms",
  "dbXrefs",
  "chromosome",
  "map_location",
  "description",
  "type_of_gene",
  "Symbol_from_nomenclature_authority",
  "Full_name_from_nomenclature_authority",
  "Nomenclature_status",
  "Other_designations",
  "Modification_date",
  "Feature_type"
)

geneids <- read.delim(
  "gene_info_medicago.txt",
  header = FALSE,
  col.names = gene_info_names,
  quote = "",
  comment.char = "",
  na.strings = "-"
)

gene2accession_names <- c(
  "tax_id",
  "GeneID",
  "status",
  "RNA_nucleotide_accession",
  "RNA_nucleotide_gi",
  "protein_accession",
  "protein_gi",
  "genomic_nucleotide_accession",
  "genomic_nucleotide_gi",
  "start_position",
  "end_position",
  "orientation",
  "assembly",
  "mature_peptide_accession",
  "mature_peptide_gi",
  "Symbol"
)

gene2accession <- read.delim(
  "gene2accession_medicago.txt",
  header = FALSE,
  col.names = gene2accession_names,
  quote = "",
  comment.char = "",
  na.strings = c("-", "NA")
)

gene2accession_clean <- gene2accession %>%
  filter(
    !is.na(RNA_nucleotide_accession),
    str_detect(
      RNA_nucleotide_accession,
      "^XM_|^NM_|^XR_|^NR_"
    )
  ) %>%
  mutate(
    transcript_id = str_remove(
      RNA_nucleotide_accession,
      "\\.[0-9]+$"
    ),
    GeneID = as.character(GeneID)
  ) %>%
  dplyr::select(
    transcript_id,
    GeneID,
    Symbol,
    protein_accession,
    genomic_nucleotide_accession
  ) %>%
  distinct()

host_results_mapped <- host_results %>%
  mutate(
    transcript_id = str_remove(
      as.character(geneid),
      "\\.[0-9]+$"
    )
  ) %>%
  left_join(
    gene2accession_clean,
    by = "transcript_id"
  )

geneids$Symbol = as.character(geneids$Symbol)

host_results_mapped = left_join(host_results_mapped,
                                geneids,
                                by = c("Symbol" = "geneid"))

library(dplyr)

host_results_clean <- host_results_mapped %>%
  transmute(
    transcript_id = transcript_id,
    GeneID,
    Symbol,
    LocusTag = LocusTag,
    description = description,
    type_of_gene = type_of_gene,
    
    S,
    S.err,
    S.sig,
    
    S.quad,
    S.quad.err,
    S.quad.sig,
    
    C,
    C.err,
    C.sig,
    
    S.padj,
    C.padj
  )

write.csv(host_results_clean,
          "08_univariate_model_medicago/medicago_model_annotated.csv")

host_gene_column <- "GeneID"
host_S_column    <- "S"
host_padj_column <- "S.padj"

significance_cutoff <- 0.05
minimum_go_size     <- 5

## Host color used in your figure
host_color <- "#5F8B70"


## =========================================================
## 3. Standardize the host results
## =========================================================

host_results_clean <- host_results_clean %>%
  transmute(
    geneid = as.character(.data[[host_gene_column]]),
    Sh     = as.numeric(.data[[host_S_column]]),
    Sh_padj = as.numeric(.data[[host_padj_column]])
  ) %>%
  filter(
    !is.na(geneid),
    geneid != "",
    !is.na(Sh),
    !is.na(Sh_padj)
  ) %>%
  mutate(
    ## Remove Ensembl version suffixes such as ".1" if present.
    ##
    ## Example:
    ## MtrunA17Chr1g0000101.1
    ## becomes
    ## MtrunA17Chr1g0000101
    geneid = str_remove(geneid, "\\.[0-9]+$")
  ) %>%
  distinct(geneid, .keep_all = TRUE)


## Significant positive and negative host genes

host_positive_genes <- host_results_mapped %>%
  filter(
    S.padj <= significance_cutoff,
    S > 0,
  ) %>%
  pull(dbXrefs) %>%
  str_remove("^EnsemblRapid:") %>%
  unique()

host_negative_genes <- host_results_mapped %>%
  filter(
    S.padj <= significance_cutoff,
    S < 0,
  ) %>%
  pull(dbXrefs) %>%
  str_remove("^EnsemblRapid:") %>%
  unique()

cat(
  "Positive Sh genes:",
  length(host_positive_genes),
  "\n"
)

# 258

cat(
  "Negative Sh genes:",
  length(host_negative_genes),
  "\n"
)

# 264

## Ensure the gene sets do not overlap

overlapping_genes <- intersect(
  host_positive_genes,
  host_negative_genes
)

if (length(overlapping_genes) > 0) {
  stop(
    "Some genes occur in both positive and negative gene sets."
  )
}


## =========================================================
## 4. Connect to Ensembl Plants
## =========================================================

library(biomaRt)

plants_mart <- useMart(
  biomart = "plants_mart",
  host = "https://plants.ensembl.org"
)

## Find the current Medicago truncatula dataset rather than
## assuming that the dataset name has remained unchanged.

medicago_datasets <- searchDatasets(
  mart = plants_mart,
  pattern = "Medicago truncatula"
)

print(medicago_datasets)


if (nrow(medicago_datasets) == 0) {
  stop(
    paste(
      "No Medicago truncatula dataset was found.",
      "Inspect listDatasets(plants_mart) manually."
    )
  )
}


## Usually this should return one Medicago truncatula dataset.
## If it returns multiple datasets, inspect the printed table
## and specify the correct row manually.

medicago_dataset_name <- medicago_datasets$dataset[1]

cat(
  "Using BioMart dataset:",
  medicago_dataset_name,
  "\n"
)



medicago_mart <- useMart(
  biomart = "plants_mart",
  dataset = medicago_dataset_name,
  host = "https://plants.ensembl.org"
)

getBM(
  attributes = "ensembl_gene_id",
  mart = medicago_mart
)

## =========================================================
## 5. Identify current GO attribute names
## =========================================================

available_attributes <- listAttributes(medicago_mart)

## Inspect likely GO attributes

go_attribute_candidates <- available_attributes %>%
  filter(
    str_detect(
      str_to_lower(paste(name, description)),
      "gene stable id|go term|gene ontology|go id|namespace"
    )
  )

print(go_attribute_candidates)


## These are the common Ensembl Plants attribute names.
## The checks below stop with a useful message if Ensembl has
## renamed one of them.

required_attributes <- c(
  "ensembl_gene_id",
  "go_id",
  "name_1006",
  "namespace_1003"
)

missing_attributes <- setdiff(
  required_attributes,
  available_attributes$name
)


if (length(missing_attributes) > 0) {
  
  message(
    "These expected attributes were not found: ",
    paste(missing_attributes, collapse = ", ")
  )
  
  message(
    paste(
      "Examine go_attribute_candidates and replace",
      "required_attributes with the current names."
    )
  )
  
  stop("BioMart GO attribute names need to be updated.")
}


## =========================================================
## 6. Download GO annotations
## =========================================================

## Use all tested host genes as the analysis universe.
## This is preferable to using every annotated gene in the
## genome because only tested genes could have entered either
## significant gene set.

host_tested_genes <- unique(host_results_clean$geneid)

library(dplyr)
library(stringr)
library(biomaRt)

ncbi_attributes <- listAttributes(medicago_mart) %>%
  filter(
    str_detect(
      str_to_lower(paste(name, description)),
      "entrez|ncbi gene"
    )
  )

ncbi_filters <- listFilters(medicago_mart) %>%
  filter(
    str_detect(
      str_to_lower(paste(name, description)),
      "entrez|ncbi gene"
    )
  )

ncbi_attributes
ncbi_filters

host_tested_genes <- host_results_clean %>%
  filter(!is.na(geneid)) %>%
  pull(geneid) %>%
  as.character() %>%
  unique()

host_go_raw <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "entrezgene_id",
    "go_id",
    "name_1006",
    "namespace_1003"
  ),
  filters = "entrezgene_id",
  values = host_tested_genes,
  mart = medicago_mart,
  uniqueRows = TRUE
)

## Query in batches. This is more robust than submitting one
## extremely large request to BioMart.

## Standardize names to match the parasite-side analysis.

host_go_background <- host_go_raw %>%
  transmute(
    geneid   = as.character(ensembl_gene_id),
    goID     = as.character(go_id),
    goNAME   = as.character(name_1006),
    goDOMAIN = as.character(namespace_1003)
  ) %>%
  filter(
    !is.na(geneid),
    geneid != "",
    !is.na(goID),
    goID != "",
    !is.na(goNAME),
    goNAME != "",
    !is.na(goDOMAIN),
    goDOMAIN != ""
  ) %>%
  mutate(
    geneid = str_remove(geneid, "\\.[0-9]+$"),
    
    ## Standardize domain labels in case BioMart returns
    ## spaces rather than underscores.
    goDOMAIN = case_when(
      str_detect(
        str_to_lower(goDOMAIN),
        "biological"
      ) ~ "biological_process",
      
      str_detect(
        str_to_lower(goDOMAIN),
        "molecular"
      ) ~ "molecular_function",
      
      str_detect(
        str_to_lower(goDOMAIN),
        "cellular"
      ) ~ "cellular_component",
      
      TRUE ~ goDOMAIN
    )
  ) %>%
  distinct(
    geneid,
    goID,
    .keep_all = TRUE
  )


cat(
  "Tested host genes:",
  length(host_tested_genes),
  "\n"
)

cat(
  "Tested genes with at least one GO annotation:",
  n_distinct(host_go_background$geneid),
  "\n"
)

cat(
  "Unique GO terms:",
  n_distinct(host_go_background$goID),
  "\n"
)


## Check gene-ID matching explicitly.

host_id_check <- tibble(
  geneid = host_tested_genes
) %>%
  mutate(
    HasGOAnnotation =
      geneid %in% host_go_background$geneid
  )

print(
  host_id_check %>%
    count(HasGOAnnotation)
)


## Optional: inspect unmatched IDs.

unmatched_host_ids <- host_id_check %>%
  filter(!HasGOAnnotation)

head(unmatched_host_ids, 20)


## =========================================================
## 7. Directly compare positive versus negative Sh genes
## =========================================================

compare_positive_negative <- function(
    positive_genes,
    negative_genes,
    go_annotations,
    min_term_size = 5
) {
  
  positive_genes <- intersect(
    unique(as.character(positive_genes)),
    unique(as.character(go_annotations$geneid))
  )
  
  negative_genes <- intersect(
    unique(as.character(negative_genes)),
    unique(as.character(go_annotations$geneid))
  )
  
  message(
    "Positive genes with GO annotations: ",
    length(positive_genes)
  )
  
  message(
    "Negative genes with GO annotations: ",
    length(negative_genes)
  )
  
  if (length(positive_genes) == 0) {
    stop("No positive genes overlap host_go_background$geneid.")
  }
  
  if (length(negative_genes) == 0) {
    stop("No negative genes overlap host_go_background$geneid.")
  }
  
  directional_genes <- union(
    positive_genes,
    negative_genes
  )
  
  directional_annotations <- go_annotations %>%
    mutate(
      geneid = as.character(geneid)
    ) %>%
    filter(geneid %in% directional_genes) %>%
    distinct(
      geneid,
      goID,
      .keep_all = TRUE
    )
  
  term_table <- directional_annotations %>%
    group_by(
      goID,
      goNAME,
      goDOMAIN
    ) %>%
    summarise(
      TermGenes = list(unique(geneid)),
      TermSize = n_distinct(geneid),
      .groups = "drop"
    ) %>%
    filter(TermSize >= min_term_size)
  
  results <- term_table %>%
    rowwise() %>%
    mutate(
      PositiveInTerm = length(
        intersect(positive_genes, TermGenes)
      ),
      
      NegativeInTerm = length(
        intersect(negative_genes, TermGenes)
      ),
      
      PositiveOutside =
        length(positive_genes) - PositiveInTerm,
      
      NegativeOutside =
        length(negative_genes) - NegativeInTerm,
      
      ContingencyTable = list(
        matrix(
          c(
            PositiveInTerm,
            PositiveOutside,
            NegativeInTerm,
            NegativeOutside
          ),
          nrow = 2,
          byrow = TRUE,
          dimnames = list(
            Direction = c("Positive", "Negative"),
            Annotation = c("InTerm", "Outside")
          )
        )
      ),
      
      ValidTable =
        all(rowSums(ContingencyTable) > 0) &&
        all(colSums(ContingencyTable) > 0),
      
      FisherResult = list(
        if (ValidTable) {
          fisher.test(
            ContingencyTable,
            alternative = "two.sided"
          )
        } else {
          NULL
        }
      ),
      
      OddsRatio = if (ValidTable) {
        unname(FisherResult$estimate)
      } else {
        NA_real_
      },
      
      PValue = if (ValidTable) {
        FisherResult$p.value
      } else {
        NA_real_
      }
    ) %>%
    ungroup() %>%
    mutate(
      FDR = p.adjust(
        PValue,
        method = "BH"
      ),
      
      DirectionBias = case_when(
        is.na(OddsRatio) ~ "Not testable",
        OddsRatio < 1 ~ "Negative S",
        OddsRatio > 1 ~ "Positive S",
        TRUE ~ "Neither"
      )
    ) %>%
    dplyr::select(
      goID,
      goNAME,
      goDOMAIN,
      TermSize,
      PositiveInTerm,
      NegativeInTerm,
      PositiveOutside,
      NegativeOutside,
      ValidTable,
      OddsRatio,
      PValue,
      FDR,
      DirectionBias
    ) %>%
    arrange(FDR)
  
  results
}


host_go_direction_comparison <-
  compare_positive_negative(
    positive_genes = host_positive_genes,
    negative_genes = host_negative_genes,
    go_annotations = host_go_background,
    min_term_size = minimum_go_size
  )


print(
  host_go_direction_comparison,
  n = 30
)


## Save the complete results.

write.csv(
  host_go_direction_comparison,
  file = "12_GO_medicago/host_GO_positive_negative_direction_comparison.csv",
  row.names = FALSE
)


## =========================================================
## 8. Check significant host GO results
## =========================================================

host_significant_go <- host_go_direction_comparison %>%
  filter(FDR <= significance_cutoff)

print(
  host_significant_go,
  n = Inf
)


cat(
  "Significant host GO terms:",
  nrow(host_significant_go),
  "\n"
)

print(
  host_significant_go %>%
    count(
      goDOMAIN,
      DirectionBias
    )
)


## =========================================================
## 9. Prepare the panel-D plotting data
## =========================================================

## IMPORTANT:
##
## Fisher's exact-test P-values and FDR values above are based
## on the original count tables.
##
## The 0.5 correction below is used only to make zero and
## infinite odds ratios finite and plottable.

host_direction_plot_data <- host_go_direction_comparison %>%
  filter(FDR <= significance_cutoff) %>%
  mutate(
    CorrectedOddsRatio =
      (
        (PositiveInTerm + 0.5) *
          (NegativeOutside + 0.5)
      ) /
      (
        (PositiveOutside + 0.5) *
          (NegativeInTerm + 0.5)
      ),
    
    EnrichmentScore =
      log2(CorrectedOddsRatio),
    
    SignificantGenes =
      PositiveInTerm +
      NegativeInTerm,
    
    ExclusiveDirection = case_when(
      PositiveInTerm == 0 &
        NegativeInTerm > 0 ~
        "Negative only",
      
      NegativeInTerm == 0 &
        PositiveInTerm > 0 ~
        "Positive only",
      
      TRUE ~ "Both directions"
    )
  )


print(
  host_direction_plot_data %>%
    dplyr::select(
      goID,
      goNAME,
      goDOMAIN,
      PositiveInTerm,
      NegativeInTerm,
      OddsRatio,
      CorrectedOddsRatio,
      EnrichmentScore,
      FDR
    ),
  n = Inf
)


## =========================================================
## 10. Plot the host directional GO comparison
## =========================================================

host_go_plot <- ggplot(
  host_direction_plot_data,
  aes(
    x = EnrichmentScore,
    y = fct_reorder(
      goNAME,
      EnrichmentScore
    )
  )
) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey60",
    linewidth = 0.5
  ) +
  
  geom_segment(
    aes(
      x = 0,
      xend = EnrichmentScore,
      yend = goNAME
    ),
    color = "grey75",
    linewidth = 0.6
  ) +
  
  geom_point(
    aes(
      size = SignificantGenes
    ),
    color = host_color,
    alpha = 0.95
  ) +
  
  scale_size_continuous(
    name = "Significant genes",
    range = c(1.5, 5)
  ) +
  
  labs(
    x = "Enrichment Score",
    y = NULL
  ) +
  
  theme_classic() +
  
  theme(
    axis.text.y = element_text(
      size = 8,
      color = "black"
    ),
    
    axis.text.x = element_text(
      size = 8,
      color = "black"
    ),
    
    axis.title.x = element_text(
      size = 11,
      face = "bold"
    ),
    
    axis.ticks = element_line(
      color = "black"
    ),
    
    legend.position = "bottom",
    
    legend.title = element_text(
      size = 8
    ),
    
    legend.text = element_text(
      size = 8
    )
  )


host_go_plot


ggsave(
  filename = "host_GO_directional_enrichment_panel_D.png",
  plot = host_go_plot,
  width = 5.2,
  height = 4.2,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "host_GO_directional_enrichment_panel_D.pdf",
  plot = host_go_plot,
  width = 5.2,
  height = 4.2,
  bg = "white"
)

## =========================================================
## 11. GO enrichment of the host brown module
## =========================================================

## This analysis asks:
##
## Are particular GO terms overrepresented among genes in the
## brown module relative to all host genes tested in the
## interaction-differential analysis?
##
## host_go_background contains:
##   geneid   = Ensembl Plants / EnsemblRapid gene ID
##   goID
##   goNAME
##   goDOMAIN
##
## Therefore, the module genes must also be converted to
## EnsemblRapid gene IDs before enrichment.


## ---------------------------------------------------------
## 11A. Read and standardize brown-module transcript IDs
## ---------------------------------------------------------

modulebrown <- read.csv(
  "coXnetwork/medicago/ModuleBrown.csv",
  stringsAsFactors = FALSE
)

## Change `x` below if the transcript-ID column has a
## different name.

module_brown_transcripts <- modulebrown %>%
  transmute(
    transcript_id = as.character(x)
  ) %>%
  mutate(
    transcript_id = str_remove(
      transcript_id,
      "^rna-"
    ),

    ## Remove accession-version suffixes rather than
    ## truncating every accession to a fixed length.
    transcript_id = str_remove(
      transcript_id,
      "\\.[0-9]+$"
    )
  ) %>%
  filter(
    !is.na(transcript_id),
    transcript_id != ""
  ) %>%
  distinct()


cat(
  "Unique brown-module transcripts:",
  nrow(module_brown_transcripts),
  "\n"
)


## ---------------------------------------------------------
## 11B. Prepare the NCBI GeneID-to-Ensembl mapping
## ---------------------------------------------------------

## gene2accession_clean was generated earlier and contains:
##
## transcript_id
## GeneID
## Symbol
## protein_accession
## genomic_nucleotide_accession

gene_info_mapping <- geneids %>%
  transmute(
    GeneID = as.character(geneid),
    GeneSymbol = as.character(Symbol),
    LocusTag = as.character(LocusTag),
    dbXrefs = as.character(dbXrefs),
    description = as.character(description),
    type_of_gene = as.character(type_of_gene)
  ) %>%
  distinct()


## ---------------------------------------------------------
## 11C. Map brown-module transcripts to Ensembl gene IDs
## ---------------------------------------------------------

module_brown_mapped <- module_brown_transcripts %>%
  left_join(
    gene2accession_clean,
    by = "transcript_id"
  ) %>%
  left_join(
    gene_info_mapping,
    by = "GeneID"
  ) %>%
  mutate(
    ## dbXrefs can contain multiple database references
    ## separated by "|". Extract only the EnsemblRapid ID.
    ensembl_gene_id = str_extract(
      dbXrefs,
      "(?<=EnsemblRapid:)[^|;]+"
    ),

    ensembl_gene_id = str_remove(
      ensembl_gene_id,
      "\\.[0-9]+$"
    )
  )


write.csv(
  module_brown_mapped,
  file = "coXnetwork/medicago/moduleBrownGenes_mapped.csv",
  row.names = FALSE
)


## Inspect mapping success.

module_brown_mapping_summary <- module_brown_mapped %>%
  summarise(
    ModuleTranscripts = n_distinct(transcript_id),

    TranscriptsWithNCBIGeneID = n_distinct(
      transcript_id[!is.na(GeneID)]
    ),

    TranscriptsWithEnsemblID = n_distinct(
      transcript_id[!is.na(ensembl_gene_id)]
    ),

    UniqueNCBIGenes = n_distinct(
      GeneID[!is.na(GeneID)]
    ),

    UniqueEnsemblGenes = n_distinct(
      ensembl_gene_id[!is.na(ensembl_gene_id)]
    )
  )

print(module_brown_mapping_summary)


## Clean module gene list.

module_brown_genes <- module_brown_mapped %>%
  filter(
    !is.na(ensembl_gene_id),
    ensembl_gene_id != ""
  ) %>%
  pull(ensembl_gene_id) %>%
  unique()


cat(
  "Brown-module genes with Ensembl IDs:",
  length(module_brown_genes),
  "\n"
)


## ---------------------------------------------------------
## 11D. Define the enrichment universe
## ---------------------------------------------------------

## Use tested host genes that have at least one GO annotation.
##
## This ensures that every gene in the Fisher tables was
## capable of being assigned to a GO term.

host_go_universe <- host_go_background %>%
  filter(
    !is.na(geneid),
    geneid != ""
  ) %>%
  pull(geneid) %>%
  unique()


## Retain only brown-module genes occurring in the tested
## and GO-annotated background.

module_brown_genes_in_universe <- intersect(
  module_brown_genes,
  host_go_universe
)


cat(
  "GO-annotated genes in host universe:",
  length(host_go_universe),
  "\n"
)

cat(
  "Brown-module genes in GO universe:",
  length(module_brown_genes_in_universe),
  "\n"
)

cat(
  "Brown-module mapping proportion:",
  round(
    length(module_brown_genes_in_universe) /
      length(module_brown_genes),
    3
  ),
  "\n"
)


if (length(module_brown_genes_in_universe) == 0) {
  stop(
    paste(
      "No brown-module genes overlap the GO background.",
      "Check the EnsemblRapid ID extraction and transcript mapping."
    )
  )
}


## Save the final gene list used for enrichment.

write.table(
  module_brown_genes_in_universe,
  file = paste0(
    "coXnetwork/medicago/",
    "moduleBrown_EnsemblGenes_in_GO_universe.txt"
  ),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)


## ---------------------------------------------------------
## 11E. GO overrepresentation function
## ---------------------------------------------------------

run_module_go_enrichment <- function(
    module_genes,
    go_annotations,
    background_genes = NULL,
    min_term_size = 5
) {

  go_annotations <- go_annotations %>%
    transmute(
      geneid = as.character(geneid),
      goID = as.character(goID),
      goNAME = as.character(goNAME),
      goDOMAIN = as.character(goDOMAIN)
    ) %>%
    filter(
      !is.na(geneid),
      geneid != "",
      !is.na(goID),
      goID != ""
    ) %>%
    distinct(
      geneid,
      goID,
      .keep_all = TRUE
    )

  if (is.null(background_genes)) {
    background_genes <- unique(
      go_annotations$geneid
    )
  }

  background_genes <- unique(
    as.character(background_genes)
  )

  ## Restrict annotations to the specified universe.

  go_annotations <- go_annotations %>%
    filter(
      geneid %in% background_genes
    )

  background_genes <- intersect(
    background_genes,
    unique(go_annotations$geneid)
  )

  module_genes <- intersect(
    unique(as.character(module_genes)),
    background_genes
  )

  nonmodule_genes <- setdiff(
    background_genes,
    module_genes
  )

  message(
    "Background genes: ",
    length(background_genes)
  )

  message(
    "Module genes: ",
    length(module_genes)
  )

  message(
    "Nonmodule background genes: ",
    length(nonmodule_genes)
  )

  if (length(module_genes) == 0) {
    stop(
      "No module genes occur in the GO background."
    )
  }

  if (length(nonmodule_genes) == 0) {
    stop(
      "The module contains the entire background."
    )
  }


  ## Create one list of annotated genes for each GO term.

  go_term_table <- go_annotations %>%
    group_by(
      goID,
      goNAME,
      goDOMAIN
    ) %>%
    summarise(
      TermGenes = list(
        unique(geneid)
      ),

      BackgroundTermSize = n_distinct(
        geneid
      ),

      .groups = "drop"
    ) %>%
    filter(
      BackgroundTermSize >= min_term_size
    )


  ## For each GO term, construct:
  ##
  ##                         In term     Not in term
  ## Brown module               a             b
  ## Remaining background       c             d
  ##
  ## The one-sided Fisher test asks whether the term is
  ## overrepresented in the brown module.

  enrichment_results <- go_term_table %>%
    rowwise() %>%
    mutate(
      ModuleInTerm = length(
        intersect(
          module_genes,
          TermGenes
        )
      ),

      ModuleOutsideTerm =
        length(module_genes) -
        ModuleInTerm,

      BackgroundInTerm = length(
        intersect(
          nonmodule_genes,
          TermGenes
        )
      ),

      BackgroundOutsideTerm =
        length(nonmodule_genes) -
        BackgroundInTerm,

      ContingencyTable = list(
        matrix(
          c(
            ModuleInTerm,
            ModuleOutsideTerm,
            BackgroundInTerm,
            BackgroundOutsideTerm
          ),
          nrow = 2,
          byrow = TRUE,
          dimnames = list(
            GeneSet = c(
              "BrownModule",
              "OtherBackground"
            ),
            Annotation = c(
              "InTerm",
              "OutsideTerm"
            )
          )
        )
      ),

      ValidTable =
        all(rowSums(ContingencyTable) > 0) &&
        all(colSums(ContingencyTable) > 0),

      FisherResult = list(
        if (ValidTable) {
          fisher.test(
            ContingencyTable,
            alternative = "greater"
          )
        } else {
          NULL
        }
      ),

      OddsRatio = if (ValidTable) {
        unname(
          FisherResult$estimate
        )
      } else {
        NA_real_
      },

      PValue = if (ValidTable) {
        FisherResult$p.value
      } else {
        NA_real_
      },

      ## Fraction of the module assigned to this term.
      ModuleGeneRatio =
        ModuleInTerm /
        length(module_genes),

      ## Fraction of the remaining background assigned to it.
      BackgroundGeneRatio =
        BackgroundInTerm /
        length(nonmodule_genes),

      FoldEnrichment = if_else(
        BackgroundGeneRatio > 0,
        ModuleGeneRatio /
          BackgroundGeneRatio,
        Inf
      )
    ) %>%
    ungroup() %>%
    mutate(
      FDR = p.adjust(
        PValue,
        method = "BH"
      )
    ) %>%
    dplyr::select(
      goID,
      goNAME,
      goDOMAIN,
      BackgroundTermSize,
      ModuleInTerm,
      ModuleOutsideTerm,
      BackgroundInTerm,
      BackgroundOutsideTerm,
      ModuleGeneRatio,
      BackgroundGeneRatio,
      FoldEnrichment,
      OddsRatio,
      PValue,
      FDR
    ) %>%
    arrange(
      FDR,
      desc(FoldEnrichment)
    )

  enrichment_results
}


## ---------------------------------------------------------
## 11F. Run brown-module GO enrichment
## ---------------------------------------------------------

module_brown_go_enrichment <- run_module_go_enrichment(
  module_genes =
    module_brown_genes_in_universe,

  go_annotations =
    host_go_background,

  background_genes =
    host_go_universe,

  min_term_size =
    minimum_go_size
)


print(
  module_brown_go_enrichment,
  n = 30
)


write.csv(
  module_brown_go_enrichment,
  file = paste0(
    "12_GO_medicago/",
    "moduleBrown_GO_enrichment_all_terms.csv"
  ),
  row.names = FALSE
)


## ---------------------------------------------------------
## 11G. Significant GO terms
## ---------------------------------------------------------

module_brown_significant_go <-
  module_brown_go_enrichment %>%
  filter(
    FDR <= significance_cutoff,
    ModuleInTerm > 0
  )


cat(
  "Significant brown-module GO terms:",
  nrow(module_brown_significant_go),
  "\n"
)


print(
  module_brown_significant_go,
  n = Inf
)


print(
  module_brown_significant_go %>%
    count(
      goDOMAIN
    )
)


write.csv(
  module_brown_significant_go,
  file = paste0(
    "12_GO_medicago/",
    "moduleBrown_GO_enrichment_significant.csv"
  ),
  row.names = FALSE
)


## ---------------------------------------------------------
## 11H. Add module genes associated with each GO term
## ---------------------------------------------------------

module_brown_go_with_genes <-
  module_brown_significant_go %>%
  left_join(
    host_go_background %>%
      filter(
        geneid %in%
          module_brown_genes_in_universe
      ) %>%
      group_by(
        goID
      ) %>%
      summarise(
        ModuleGenes = paste(
          sort(unique(geneid)),
          collapse = ";"
        ),
        .groups = "drop"
      ),
    by = "goID"
  )


write.csv(
  module_brown_go_with_genes,
  file = paste0(
    "12_GO_medicago/",
    "moduleBrown_GO_enrichment_significant_with_genes.csv"
  ),
  row.names = FALSE
)


## ---------------------------------------------------------
## 11I. Plot significant module-brown GO enrichment
## ---------------------------------------------------------

## Select the strongest terms for plotting so that a large
## number of redundant GO terms does not make the figure
## unreadable.

number_terms_to_plot <- 20

module_brown_go_plot_data <-
  module_brown_significant_go %>%
  arrange(
    FDR,
    desc(FoldEnrichment)
  ) %>%
  slice_head(
    n = number_terms_to_plot
  ) %>%
  mutate(
    NegativeLog10FDR =
      -log10(FDR),

    PlotLabel = paste0(
      goNAME,
      " (",
      goID,
      ")"
    )
  )


module_brown_go_plot <- ggplot(
  module_brown_go_plot_data,
  aes(
    x = FoldEnrichment,
    y = fct_reorder(
      PlotLabel,
      FoldEnrichment
    )
  )
) +

  geom_segment(
    aes(
      x = 1,
      xend = FoldEnrichment,
      yend = PlotLabel
    ),
    color = "grey75",
    linewidth = 0.6
  ) +

  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    color = "grey60",
    linewidth = 0.5
  ) +

  geom_point(
    aes(
      size = ModuleInTerm,
      alpha = NegativeLog10FDR
    ),
    color = host_color
  ) +

  scale_size_continuous(
    name = "Module genes",
    range = c(1.5, 6)
  ) +

  scale_alpha_continuous(
    name = expression(
      -log[10]("FDR")
    ),
    range = c(0.55, 1)
  ) +

  labs(
    x = "Fold enrichment",
    y = NULL
  ) +

  theme_classic() +

  theme(
    axis.text.y = element_text(
      size = 8,
      color = "black"
    ),

    axis.text.x = element_text(
      size = 8,
      color = "black"
    ),

    axis.title.x = element_text(
      size = 11,
      face = "bold"
    ),

    legend.position = "bottom",

    legend.title = element_text(
      size = 8
    ),

    legend.text = element_text(
      size = 8
    )
  )


module_brown_go_plot


ggsave(
  filename = paste0(
    "12_GO_medicago/",
    "moduleBrown_GO_enrichment.png"
  ),
  plot = module_brown_go_plot,
  width = 6.2,
  height = 5.2,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = paste0(
    "12_GO_medicago/",
    "moduleBrown_GO_enrichment.pdf"
  ),
  plot = module_brown_go_plot,
  width = 6.2,
  height = 5.2,
  bg = "white"
)

## =========================================================
## 11. GO enrichment of the host brown module
## =========================================================

## This analysis asks:
##
## Are particular GO terms overrepresented among genes in the
## brown module relative to all host genes tested in the
## interaction-differential analysis?
##
## host_go_background contains:
##   geneid   = Ensembl Plants / EnsemblRapid gene ID
##   goID
##   goNAME
##   goDOMAIN
##
## Therefore, the module genes must also be converted to
## EnsemblRapid gene IDs before enrichment.


## ---------------------------------------------------------
## 11A. Read and standardize brown-module transcript IDs
## ---------------------------------------------------------

modulebrown <- read.csv(
  "coXnetwork/medicago/ModuleBrown.csv",
  stringsAsFactors = FALSE
)

## Change `x` below if the transcript-ID column has a
## different name.

module_brown_transcripts <- modulebrown %>%
  transmute(
    transcript_id = as.character(x)
  ) %>%
  mutate(
    transcript_id = str_remove(
      transcript_id,
      "^rna-"
    ),
    
    ## Remove accession-version suffixes rather than
    ## truncating every accession to a fixed length.
    transcript_id = str_remove(
      transcript_id,
      "\\.[0-9]+$"
    )
  ) %>%
  filter(
    !is.na(transcript_id),
    transcript_id != ""
  ) %>%
  distinct()


cat(
  "Unique brown-module transcripts:",
  nrow(module_brown_transcripts),
  "\n"
)


## ---------------------------------------------------------
## 11B. Prepare the NCBI GeneID-to-Ensembl mapping
## ---------------------------------------------------------

## gene2accession_clean was generated earlier and contains:
##
## transcript_id
## GeneID
## Symbol
## protein_accession
## genomic_nucleotide_accession

gene_info_mapping <- geneids %>%
  transmute(
    GeneID = as.character(geneid),
    GeneSymbol = as.character(Symbol),
    LocusTag = as.character(LocusTag),
    dbXrefs = as.character(dbXrefs),
    description = as.character(description),
    type_of_gene = as.character(type_of_gene)
  ) %>%
  distinct()


## ---------------------------------------------------------
## 11C. Map brown-module transcripts to Ensembl gene IDs
## ---------------------------------------------------------

module_brown_mapped <- module_brown_transcripts %>%
  left_join(
    gene2accession_clean,
    by = "transcript_id"
  ) %>%
  left_join(
    gene_info_mapping,
    by = c("Symbol" = "GeneID")
  ) %>%
  mutate(
    ## dbXrefs can contain multiple database references
    ## separated by "|". Extract only the EnsemblRapid ID.
    ensembl_gene_id = str_extract(
      dbXrefs,
      "(?<=EnsemblRapid:)[^|;]+"
    ),
    
    ensembl_gene_id = str_remove(
      ensembl_gene_id,
      "\\.[0-9]+$"
    )
  )


write.csv(
  module_brown_mapped,
  file = "coXnetwork/medicago/moduleBrownGenes_mapped.csv",
  row.names = FALSE
)


## Inspect mapping success.

module_brown_mapping_summary <- module_brown_mapped %>%
  summarise(
    ModuleTranscripts = n_distinct(transcript_id),
    
    TranscriptsWithNCBIGeneID = n_distinct(
      transcript_id[!is.na(GeneID)]
    ),
    
    TranscriptsWithEnsemblID = n_distinct(
      transcript_id[!is.na(ensembl_gene_id)]
    ),
    
    UniqueNCBIGenes = n_distinct(
      GeneID[!is.na(GeneID)]
    ),
    
    UniqueEnsemblGenes = n_distinct(
      ensembl_gene_id[!is.na(ensembl_gene_id)]
    )
  )

print(module_brown_mapping_summary)


## Clean module gene list.

module_brown_genes <- module_brown_mapped %>%
  filter(
    !is.na(ensembl_gene_id),
    ensembl_gene_id != ""
  ) %>%
  pull(ensembl_gene_id) %>%
  unique()


cat(
  "Brown-module genes with Ensembl IDs:",
  length(module_brown_genes),
  "\n"
)


## ---------------------------------------------------------
## 11D. Define the enrichment universe
## ---------------------------------------------------------

## Use tested host genes that have at least one GO annotation.
##
## This ensures that every gene in the Fisher tables was
## capable of being assigned to a GO term.

host_go_universe <- host_go_background %>%
  filter(
    !is.na(geneid),
    geneid != ""
  ) %>%
  pull(geneid) %>%
  unique()


## Retain only brown-module genes occurring in the tested
## and GO-annotated background.

module_brown_genes_in_universe <- intersect(
  module_brown_genes,
  host_go_universe
)


cat(
  "GO-annotated genes in host universe:",
  length(host_go_universe),
  "\n"
)

cat(
  "Brown-module genes in GO universe:",
  length(module_brown_genes_in_universe),
  "\n"
)

cat(
  "Brown-module mapping proportion:",
  round(
    length(module_brown_genes_in_universe) /
      length(module_brown_genes),
    3
  ),
  "\n"
)


if (length(module_brown_genes_in_universe) == 0) {
  stop(
    paste(
      "No brown-module genes overlap the GO background.",
      "Check the EnsemblRapid ID extraction and transcript mapping."
    )
  )
}


## Save the final gene list used for enrichment.

write.table(
  module_brown_genes_in_universe,
  file = paste0(
    "coXnetwork/medicago/",
    "moduleBrown_EnsemblGenes_in_GO_universe.txt"
  ),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)


## ---------------------------------------------------------
## 11E. GO overrepresentation function
## ---------------------------------------------------------

run_module_go_enrichment <- function(
    module_genes,
    go_annotations,
    background_genes = NULL,
    min_term_size = 5
) {
  
  go_annotations <- go_annotations %>%
    transmute(
      geneid = as.character(geneid),
      goID = as.character(goID),
      goNAME = as.character(goNAME),
      goDOMAIN = as.character(goDOMAIN)
    ) %>%
    filter(
      !is.na(geneid),
      geneid != "",
      !is.na(goID),
      goID != ""
    ) %>%
    distinct(
      geneid,
      goID,
      .keep_all = TRUE
    )
  
  if (is.null(background_genes)) {
    background_genes <- unique(
      go_annotations$geneid
    )
  }
  
  background_genes <- unique(
    as.character(background_genes)
  )
  
  ## Restrict annotations to the specified universe.
  
  go_annotations <- go_annotations %>%
    filter(
      geneid %in% background_genes
    )
  
  background_genes <- intersect(
    background_genes,
    unique(go_annotations$geneid)
  )
  
  module_genes <- intersect(
    unique(as.character(module_genes)),
    background_genes
  )
  
  nonmodule_genes <- setdiff(
    background_genes,
    module_genes
  )
  
  message(
    "Background genes: ",
    length(background_genes)
  )
  
  message(
    "Module genes: ",
    length(module_genes)
  )
  
  message(
    "Nonmodule background genes: ",
    length(nonmodule_genes)
  )
  
  if (length(module_genes) == 0) {
    stop(
      "No module genes occur in the GO background."
    )
  }
  
  if (length(nonmodule_genes) == 0) {
    stop(
      "The module contains the entire background."
    )
  }
  
  
  ## Create one list of annotated genes for each GO term.
  
  go_term_table <- go_annotations %>%
    group_by(
      goID,
      goNAME,
      goDOMAIN
    ) %>%
    summarise(
      TermGenes = list(
        unique(geneid)
      ),
      
      BackgroundTermSize = n_distinct(
        geneid
      ),
      
      .groups = "drop"
    ) %>%
    filter(
      BackgroundTermSize >= min_term_size
    )
  
  
  ## For each GO term, construct:
  ##
  ##                         In term     Not in term
  ## Brown module               a             b
  ## Remaining background       c             d
  ##
  ## The one-sided Fisher test asks whether the term is
  ## overrepresented in the brown module.
  
  enrichment_results <- go_term_table %>%
    rowwise() %>%
    mutate(
      ModuleInTerm = length(
        intersect(
          module_genes,
          TermGenes
        )
      ),
      
      ModuleOutsideTerm =
        length(module_genes) -
        ModuleInTerm,
      
      BackgroundInTerm = length(
        intersect(
          nonmodule_genes,
          TermGenes
        )
      ),
      
      BackgroundOutsideTerm =
        length(nonmodule_genes) -
        BackgroundInTerm,
      
      ContingencyTable = list(
        matrix(
          c(
            ModuleInTerm,
            ModuleOutsideTerm,
            BackgroundInTerm,
            BackgroundOutsideTerm
          ),
          nrow = 2,
          byrow = TRUE,
          dimnames = list(
            GeneSet = c(
              "BrownModule",
              "OtherBackground"
            ),
            Annotation = c(
              "InTerm",
              "OutsideTerm"
            )
          )
        )
      ),
      
      ValidTable =
        all(rowSums(ContingencyTable) > 0) &&
        all(colSums(ContingencyTable) > 0),
      
      FisherResult = list(
        if (ValidTable) {
          fisher.test(
            ContingencyTable,
            alternative = "greater"
          )
        } else {
          NULL
        }
      ),
      
      OddsRatio = if (ValidTable) {
        unname(
          FisherResult$estimate
        )
      } else {
        NA_real_
      },
      
      PValue = if (ValidTable) {
        FisherResult$p.value
      } else {
        NA_real_
      },
      
      ## Fraction of the module assigned to this term.
      ModuleGeneRatio =
        ModuleInTerm /
        length(module_genes),
      
      ## Fraction of the remaining background assigned to it.
      BackgroundGeneRatio =
        BackgroundInTerm /
        length(nonmodule_genes),
      
      FoldEnrichment = if_else(
        BackgroundGeneRatio > 0,
        ModuleGeneRatio /
          BackgroundGeneRatio,
        Inf
      )
    ) %>%
    ungroup() %>%
    mutate(
      FDR = p.adjust(
        PValue,
        method = "BH"
      )
    ) %>%
    dplyr::select(
      goID,
      goNAME,
      goDOMAIN,
      BackgroundTermSize,
      ModuleInTerm,
      ModuleOutsideTerm,
      BackgroundInTerm,
      BackgroundOutsideTerm,
      ModuleGeneRatio,
      BackgroundGeneRatio,
      FoldEnrichment,
      OddsRatio,
      PValue,
      FDR
    ) %>%
    arrange(
      FDR,
      desc(FoldEnrichment)
    )
  
  enrichment_results
}


## ---------------------------------------------------------
## 11F. Run brown-module GO enrichment
## ---------------------------------------------------------

module_brown_go_enrichment <- run_module_go_enrichment(
  module_genes =
    module_brown_genes_in_universe,
  
  go_annotations =
    host_go_background,
  
  background_genes =
    host_go_universe,
  
  min_term_size =
    minimum_go_size
)


print(
  module_brown_go_enrichment,
  n = 30
)


write.csv(
  module_brown_go_enrichment,
  file = paste0(
    "12_GO_medicago/",
    "moduleBrown_GO_enrichment_all_terms.csv"
  ),
  row.names = FALSE
)


## ---------------------------------------------------------
## 11G. Significant GO terms
## ---------------------------------------------------------

module_brown_significant_go <-
  module_brown_go_enrichment %>%
  filter(
    FDR <= significance_cutoff,
    ModuleInTerm > 0
  )


cat(
  "Significant brown-module GO terms:",
  nrow(module_brown_significant_go),
  "\n"
)


print(
  module_brown_significant_go,
  n = Inf
)


print(
  module_brown_significant_go %>%
    count(
      goDOMAIN
    )
)


write.csv(
  module_brown_significant_go,
  file = paste0(
    "12_GO_medicago/",
    "moduleBrown_GO_enrichment_significant.csv"
  ),
  row.names = FALSE
)


## ---------------------------------------------------------
## 11H. Add module genes associated with each GO term
## ---------------------------------------------------------

module_brown_go_with_genes <-
  module_brown_go_enrichment %>%
  left_join(
    host_go_background %>%
      filter(
        geneid %in%
          module_brown_genes_in_universe
      ) %>%
      group_by(
        goID
      ) %>%
      summarise(
        ModuleGenes = paste(
          sort(unique(geneid)),
          collapse = ";"
        ),
        .groups = "drop"
      ),
    by = "goID"
  )


write.csv(
  module_brown_go_with_genes,
  file = paste0(
    "12_GO_medicago/",
    "moduleBrown_GO_enrichment_significant_with_genes.csv"
  ),
  row.names = FALSE
)


## ---------------------------------------------------------
## 11I. Plot significant module-brown GO enrichment
## ---------------------------------------------------------

## Select the strongest terms for plotting so that a large
## number of redundant GO terms does not make the figure
## unreadable.

number_terms_to_plot <- 20

module_brown_go_plot_data <-
  module_brown_go_enrichment %>%
  arrange(
    FDR,
    desc(FoldEnrichment)
  ) %>%
  slice_head(
    n = number_terms_to_plot
  ) %>%
  mutate(
    NegativeLog10FDR =
      -log10(FDR),
    
    PlotLabel = paste0(
      goNAME,
      " (",
      goID,
      ")"
    )
  )


module_brown_go_plot <- ggplot(
  module_brown_go_plot_data,
  aes(
    x = FoldEnrichment,
    y = fct_reorder(
      PlotLabel,
      FoldEnrichment
    )
  )
) +
  
  geom_segment(
    aes(
      x = 1,
      xend = FoldEnrichment,
      yend = PlotLabel
    ),
    color = "grey75",
    linewidth = 0.6
  ) +
  
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    color = "grey60",
    linewidth = 0.5
  ) +
  
  geom_point(
    aes(
      size = ModuleInTerm,
      alpha = NegativeLog10FDR
    ),
    color = host_color
  ) +
  
  scale_size_continuous(
    name = "Module genes",
    range = c(1.5, 6)
  ) +
  
  scale_alpha_continuous(
    name = expression(
      -log[10]("FDR")
    ),
    range = c(0.55, 1)
  ) +
  
  labs(
    x = "Fold enrichment",
    y = NULL
  ) +
  
  theme_classic() +
  
  theme(
    axis.text.y = element_text(
      size = 8,
      color = "black"
    ),
    
    axis.text.x = element_text(
      size = 8,
      color = "black"
    ),
    
    axis.title.x = element_text(
      size = 11,
      face = "bold"
    ),
    
    legend.position = "bottom",
    
    legend.title = element_text(
      size = 8
    ),
    
    legend.text = element_text(
      size = 8
    )
  )


module_brown_go_plot


ggsave(
  filename = paste0(
    "12_GO_medicago/",
    "moduleBrown_GO_enrichment.png"
  ),
  plot = module_brown_go_plot,
  width = 6,
  height = 5,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = paste0(
    "12_GO_medicago/",
    "moduleBrown_GO_enrichment.pdf"
  ),
  plot = module_brown_go_plot,
  width = 6.2,
  height = 5.2,
  bg = "white"
)
