############################################################
## Plant GO analysis
## Host: Medicago truncatula
##
## Analyses:
##   1. ORA of significant positive-S genes
##   2. ORA of significant negative-S genes
##   3. Preranked GSEA using S / S.err
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
library(fgsea)


## =========================================================
## 2. Settings
## =========================================================

significance_cutoff <- 0.05

## ORA term-size limits
minimum_go_size_ora <- 5
maximum_go_size_ora <- 500

## GSEA term-size limits
minimum_go_size_gsea <- 10
maximum_go_size_gsea <- 500

host_color <- "#5F8B70"

dir.create(
  "12_GO_medicago",
  showWarnings = FALSE,
  recursive = TRUE
)


## =========================================================
## 3. Import model results
## =========================================================

host_results <- read.csv(
  "08_univariate_model_medicago/medicagoModel.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

host_results = host_results[,2:13]

## Convert transcript labels such as:
## rna-XM_003595123.4
##
## to:
## XM_003595123

host_results <- host_results %>%
  mutate(
    transcript_id = geneid %>%
      str_remove("^rna-") %>%
      str_remove("\\.[0-9]+$")
  )


## =========================================================
## 4. Import NCBI gene2accession
## =========================================================

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
  na.strings = c("-", "NA"),
  stringsAsFactors = FALSE
)

gene2accession_clean <- gene2accession %>%
  filter(
    !is.na(RNA_nucleotide_accession),
    str_detect(
      RNA_nucleotide_accession,
      "^(XM_|NM_|XR_|NR_)"
    )
  ) %>%
  transmute(
    transcript_id = str_remove(
      RNA_nucleotide_accession,
      "\\.[0-9]+$"
    ),
    GeneID = as.character(GeneID),
    Symbol = as.character(Symbol),
    protein_accession,
    genomic_nucleotide_accession
  ) %>%
  distinct()


## =========================================================
## 5. Import NCBI gene_info correctly
## =========================================================

## NCBI gene_info starts:
## tax_id, GeneID, Symbol, LocusTag, ...

gene_info_names <- c(
  "tax_id",
  "GeneID",
  "Symbol",
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

gene_info <- read.delim(
  "gene_info_medicago.txt",
  header = FALSE,
  col.names = gene_info_names,
  quote = "",
  comment.char = "",
  na.strings = c("-", ""),
  stringsAsFactors = FALSE
) %>%
  mutate(
    GeneID = as.character(GeneID)
  )


## =========================================================
## 6. Map transcript IDs to NCBI and Ensembl Rapid IDs
## =========================================================

host_results_mapped <- host_results %>%
  left_join(
    gene2accession_clean,
    by = "transcript_id"
  ) %>%
  left_join(
    gene_info %>%
      dplyr::select(
        GeneID,
        LocusTag,
        dbXrefs,
        description,
        type_of_gene
      ),
    by = "GeneID"
  ) %>%
  mutate(
    ## Extract only the Ensembl Rapid identifier from dbXrefs.
    ##
    ## Example:
    ## "EnsemblRapid:gene238|OtherDB:123"
    ## becomes:
    ## "gene238"
    
    ensembl_gene_id = str_extract(
      dbXrefs,
      "(?<=EnsemblRapid:)gene[[:alnum:]_.-]+"
    ),
    
    S = as.numeric(S),
    S.err = as.numeric(S.err),
    S.padj = as.numeric(S.padj)
  )


## Check how many model genes obtained Ensembl IDs

host_results_mapped %>%
  summarise(
    model_rows = n(),
    rows_with_ncbi_geneid = sum(!is.na(GeneID)),
    rows_with_ensembl_id = sum(!is.na(ensembl_gene_id)),
    fraction_with_ensembl_id = mean(!is.na(ensembl_gene_id))
  ) %>%
  print()


## =========================================================
## 7. Collapse to one result per Ensembl gene
## =========================================================

## Multiple transcript accessions may map to one Ensembl gene.
##
## For GSEA, each gene must occur only once. Here, the row with
## the largest absolute S / S.err statistic is retained.

host_gene_results <- host_results_mapped %>%
  filter(
    !is.na(ensembl_gene_id),
    ensembl_gene_id != "",
    !is.na(S),
    !is.na(S.err),
    S.err > 0,
    !is.na(S.padj)
  ) %>%
  mutate(
    rank_statistic = S / S.err
  ) %>%
  group_by(ensembl_gene_id) %>%
  slice_max(
    order_by = abs(rank_statistic),
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()


## Save the usable gene-level model table

write.csv(
  host_gene_results,
  "08_univariate_model_medicago/medicago_model_ensembl_annotated.csv",
  row.names = FALSE
)


cat(
  "Unique mapped Ensembl genes:",
  n_distinct(host_gene_results$ensembl_gene_id),
  "\n"
)


## =========================================================
## 8. Define tested, positive, and negative gene sets
## =========================================================

host_tested_genes <- host_gene_results %>%
  pull(ensembl_gene_id) %>%
  unique()

host_positive_genes <- host_gene_results %>%
  filter(
    S.padj <= significance_cutoff,
    S > 0
  ) %>%
  pull(ensembl_gene_id) %>%
  unique()

host_negative_genes <- host_gene_results %>%
  filter(
    S.padj <= significance_cutoff,
    S < 0
  ) %>%
  pull(ensembl_gene_id) %>%
  unique()


cat(
  "Tested genes:",
  length(host_tested_genes),
  "\n"
)

cat(
  "Significant positive-S genes:",
  length(host_positive_genes),
  "\n"
)

cat(
  "Significant negative-S genes:",
  length(host_negative_genes),
  "\n"
)


if (
  length(
    intersect(
      host_positive_genes,
      host_negative_genes
    )
  ) > 0
) {
  stop("Positive and negative gene sets overlap.")
}


## =========================================================
## 9. Connect to Ensembl Plants BioMart
## =========================================================

plants_mart <- useMart(
  biomart = "plants_mart",
  host = "https://plants.ensembl.org"
)

medicago_datasets <- searchDatasets(
  mart = plants_mart,
  pattern = "Medicago truncatula"
)

print(medicago_datasets)

if (nrow(medicago_datasets) == 0) {
  stop("No Medicago truncatula BioMart dataset was found.")
}

medicago_dataset_name <- as.character(
  medicago_datasets$dataset[1]
)

medicago_mart <- useMart(
  biomart = "plants_mart",
  dataset = medicago_dataset_name,
  host = "https://plants.ensembl.org"
)


## =========================================================
## 10. Verify GO attributes
## =========================================================

available_attributes <- listAttributes(medicago_mart)

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
  stop(
    "Missing BioMart attributes: ",
    paste(missing_attributes, collapse = ", ")
  )
}


## =========================================================
## 11. Download GO annotations in batches
## =========================================================

get_biomart_go_batches <- function(
    genes,
    mart,
    batch_size = 500
) {
  
  batches <- split(
    genes,
    ceiling(seq_along(genes) / batch_size)
  )
  
  map_dfr(
    batches,
    function(batch_genes) {
      
      message(
        "Querying ",
        length(batch_genes),
        " genes..."
      )
      
      getBM(
        attributes = required_attributes,
        filters = "ensembl_gene_id",
        values = batch_genes,
        mart = mart,
        uniqueRows = TRUE
      )
    }
  )
}


host_go_raw <- get_biomart_go_batches(
  genes = host_tested_genes,
  mart = medicago_mart,
  batch_size = 500
)


## =========================================================
## 12. Standardize GO annotations
## =========================================================

host_go_background <- host_go_raw %>%
  transmute(
    geneid = as.character(ensembl_gene_id),
    goID = as.character(go_id),
    goNAME = as.character(name_1006),
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


## The actual enrichment universe is all tested genes having
## at least one usable GO annotation.

host_annotated_universe <- intersect(
  host_tested_genes,
  unique(host_go_background$geneid)
)

host_positive_annotated <- intersect(
  host_positive_genes,
  host_annotated_universe
)

host_negative_annotated <- intersect(
  host_negative_genes,
  host_annotated_universe
)


cat(
  "Tested genes with GO:",
  length(host_annotated_universe),
  "\n"
)

cat(
  "Positive-S genes with GO:",
  length(host_positive_annotated),
  "\n"
)

cat(
  "Negative-S genes with GO:",
  length(host_negative_annotated),
  "\n"
)


if (length(host_positive_annotated) == 0) {
  stop("No significant positive-S genes have GO annotations.")
}

if (length(host_negative_annotated) == 0) {
  stop("No significant negative-S genes have GO annotations.")
}


## =========================================================
## 13. ORA function
## =========================================================

run_go_ora <- function(
    selected_genes,
    universe_genes,
    go_annotations,
    direction_label,
    min_term_size = 5,
    max_term_size = 500
) {
  
  selected_genes <- intersect(
    unique(as.character(selected_genes)),
    universe_genes
  )
  
  universe_genes <- unique(
    as.character(universe_genes)
  )
  
  go_annotations_use <- go_annotations %>%
    filter(geneid %in% universe_genes) %>%
    distinct(geneid, goID, .keep_all = TRUE)
  
  term_table <- go_annotations_use %>%
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
    filter(
      TermSize >= min_term_size,
      TermSize <= max_term_size
    )
  
  results <- term_table %>%
    rowwise() %>%
    mutate(
      SelectedInTerm = length(
        intersect(
          selected_genes,
          TermGenes
        )
      ),
      
      SelectedOutside =
        length(selected_genes) -
        SelectedInTerm,
      
      BackgroundInTerm =
        TermSize -
        SelectedInTerm,
      
      BackgroundOutside =
        length(universe_genes) -
        length(selected_genes) -
        BackgroundInTerm,
      
      ValidTable =
        SelectedInTerm > 0 &&
        BackgroundOutside >= 0,
      
      FisherResult = list(
        if (ValidTable) {
          fisher.test(
            matrix(
              c(
                SelectedInTerm,
                SelectedOutside,
                BackgroundInTerm,
                BackgroundOutside
              ),
              nrow = 2,
              byrow = TRUE
            ),
            alternative = "greater"
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
      
      FoldEnrichment =
        (
          SelectedInTerm /
            length(selected_genes)
        ) /
        (
          TermSize /
            length(universe_genes)
        ),
      
      GeneRatio =
        SelectedInTerm /
        length(selected_genes),
      
      BackgroundRatio =
        TermSize /
        length(universe_genes),
      
      Direction = direction_label
    ) %>%
    dplyr::select(
      goID,
      goNAME,
      goDOMAIN,
      Direction,
      TermSize,
      SelectedInTerm,
      SelectedOutside,
      BackgroundInTerm,
      BackgroundOutside,
      GeneRatio,
      BackgroundRatio,
      FoldEnrichment,
      OddsRatio,
      PValue,
      FDR
    ) %>%
    arrange(FDR, desc(FoldEnrichment))
  
  results
}


## =========================================================
## 14. Run positive and negative ORA separately
## =========================================================

host_positive_ora <- run_go_ora(
  selected_genes = host_positive_annotated,
  universe_genes = host_annotated_universe,
  go_annotations = host_go_background,
  direction_label = "Positive S",
  min_term_size = minimum_go_size_ora,
  max_term_size = maximum_go_size_ora
)

host_negative_ora <- run_go_ora(
  selected_genes = host_negative_annotated,
  universe_genes = host_annotated_universe,
  go_annotations = host_go_background,
  direction_label = "Negative S",
  min_term_size = minimum_go_size_ora,
  max_term_size = maximum_go_size_ora
)

host_ora_results <- bind_rows(
  host_positive_ora,
  host_negative_ora
)


write.csv(
  host_positive_ora,
  "12_GO_medicago/host_positive_S_ORA.csv",
  row.names = FALSE
)

write.csv(
  host_negative_ora,
  "12_GO_medicago/host_negative_S_ORA.csv",
  row.names = FALSE
)

write.csv(
  host_ora_results,
  "12_GO_medicago/host_positive_negative_ORA_combined.csv",
  row.names = FALSE
)


## Inspect significant ORA terms

host_ora_significant <- host_ora_results %>%
  filter(FDR <= significance_cutoff)

print(
  host_ora_significant,
  n = Inf
)


## =========================================================
## 15. Plot ORA results
## =========================================================

## Select up to 10 terms per direction, after FDR filtering.

host_ora_plot_data <- host_ora_results %>%
  filter(
    FDR <= significance_cutoff,
    is.finite(FoldEnrichment),
    FoldEnrichment > 0
  ) %>%
  group_by(Direction) %>%
  slice_min(
    order_by = FDR,
    n = 10,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  mutate(
    SignedLog10FDR = case_when(
      Direction == "Positive S" ~ -log10(FDR),
      Direction == "Negative S" ~  log10(FDR)
    ),
    
    PlotLabel = paste0(
      goNAME,
      " [",
      goDOMAIN,
      "]"
    )
  )


host_ora_plot <- ggplot(
  host_ora_plot_data,
  aes(
    x = SignedLog10FDR,
    y = fct_reorder(
      goNAME,
      SignedLog10FDR
    )
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey60"
  ) +
  geom_segment(
    aes(
      x = 0,
      xend = SignedLog10FDR,
      yend = goNAME
    ),
    linewidth = 0.6,
    color = "grey75"
  ) +
  geom_point(
    aes(
      size = SelectedInTerm,
    ),
    color = host_color
  ) +
  scale_x_continuous(
    labels = function(x) abs(x)
  ) +
  scale_size_continuous(
    name = "Significant genes",
    range = c(2, 6)
  ) +
  labs(
    x = expression(
      "Directional significance (" *
        -log[10] *
        " FDR)"
    ),
    y = NULL,
    shape = NULL
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(
      color = "black"
    ),
    axis.title.x = element_text(
      face = "bold"
    ),
    legend.position = "bottom"
  )

host_ora_plot


ggsave(
  "12_GO_medicago/host_GO_ORA_positive_negative.png",
  host_ora_plot,
  width = 5,
  height = 4,
  dpi = 600,
  bg = "white"
)

ggsave(
  "12_GO_medicago/host_GO_ORA_positive_negative.pdf",
  host_ora_plot,
  width = 5.5,
  height = 4,
  bg = "white"
)


## =========================================================
## 16. Build GO pathways for preranked GSEA
## =========================================================

## fgsea requires a named list:
## GO term -> vector of gene IDs

go_pathways <- host_go_background %>%
  filter(
    geneid %in% host_annotated_universe
  ) %>%
  distinct(
    goID,
    geneid
  ) %>%
  group_by(goID) %>%
  summarise(
    genes = list(unique(geneid)),
    .groups = "drop"
  ) %>%
  deframe()


## GO metadata for attaching names and domains afterward

go_metadata <- host_go_background %>%
  distinct(
    goID,
    goNAME,
    goDOMAIN
  )


## =========================================================
## 17. Construct the ranked gene statistic
## =========================================================

## Recommended ranking:
## signed model statistic = S / S.err
##
## Positive values mean positive S.
## Negative values mean negative S.

host_gene_ranks <- host_gene_results %>%
  filter(
    ensembl_gene_id %in% host_annotated_universe,
    is.finite(rank_statistic)
  ) %>%
  select(
    ensembl_gene_id,
    rank_statistic
  ) %>%
  distinct(
    ensembl_gene_id,
    .keep_all = TRUE
  )


## Named numeric vector required by fgsea

host_rank_vector <- host_gene_ranks$rank_statistic
names(host_rank_vector) <- host_gene_ranks$ensembl_gene_id

host_rank_vector <- sort(
  host_rank_vector,
  decreasing = TRUE
)


## Break exact ties reproducibly with a negligible offset.
## This does not meaningfully change the rankings.

set.seed(2026)

host_rank_vector <- host_rank_vector +
  runif(
    length(host_rank_vector),
    min = -1e-10,
    max = 1e-10
  )

host_rank_vector <- sort(
  host_rank_vector,
  decreasing = TRUE
)


## =========================================================
## 18. Run preranked GSEA
## =========================================================

set.seed(2026)

host_gsea_raw <- fgsea(
  pathways = go_pathways,
  stats = host_rank_vector,
  minSize = minimum_go_size_gsea,
  maxSize = maximum_go_size_gsea,
  eps = 0
)


host_gsea_results <- as.data.frame(
  host_gsea_raw
) %>%
  rename(
    goID = pathway,
    PValue = pval,
    FDR = padj,
    ES = ES,
    NES = NES,
    TermSize = size
  ) %>%
  left_join(
    go_metadata,
    by = "goID"
  ) %>%
  mutate(
    Direction = case_when(
      NES > 0 ~ "Positive S",
      NES < 0 ~ "Negative S",
      TRUE ~ "Neutral"
    ),
    
    leadingEdge = map_chr(
      leadingEdge,
      ~ paste(.x, collapse = ";")
    )
  ) %>%
  select(
    goID,
    goNAME,
    goDOMAIN,
    Direction,
    TermSize,
    ES,
    NES,
    PValue,
    FDR,
    leadingEdge,
    everything()
  ) %>%
  arrange(FDR)


write.csv(
  host_gsea_results,
  "12_GO_medicago/host_GO_preranked_GSEA.csv",
  row.names = FALSE
)


print(
  host_gsea_results %>%
    filter(FDR <= significance_cutoff),
  n = Inf
)


## =========================================================
## 19. Plot significant GSEA terms
## =========================================================

host_gsea_plot_data <- host_gsea_results %>%
  filter(
    FDR <= significance_cutoff,
    is.finite(NES)
  ) %>%
  group_by(Direction) %>%
  slice_min(
    order_by = FDR,
    n = 10,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  mutate(
    PlotLabel = paste0(
      goNAME,
      " [",
      goDOMAIN,
      "]"
    )
  )


host_gsea_plot <- ggplot(
  host_gsea_plot_data,
  aes(
    x = NES,
    y = fct_reorder(
      PlotLabel,
      NES
    )
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey60"
  ) +
  geom_segment(
    aes(
      x = 0,
      xend = NES,
      yend = PlotLabel
    ),
    linewidth = 0.6,
    color = "grey75"
  ) +
  geom_point(
    aes(
      size = TermSize,
      shape = Direction
    ),
    color = host_color
  ) +
  scale_size_continuous(
    name = "GO term size",
    range = c(2, 6)
  ) +
  labs(
    x = "Normalized enrichment score",
    y = NULL,
    shape = NULL
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(
      color = "black"
    ),
    axis.title.x = element_text(
      face = "bold"
    ),
    legend.position = "bottom"
  )

host_gsea_plot


ggsave(
  "12_GO_medicago/host_GO_preranked_GSEA.png",
  host_gsea_plot,
  width = 6.5,
  height = 5,
  dpi = 600,
  bg = "white"
)

ggsave(
  "12_GO_medicago/host_GO_preranked_GSEA.pdf",
  host_gsea_plot,
  width = 6.5,
  height = 5,
  bg = "white"
)


## =========================================================
## 20. Optional: examine agreement between ORA and GSEA
## =========================================================

host_ora_gsea_comparison <- host_ora_results %>%
  filter(FDR <= 0.10) %>%
  select(
    goID,
    goNAME,
    goDOMAIN,
    ORA_Direction = Direction,
    ORA_FDR = FDR,
    FoldEnrichment
  ) %>%
  full_join(
    host_gsea_results %>%
      filter(FDR <= 0.10) %>%
      select(
        goID,
        GSEA_Direction = Direction,
        GSEA_FDR = FDR,
        NES
      ),
    by = "goID"
  ) %>%
  arrange(
    pmin(
      ORA_FDR,
      GSEA_FDR,
      na.rm = TRUE
    )
  )

write.csv(
  host_ora_gsea_comparison,
  "12_GO_medicago/host_GO_ORA_GSEA_comparison.csv",
  row.names = FALSE
)

print(
  host_ora_gsea_comparison,
  n = 50
)