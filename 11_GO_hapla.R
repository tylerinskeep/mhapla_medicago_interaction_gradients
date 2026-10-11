#### GO ENRICHMENT OF POSITIVE AND NEGATIVE S ####
## Tyler Inskeep, Groen Lab, UC Riverside ##

### LOADING PACKAGES ###
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)

### SETWD ###
setwd("~/Desktop/PhD/04_projects/04_mhapla_infection_intensity")

source("00_R/plot_utils.R")

### IMPORT DATA ###
model <- read.csv(
  "selectionModel/hapla_selectionModel.csv",
  stringsAsFactors = FALSE
)

go <- read.delim(
  "mHaplaGO.txt",
  sep = ",",
  stringsAsFactors = FALSE
)

colnames(go) <- c(
  "genomeProj",
  "geneID",
  "geneid",
  "GENE.NAME",
  "GENE.DESC",
  "biotype",
  "goID",
  "goNAME",
  "goDEF",
  "goEVIDENCE",
  "goDOMAIN"
)

model <- model %>%
  mutate(
    geneid = substr(
      geneid,
      1,
      nchar(geneid) - 5
    )
  )

sum(model$geneid %in% go$geneid)
length(unique(model$geneid))

go_clean <- go %>%
  filter(
    !is.na(geneid),
    !is.na(goID),
    goID != ""
  ) %>%
  distinct(
    geneid,
    goID,
    goNAME,
    goDOMAIN
  )

tested_genes <- model %>%
  filter(
    !is.na(S),
    !is.na(S.padj)
  ) %>%
  distinct(geneid)

go_background <- go_clean %>%
  semi_join(
    tested_genes,
    by = "geneid"
  )

background_genes <- unique(go_background$geneid)

length(background_genes) # 3720

positive_genes <- model %>%
  filter(
    S.padj <= 0.05,
    S > 0,
    geneid %in% background_genes
  ) %>%
  distinct(geneid) %>%
  pull(geneid)

negative_genes <- model %>%
  filter(
    S.padj <= 0.05,
    S < 0,
    geneid %in% background_genes
  ) %>%
  distinct(geneid) %>%
  pull(geneid)

length(positive_genes)
length(negative_genes)
length(background_genes)

run_go_enrichment <- function(
    target_genes,
    go_annotations,
    background_genes,
    min_term_size = 5,
    max_term_size = Inf
) {
  
  target_genes <- intersect(
    unique(target_genes),
    background_genes
  )
  
  target_size <- length(target_genes)
  background_size <- length(background_genes)
  
  go_term_info <- go_annotations %>%
    filter(geneid %in% background_genes) %>%
    group_by(
      goID,
      goNAME,
      goDOMAIN
    ) %>%
    summarise(
      term_genes = list(unique(geneid)),
      TermSize = n_distinct(geneid),
      .groups = "drop"
    ) %>%
    filter(
      TermSize >= min_term_size,
      TermSize <= max_term_size
    )
  
  enrichment_results <- go_term_info %>%
    rowwise() %>%
    mutate(
      TargetHits = length(
        intersect(
          target_genes,
          term_genes
        )
      ),
      
      ExpectedHits =
        target_size * TermSize / background_size,
      
      a = TargetHits,
      
      b = target_size - TargetHits,
      
      c = TermSize - TargetHits,
      
      d = background_size - a - b - c,
      
      OddsRatio = ifelse(
        b == 0 | c == 0,
        NA_real_,
        (a * d) / (b * c)
      ),
      
      PValue = fisher.test(
        matrix(
          c(a, b, c, d),
          nrow = 2,
          byrow = TRUE
        ),
        alternative = "greater"
      )$p.value,
      
      TargetGenes = paste(
        sort(
          intersect(
            target_genes,
            term_genes
          )
        ),
        collapse = ";"
      )
    ) %>%
    ungroup() %>%
    mutate(
      FDR = p.adjust(
        PValue,
        method = "BH"
      ),
      
      GeneRatio = TargetHits / target_size,
      
      BackgroundRatio = TermSize / background_size,
      
      FoldEnrichment =
        GeneRatio / BackgroundRatio
    ) %>%
    select(
      goID,
      goNAME,
      goDOMAIN,
      TargetHits,
      TermSize,
      ExpectedHits,
      GeneRatio,
      BackgroundRatio,
      FoldEnrichment,
      OddsRatio,
      PValue,
      FDR,
      TargetGenes
    ) %>%
    arrange(
      FDR,
      desc(FoldEnrichment)
    )
  
  return(enrichment_results)
}

go_positive <- run_go_enrichment(
  target_genes = positive_genes,
  go_annotations = go_background,
  background_genes = background_genes,
  min_term_size = 5
)

go_negative <- run_go_enrichment(
  target_genes = negative_genes,
  go_annotations = go_background,
  background_genes = background_genes,
  min_term_size = 5
)

go_positive %>%
  filter(FDR <= 0.05)

go_negative %>%
  filter(FDR <= 0.05)

positive_BP <- go_positive %>%
  filter(
    goDOMAIN == "biological_process",
    FDR <= 0.05
  )

positive_MF <- go_positive %>%
  filter(
    goDOMAIN == "molecular_function",
    FDR <= 0.05
  )

positive_CC <- go_positive %>%
  filter(
    goDOMAIN == "cellular_component",
    FDR <= 0.05
  )

negative_BP <- go_negative %>%
  filter(
    goDOMAIN == "biological_process",
    FDR <= 0.05
  )

negative_MF <- go_negative %>%
  filter(
    goDOMAIN == "molecular_function",
    FDR <= 0.05
  )

negative_CC <- go_negative %>%
  filter(
    goDOMAIN == "cellular_component",
    FDR <= 0.05
  )

go_enrichment_all <- bind_rows(
  go_positive,
  go_negative
)

write.table(
  go_enrichment_all,
  "11_GO_hapla/positive_negative_S_GO_enrichment.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

go_positive <- go_positive %>%
  mutate(Direction = "Positive S")

go_negative <- go_negative %>%
  mutate(Direction = "Negative S")

go_enrichment_all <- bind_rows(
  go_positive,
  go_negative
)

plot_data <- go_enrichment_all %>%
  filter(
    FDR <= 0.05,
    goDOMAIN == "biological_process"
  ) %>%
  group_by(Direction) %>%
  slice_min(
    order_by = FDR,
    n = 15,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  mutate(
    log10FDR = -log10(FDR),
    Direction = factor(
      Direction,
      levels = c(
        "Negative S",
        "Positive S"
      )
    )
  )

plot_go <- ggplot(
  plot_data,
  aes(
    x = FoldEnrichment,
    y = reorder(goNAME, FoldEnrichment),
    size = TargetHits,
    fill = log10FDR
  )
) +
  geom_point(
    shape = 21
  ) +
  facet_wrap(
    ~ Direction,
    scales = "free_y"
  ) +
  scale_fill_viridis(
    name = expression(-log[10]("FDR"))
  ) +
  labs(
    x = "Fold enrichment",
    y = NULL,
    size = "Gene count"
  ) +
  theme_classic() +
  theme(
    strip.text = element_text(
      face = "bold"
    )
  )

plot_go

ggsave(
  "11_GO_hapla/positive_negative_S_GO_enrichment.png",
  plot_go,
  width = 10,
  height = 7,
  dpi = 300,
  bg = "white"
)

compare_positive_negative <- function(
    positive_genes,
    negative_genes,
    go_annotations,
    min_term_size = 5
) {
  
  directional_genes <- union(
    positive_genes,
    negative_genes
  )
  
  term_info <- go_annotations %>%
    filter(geneid %in% directional_genes) %>%
    group_by(
      goID,
      goNAME,
      goDOMAIN
    ) %>%
    summarise(
      term_genes = list(unique(geneid)),
      TermSize = n_distinct(geneid),
      .groups = "drop"
    ) %>%
    filter(
      TermSize >= min_term_size
    )
  
  results <- term_info %>%
    rowwise() %>%
    mutate(
      PositiveInTerm = length(
        intersect(
          positive_genes,
          term_genes
        )
      ),
      
      NegativeInTerm = length(
        intersect(
          negative_genes,
          term_genes
        )
      ),
      
      PositiveOutside =
        length(positive_genes) -
        PositiveInTerm,
      
      NegativeOutside =
        length(negative_genes) -
        NegativeInTerm,
      
      OddsRatio = fisher.test(
        matrix(
          c(
            PositiveInTerm,
            PositiveOutside,
            NegativeInTerm,
            NegativeOutside
          ),
          nrow = 2,
          byrow = TRUE
        )
      )$estimate,
      
      PValue = fisher.test(
        matrix(
          c(
            PositiveInTerm,
            PositiveOutside,
            NegativeInTerm,
            NegativeOutside
          ),
          nrow = 2,
          byrow = TRUE
        )
      )$p.value
    ) %>%
    ungroup() %>%
    mutate(
      FDR = p.adjust(
        PValue,
        method = "BH"
      ),
      
      DirectionBias = case_when(
        OddsRatio > 1 ~ "Positive S",
        OddsRatio < 1 ~ "Negative S",
        TRUE ~ "Neither"
      )
    ) %>%
    arrange(FDR)
  
  return(results)
}

go_direction_comparison <- compare_positive_negative(
  positive_genes = positive_genes,
  negative_genes = negative_genes,
  go_annotations = go_background,
  min_term_size = 5
)

go_direction_comparison <- compare_positive_negative(
  positive_genes = positive_genes,
  negative_genes = negative_genes,
  go_annotations = go_background,
  min_term_size = 5
)

library(dplyr)
library(ggplot2)
library(forcats)

direction_plot_data <- go_direction_comparison %>%
  filter(
    FDR <= 0.05,
    is.finite(OddsRatio),
    OddsRatio > 0
  ) %>%
  mutate(
    log2_OR = log2(OddsRatio),
    GO_direction = case_when(
      log2_OR < 0 ~ "Negative S",
      log2_OR > 0 ~ "Positive S",
      TRUE ~ "No bias"
    ),
    log10_FDR = -log10(FDR),
    TotalHits = PositiveInTerm + NegativeInTerm
  ) %>%
  arrange(log2_OR)

direction_plot <- ggplot(
  direction_plot_data,
  aes(
    x = log2_OR,
    y = fct_reorder(goNAME, log2_OR)
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_segment(
    aes(
      x = 0,
      xend = log2_OR,
      yend = goNAME
    ),
    linewidth = 0.7
  ) +
  geom_point(
    aes(
      size = TotalHits,
      fill = log10_FDR
    ),
    shape = 21
  ) +
  scale_fill_viridis_c(
    name = expression(-log[10]("FDR"))
  ) +
  labs(
    x = expression(
      log[2]("odds ratio: positive " * S * " versus negative " * S)
    ),
    y = NULL,
    size = "Genes"
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(
      face = "bold",
      size = 12
    ),
    legend.title = element_text(face = "bold")
  )

direction_plot

compare_positive_negative <- function(
    positive_genes,
    negative_genes,
    go_annotations,
    min_term_size = 5
) {
  
  directional_genes <- union(
    positive_genes,
    negative_genes
  )
  
  term_info <- go_annotations %>%
    filter(geneid %in% directional_genes) %>%
    group_by(
      goID,
      goNAME,
      goDOMAIN
    ) %>%
    summarise(
      term_genes = list(unique(geneid)),
      TermSize = n_distinct(geneid),
      .groups = "drop"
    ) %>%
    filter(TermSize >= min_term_size)
  
  results <- term_info %>%
    rowwise() %>%
    mutate(
      PositiveInTerm = length(
        intersect(positive_genes, term_genes)
      ),
      
      NegativeInTerm = length(
        intersect(negative_genes, term_genes)
      ),
      
      PositiveOutside =
        length(positive_genes) - PositiveInTerm,
      
      NegativeOutside =
        length(negative_genes) - NegativeInTerm,
      
      fisher_result = list(
        fisher.test(
          matrix(
            c(
              PositiveInTerm,
              PositiveOutside,
              NegativeInTerm,
              NegativeOutside
            ),
            nrow = 2,
            byrow = TRUE
          ),
          alternative = "two.sided"
        )
      ),
      
      OddsRatio = unname(fisher_result$estimate),
      
      PValue = fisher_result$p.value
    ) %>%
    ungroup() %>%
    select(-fisher_result, -term_genes) %>%
    mutate(
      FDR = p.adjust(PValue, method = "BH"),
      
      DirectionBias = case_when(
        OddsRatio < 1 ~ "Negative S",
        OddsRatio > 1 ~ "Positive S",
        TRUE ~ "Neither"
      )
    ) %>%
    arrange(FDR)
  
  results
}

go_direction_comparison <- compare_positive_negative(
  positive_genes = positive_genes,
  negative_genes = negative_genes,
  go_annotations = go_background,
  min_term_size = 5
)

direction_plot_data <- go_direction_comparison %>%
  mutate(
    CorrectedOR =
      (
        (PositiveInTerm + 0.5) *
          (NegativeOutside + 0.5)
      ) /
      (
        (PositiveOutside + 0.5) *
          (NegativeInTerm + 0.5)
      ),
    
    log2_OR = log2(CorrectedOR),
    
    TotalHits =
      PositiveInTerm + NegativeInTerm,
    
    log10_FDR =
      -log10(pmax(FDR, 1e-300)),
    
    ExclusiveDirection = case_when(
      PositiveInTerm == 0 &
        NegativeInTerm > 0 ~
        "Negative S only",
      
      NegativeInTerm == 0 &
        PositiveInTerm > 0 ~
        "Positive S only",
      
      TRUE ~ "Both directions"
    )
  ) %>%
  filter(FDR <= 0.05)

library(forcats)
library(ggplot2)

direction_plot <- ggplot(
  direction_plot_data,
  aes(
    x = log2_OR,
    y = fct_reorder(goNAME, log2_OR)
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey55"
  ) +
  geom_segment(
    aes(
      x = 0,
      xend = log2_OR,
      yend = goNAME
    ),
    color = "grey75",
    linewidth = 0.7
  ) +
  geom_point(
    aes(size = TotalHits),
    shape = 21,
    fill = "#593745",
    color = "black",
    stroke = 0.3
  ) +
  scale_shape_manual(
    values = c(
      "Negative S only" = 17,
      "Positive S only" = 17,
      "Both directions" = 16
    ),
    name = "GO-term membership"
  ) +
  labs(
    x = "Enrichment Score",
    y = NULL,
    size = "Significant genes"
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(
      face = "bold",
      size = 12
    ),
    legend.position = 'bottom'
  )

direction_plot

ggsave("11_GO_hapla/enrichment_score_GO.pdf",
       height = 4, width = 5,
       dpi = 600, bg = "transparent")
