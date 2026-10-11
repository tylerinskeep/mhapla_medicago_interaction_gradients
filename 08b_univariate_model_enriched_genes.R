library(dplyr)
library(ggplot2)
library(ggrepel)

volcano_data <- host_results_mapped %>%
  mutate(
    Highlight = case_when(
      LocusTag == "MtrunA17_Chr1g0192201" ~ "TCTP",
      LocusTag == "MtrunA17_Chr5g0394891" ~ "OPR2",
      TRUE ~ "Other"
    ),
    label = if_else(Highlight == "Other", NA_character_, Highlight)
  ) %>%
  filter(!is.na(S), !is.na(S.padj), S.padj > 0)

ggplot(volcano_data, aes(x = S, y = -log10(S.sig))) +
  geom_point(
    aes(color = Highlight),
    size = 1.5,
    alpha = 0.7
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "grey60"
  ) +
  geom_text_repel(
    aes(label = label,
        color = label),
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.3,
    min.segment.length = 0
  ) +
  scale_color_manual(
    values = c(
      "Other" = "grey80",
      "TCTP" = "#8A4F3D",
      "OPR2" = "#3E5F73"
    )
  ) +
  labs(
    x = "Interaction differential (S)",
    y = expression(-log[10]("P")),
    color = NULL
  ) +
  theme_tyler() +
  theme(legend.position = "none")

ggsave("08_univariate_model_medicago/volcano_opr2_tctp.pdf",
       height = 4, width = 4,
       dpi = 600,
       bg = "transparent")
