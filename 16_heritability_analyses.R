#### Heritability of Intection Intensity Relative to Transcript Abundance ####
## Tyler Inskeep, Groen Lab, UC Riverside ##
# 27 July 2026 #


### SETWD ###
setwd("~/Desktop/PhD/04_projects/04_mhapla_infection_intensity")

source("00_R/plot_utils.R")

# Importing data
host_heritability = read.csv("16_heritability_analyses/host_transcript_heritability_results.csv")
parasite_heritability = read.csv("16_heritability_analyses/parasite_transcript_heritability_results.csv")
parasite_heritability_v2 = read.csv("16_heritability_analyses/parasite_transcript_heritability_80_samples.csv")

gwas = read.delim("16_heritability_analyses/mh10_gwas.assoc.txt")

infection_intensity_h2 <- 0.604

library(tidyverse)
library(patchwork)


### COLORS #####################################################################

host_color <- "#3A6547"
parasite_color <- "#593745"
infection_color <- "#B86740"
gray_color <- "#8A8A8A"

### INFECTION-INTENSITY HERITABILITY ###########################################

# Replace this with the estimate you decide to report.
infection_intensity_h2 <- 0.604

### CLEAN HERITABILITY DATA ####################################################

host_h2 <- host_heritability %>%
  filter(
    converged,
    is.finite(h2),
    h2 >= 0,
    h2 <= 1
  ) %>%
  mutate(species = "Host")

parasite_h2 <- parasite_heritability_v2 %>%
  filter(
    converged,
    is.finite(h2),
    h2 >= 0,
    h2 <= 1
  ) %>%
  mutate(species = "Parasite")

cat(
  "Host transcripts:", nrow(host_h2), "\n",
  "Parasite transcripts:", nrow(parasite_h2), "\n"
)

host_histogram <- ggplot(
  host_h2,
  aes(x = h2)
) +
  geom_histogram(
    binwidth = 0.025,
    boundary = 0,
    fill = host_color,
    color = "white",
    linewidth = 0.2
  ) +
  geom_vline(
    xintercept = infection_intensity_h2,
    color = infection_color,
    linewidth = 1
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    title = "Host transcripts",
    x = NULL,
    y = "Number of transcripts"
  ) +
  theme_tyler()

parasite_histogram <- ggplot(
  parasite_h2,
  aes(x = h2)
) +
  geom_histogram(
    binwidth = 0.025,
    boundary = 0,
    fill = parasite_color,
    color = "white",
    linewidth = 0.2
  ) +
  geom_vline(
    xintercept = infection_intensity_h2,
    color = infection_color,
    linewidth = 1
  ) +
  annotate(
    "text",
    x = infection_intensity_h2,
    y = Inf,
    label = expression(
      italic(h)^2["infection intensity"] == 0.604
    ),
    color = infection_color,
    angle = 360,
    hjust = 0,
    vjust = -0.5,
    size = 3
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    title = "Parasite transcripts",
    x = expression("Transcript heritability (" * italic(h)^2 * ")"),
    y = ""
  ) +
  theme_tyler()

heritability_plot <- host_histogram / parasite_histogram

heritability_plot

parasite_histogram

ggsave(
  "16_heritability_analyses/transcript_heritability_histograms.pdf",
  heritability_plot,
  width = 5,
  height = 4,
  dpi = 600
)

gwas_plot_data <- gwas %>%
  transmute(
    chr = as.character(chr),
    ps = as.numeric(ps),
    p_wald = as.numeric(p_wald)
  ) %>%
  filter(
    !is.na(chr),
    is.finite(ps),
    is.finite(p_wald),
    p_wald > 0,
    p_wald <= 1
  )

contig_order <- gwas_plot_data %>%
  distinct(chr) %>%
  mutate(
    contig_number = as.numeric(
      stringr::str_extract(chr, "(?<=Contig)\\d+")
    )
  ) %>%
  arrange(contig_number, chr) %>%
  pull(chr)

gwas_plot_data <- gwas_plot_data %>%
  mutate(
    chr = factor(chr, levels = contig_order)
  ) %>%
  arrange(chr, ps)

values = c(
  gray = gray_color,
  purple = parasite_color
)


contig_info <- gwas_plot_data %>%
  group_by(chr) %>%
  summarize(
    contig_length = max(ps, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(chr) %>%
  mutate(
    contig_index = row_number(),
    contig_start = lag(
      cumsum(contig_length),
      default = 0
    ),
    contig_group = if_else(
      contig_index %% 2 == 0,
      "purple",
      "gray"
    )
  )

gwas_plot_data <- gwas_plot_data %>%
  select(
    -any_of(c(
      "contig_length",
      "contig_index",
      "contig_start",
      "contig_group",
      "contig_length.x",
      "contig_length.y",
      "contig_index.x",
      "contig_index.y",
      "contig_start.x",
      "contig_start.y",
      "contig_group.x",
      "contig_group.y",
      "cumulative_position",
      "log10_p"
    ))
  ) %>%
  left_join(
    contig_info,
    by = "chr"
  ) %>%
  mutate(
    cumulative_position = ps + contig_start,
    log10_p = -log10(p_wald)
  )

manhattan_plot <- ggplot(
  gwas_plot_data,
  aes(
    x = cumulative_position,
    y = log10_p,
    color = contig_group
  )
) +
  geom_point(
    size = 1.25,
    alpha = 0.8
  ) +
  scale_color_manual(
    values = c(
      gray = "gray",
      purple = "#593745"
    )
  ) +
  scale_x_continuous(
    breaks = NULL,
    expand = expansion(mult = c(0.005, 0.005))
  ) +
  labs(
    x = "M. hapla contigs",
    y = expression(-log[10](italic(P))),
    title = "GWAS of infection intensity"
  ) +
  theme_tyler() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

manhattan_plot

ggsave(
  "16_heritability_analyses/infection_intensity_manhattan.pdf",
  manhattan_plot,
  width = 6,
  height = 2,
  dpi = 600
)


parasite_h2_v2 = parasite_h2 %>% arrange(desc(variance_logCPM))
ggplot(parasite_h2,
       aes(x = variance_logCPM,
           y = h2)) +
  geom_point(alpha = 0.2)

ggplot(parasite_h2,
       aes(x = mean_logCPM,
           y = variance_logCPM)) +
  geom_point(alpha = 0.2)


parasite_h2_v3 = parasite_h2_v2[1:2169,]
ggplot(parasite_h2_v3,
       aes(x = variance_logCPM,
           y = h2)) +
  geom_point(alpha = 0.2)

ggplot(parasite_h2_v3,
       aes(x = mean_logCPM,
           y = variance_logCPM)) +
  geom_point(alpha = 0.2)



heritability = left_join(parasite_heritability,
                         parasite_heritability_v2,
                         by = "transcript")

cor(heritability$h2.x, heritability$h2.y)
