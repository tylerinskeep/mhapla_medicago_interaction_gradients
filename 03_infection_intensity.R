#### Estimating Infection Intensity in Medicago truncatula ####
## Tyler Inskeep, Groen Lab, UC Riverside ##
# 10 April 2025 #

library(forcats)
library(purrr)
library(stringr)

### SETWD ###
setwd("~/Desktop/PhD/04_projects/04_mhapla_infection_intensity")

source("00_R/plot_utils.R")

# Testing importing one data
file <- read.delim(
  
  "01_alignment_summary/GC55_idxstats.tsv",
  
  header = FALSE,
  
  col.names = c("contig", "length", "mapped", "unmapped")
  
)

# Annotate contigs by genome origin
# Separate unmapped count
unmapped_count <- file %>%
  filter(origin == "unmapped") %>%
  pull(unmapped)

if (length(unmapped_count) == 0) unmapped_count <- 0  # failsafe

# Summarize mapped reads (excluding "*")
summary <- file %>%
  filter(origin != "unmapped") %>%
  group_by(origin) %>%
  summarize(mapped_reads = sum(mapped), .groups = "drop")

# Reshape and calculate percentages
summary_wide <- summary %>%
  pivot_wider(names_from = origin, values_from = mapped_reads, values_fill = 0) %>%
  mutate(
    unmapped = unmapped_count,
    total_mapped = nematode + host,
    total_reads = total_mapped + unmapped,
    host_mapping_rate = host / total_reads * 100,
    nematode_mapping_rate = nematode / total_reads * 100,
    infection_intensity = nematode / total_mapped * 100,
    sample = "GC55"
  )

library(tidyverse)

# Path to idxstats files
stats_path <- "01_alignment_summary"

# List all *_idxstats.tsv files
files <- list.files(stats_path, pattern = "_idxstats.tsv$", full.names = TRUE)

# Function to process each sample file
process_idxstats <- function(file_path) {
  sample_name <- str_remove(basename(file_path), "_idxstats.tsv")
  
  df <- read.delim(file_path,
                   header=FALSE,
                 col.names = c("contig", "length", "mapped", "unmapped")
                 )
  
  # Annotate origin
  df <- df %>%
    mutate(origin = case_when(
      str_starts(contig, "Mh") ~ "nematode",
      str_starts(contig, "NC") | str_starts(contig, "NW") ~ "host",
      contig == "*" ~ "unmapped",
      TRUE ~ "unknown"
    ))
  
  # Extract unmapped count from the "*" row
  unmapped_count <- df %>%
    filter(origin == "unmapped") %>%
    pull(unmapped)
  if (length(unmapped_count) == 0) unmapped_count <- 0
  
  # Sum mapped reads per origin (excluding unmapped)
  summary <- df %>%
    filter(origin != "unmapped") %>%
    group_by(origin) %>%
    summarize(mapped_reads = sum(mapped), .groups = "drop")
  
  # Pull reads by name
  host_reads     <- summary %>% filter(origin == "host") %>% pull(mapped_reads)
  nematode_reads <- summary %>% filter(origin == "nematode") %>% pull(mapped_reads)
  
  if (length(host_reads) == 0) host_reads <- 0
  if (length(nematode_reads) == 0) nematode_reads <- 0
  
  total_mapped   <- host_reads + nematode_reads
  total_reads    <- total_mapped + unmapped_count
  
  host_mapping_rate     <- host_reads / total_reads * 100
  nematode_mapping_rate <- nematode_reads / total_reads * 100
  infection_intensity   <- nematode_reads / total_mapped * 100
  
  # Return summary row
  tibble(
    sample = sample_name,
    host = host_reads,
    nematode = nematode_reads,
    unmapped = unmapped_count,
    total_mapped = total_mapped,
    total_reads = total_reads,
    host_mapping_rate = host_mapping_rate,
    nematode_mapping_rate = nematode_mapping_rate,
    infection_intensity = infection_intensity
  )
}

# Loop through all samples
all_summaries <- map_dfr(files, process_idxstats)

# Save to CSV
write_csv(all_summaries,
          "align_summary/mapping_summary.csv")

print("✅ All mapping summaries written to mapping_summary.csv")

infectionIntensity = read.csv("03_infection_intensity/infectionIntensity.csv")

# Sort by infection_intensity
infectionIntensity <- infectionIntensity %>%
  arrange(desc(infection_intensity)) %>%
  mutate(sample = fct_inorder(sample))

# Basic barplot
plot1 = ggplot(infectionIntensity, aes(y = sample,
                          x = infection_intensity)) +
  geom_col(fill = "#533943") +
  labs(
    y = "Inbred Line",
    x = "Infection Intensity"
  ) +
  theme_tyler() +
  theme(axis.text.y = element_blank())
plot1

ggsave("03_infection_intensity/infection_intensity.pdf",
       plot1,
       dpi = 600,
       height = 4, width = 4
       )

infectionIntensity$total_reads = as.numeric(infectionIntensity$total_reads)
infectionIntensity$infection_intensity = as.numeric(infectionIntensity$infection_intensity)

model = cor.test(infectionIntensity$infection_intensity,
                 infectionIntensity$total_reads,
                 method = "spearman")
model

#Spearman's rank correlation rho
#data:  infectionIntensity$infection_intensity and infectionIntensity$total_reads
#S = 172036, p-value = 0.1042
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#       rho 
#-0.1668204 

infectionIntensity$total_reads = infectionIntensity$total_reads  / 1000000
infectionIntensity$host = infectionIntensity$host  / 1000000
infectionIntensity$nematode = infectionIntensity$nematode  / 1000000

plot2 = ggplot(infectionIntensity,
               aes(x = total_reads,
                   y = infection_intensity)) +
  geom_point() +
  geom_smooth(method = "lm",
              color = "#B86740") +
  labs(x = "Library Size (Million Reads)",
       y = "Infection Intensity"
  ) +
  theme_tyler()
  
plot2

model = cor.test(infectionIntensity$nematode,
                 infectionIntensity$total_reads,
                 method = "spearman")
model

ggsave("03_infection_intensity/infection_intensity_vs_depth.pdf",
       plot2,
       dpi = 600,
       height = 4, width = 4
)

plot2 = ggplot(infectionIntensity,
               aes(x = total_reads,
                   y = host)) +
  geom_point() +
  geom_smooth(method = "lm",
              color = "#3A6547") +
  labs(x = "Library Size (Million Reads)",
       y = "Host reads (Million)"
  ) +
  theme_tyler()

plot2

ggsave("03_infection_intensity/medicago_vs_depth.pdf",
       plot2,
       dpi = 600,
       height = 4, width = 4
)

plot2 = ggplot(infectionIntensity,
               aes(x = total_reads,
                   y = nematode)) +
  geom_point() +
  geom_smooth(method = "lm",
              color = "#593745") +
  labs(x = "",
       y = "Parasite reads (Million)"
  ) +
  theme_tyler() 

plot2

ggsave("03_infection_intensity/nematode_vs_depth.pdf",
       plot2,
       dpi = 600,
       height = 4, width = 4
)
