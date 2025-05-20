### Regression Analysis of CUT&RUN-qPCR
# author: Vera Laub
# last edited: 2025-05-20
# input: .xlsx of Ct values
# output: .xslx file + plots of enrichment over IgG

# Load libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(writexl)

### PRE-PROCESSING
# Load data
data <- read_excel("~/Documents/postdoc/results/pioneering_activity_project/CUT&RUN/qPCR/2025-05-01_C&R_trial12+13_H3K27me3_H3K9me3_PBX1/2025-05-01_CaR_trials12-13_Vera_analysis.xlsx", 
                   sheet = "results_R")
#data <- data %>% filter(`Sample Name` != "H2O")

# Set up cell content in variables for analysis
samples <- data$`Sample Name`
targets <- data$`Target Name`
ct_values <- data$CT
ct_mean <- data$`Ct Mean`


### ANALYSIS
# Step 1: Filter rows close to Ct Mean (optional, from before)
filtered_data <- data %>%
  filter(abs(CT - `Ct Mean`) < 1)
filtered_samples <- filtered_data$`Sample Name`
filtered_targets <- filtered_data$`Target Name`
filtered_ct_values <- filtered_data$CT
filtered_ct_mean <- filtered_data$`Ct Mean`

# Step 2: Get IgG CT values per Target (average if needed)
igg_ct <- filtered_data %>%
  filter(`Sample Name` == "IgG 1:50 S C&R13") %>%
  group_by(`Target Name`) %>%
  summarise(igg_ct = mean(CT), .groups = "drop")

# Step 3: Join IgG Ct to non-IgG samples
enrichment_data <- filtered_data %>%
  filter(`Sample Name` != "IgG 1:50 S C&R13") %>%
  left_join(igg_ct, by = "Target Name") %>%
  mutate(
    delta_ct = CT - igg_ct,
    fold_enrichment = 2^(igg_ct - CT)
  )
   
 
### OUTPUT
# Basic Enrichment Bar Plot
# First, assign your plot to a variable
enrichment_plot <- ggplot(enrichment_data, aes(x = `Sample Name`, y = fold_enrichment, fill = `Target Name`)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  annotate("text", x = Inf, y = 1, label = "FE = 1", vjust = -2, hjust = 1.1, color = "red", size = 3.5) +
  labs(title = "Enrichment Over IgG",
       x = "Sample",
       y = "Fold Enrichment (vs IgG)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Then save it
ggsave(
  filename = "2025-05-01_C&R_qPCR_enrichment_plot.png",  # or .pdf
  plot = enrichment_plot,
  path = "~/Documents/postdoc/results/pioneering_activity_project/CUT&RUN/qPCR/2025-05-01_C&R_trial12+13_H3K27me3_H3K9me3_PBX1/",
  width = 10, height = 6, dpi = 300  # adjust as needed
)
