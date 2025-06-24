### Analysis of CUT&RUN-qPCR
# author: Vera Laub
# last edited: 2025-06-23
# input: .xlsx of Ct values
# output: plots of enrichment over IgG for C&R experiment, Regression analysis for Primertrial

# Load libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(writexl)


##### ANALYSIS H3K27me3 and PBX1 C&R

### PRE-PROCESSING
# Load data
data <- read_excel("~/Documents/postdoc/results/pioneering_activity_project/qPCR/2025-06-20_C&R14+15_Primertrial/2025-06-20 CR13-14 Primertrial.xlsx", 
                   sheet = "results_R")
data <- data %>% filter(`Sample Name` != "H2O")
data <- data %>% filter(`Sample Name` != "ChIP DNA 1:100")
data <- data %>% filter(`Sample Name` != "ChIP DNA 1:250")
data <- data %>% filter(`Sample Name` != "ChIP DNA 1:500")
data <- data %>% filter(`Sample Name` != "IgG 30min E10.5HL")
data <- data %>% filter(`Sample Name` != "H3K9me3 15min E10.5HL")
data <- data %>% filter(`Sample Name` != "H3K9me3 30min E10.5HL")
data <- data %>% filter(`Sample Name` != "H3K9me3 AmpureBeads E10.5HL")

# Set up cell content in variables for analysis
samples <- data$`Sample Name`
targets <- data$`Target Name`
ct_values <- data$CT
ct_mean <- data$`Ct Mean`

# Step 1: Filter rows close to Ct Mean (optional, from before)
filtered_data <- data %>%
  filter(abs(CT - `Ct Mean`) < 1)
filtered_samples <- filtered_data$`Sample Name`
filtered_targets <- filtered_data$`Target Name`
filtered_ct_values <- filtered_data$CT
filtered_ct_mean <- filtered_data$`Ct Mean`

# Step 2: Get IgG CT values per Target (average if needed)
igg_ct <- filtered_data %>%
  filter(`Sample Name` == "IgG 15min E10.5HL") %>%
  group_by(`Target Name`) %>%
  summarise(igg_ct = mean(CT), .groups = "drop")

# Step 3: Join IgG Ct to non-IgG samples
enrichment_data <- filtered_data %>%
  filter(`Sample Name` != "IgG 15min E10.5HL") %>%
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
  labs(title = "Enrichment Over IgG 15min E10.5HL",
       x = "Sample",
       y = "Fold Enrichment (vs IgG)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Then save it
ggsave(
  filename = "2025-06-23_C&R14+15_H3K27me3+PBX1_qPCR_enrichment_plot.png",  # or .pdf
  plot = enrichment_plot,
  path = "~/Documents/postdoc/results/pioneering_activity_project/qPCR/2025-06-20_C&R14+15_Primertrial/",
  width = 10, height = 6, dpi = 300  # adjust as needed
)


##### ANALYSIS H3K9me3 C&R

### PRE-PROCESSING
# Load data
data <- read_excel("~/Documents/postdoc/results/pioneering_activity_project/qPCR/2025-06-20_C&R14+15_Primertrial/2025-06-20 CR13-14 Primertrial.xlsx", 
                   sheet = "results_R")
data <- data %>% filter(`Sample Name` != "H2O")
data <- data %>% filter(`Sample Name` != "ChIP DNA 1:100")
data <- data %>% filter(`Sample Name` != "ChIP DNA 1:250")
data <- data %>% filter(`Sample Name` != "ChIP DNA 1:500")
data <- data %>% filter(`Sample Name` != "IgG 15min E10.5HL")
data <- data %>% filter(`Sample Name` != "H3K27me3 15min E10.5HL")
data <- data %>% filter(`Sample Name` != "H3K27me3 30min E10.5HL")
data <- data %>% filter(`Sample Name` != "H3K27me3 AmpureBeads E10.5HL")
data <- data %>% filter(`Sample Name` != "PBX1 15 min E10.5HL")
data <- data %>% filter(`Sample Name` != "PBX1 30 min E10.5HL")

# Set up cell content in variables for analysis
samples <- data$`Sample Name`
targets <- data$`Target Name`
ct_values <- data$CT
ct_mean <- data$`Ct Mean`

# Step 1: Filter rows close to Ct Mean (optional, from before)
filtered_data <- data %>%
  filter(abs(CT - `Ct Mean`) < 1)
filtered_samples <- filtered_data$`Sample Name`
filtered_targets <- filtered_data$`Target Name`
filtered_ct_values <- filtered_data$CT
filtered_ct_mean <- filtered_data$`Ct Mean`

# Step 2: Get IgG CT values per Target (average if needed)
igg_ct <- filtered_data %>%
  filter(`Sample Name` == "IgG 30min E10.5HL") %>%
  group_by(`Target Name`) %>%
  summarise(igg_ct = mean(CT), .groups = "drop")

# Step 3: Join IgG Ct to non-IgG samples
enrichment_data <- filtered_data %>%
  filter(`Sample Name` != "IgG 30min E10.5HL") %>%
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
  labs(title = "Enrichment Over IgG 30min E10.5HL",
       x = "Sample",
       y = "Fold Enrichment (vs IgG)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Then save it
ggsave(
  filename = "2025-06-23_C&R14+15_H3K9me3_qPCR_enrichment_plot.png",  # or .pdf
  plot = enrichment_plot,
  path = "~/Documents/postdoc/results/pioneering_activity_project/qPCR/2025-06-20_C&R14+15_Primertrial/",
  width = 10, height = 6, dpi = 300  # adjust as needed
)

################################################################################
### Regression Analysis of qPCR to determine primer efficiency
# author: Vera Laub
# last edited: 2025-06-23
# input: .xls of Ct values
# output: Excel file containing summary statistics

# Load libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(writexl)

# Load data
data <- read_excel("~/Documents/postdoc/results/pioneering_activity_project/qPCR/2025-06-20_C&R14+15_Primertrial/2025-06-20 CR13-14 Primertrial.xlsx", 
                   sheet = "results_R")
data <- data %>% filter(`Sample Name` != "H2O")
data <- data %>% filter(`Sample Name` != "IgG 15min E10.5HL")
data <- data %>% filter(`Sample Name` != "H3K27me3 15min E10.5HL")
data <- data %>% filter(`Sample Name` != "H3K27me3 30min E10.5HL")
data <- data %>% filter(`Sample Name` != "H3K27me3 AmpureBeads E10.5HL")
data <- data %>% filter(`Sample Name` != "PBX1 15 min E10.5HL")
data <- data %>% filter(`Sample Name` != "PBX1 30 min E10.5HL")
data <- data %>% filter(`Sample Name` != "IgG 30min E10.5HL")
data <- data %>% filter(`Sample Name` != "H3K9me3 15min E10.5HL")
data <- data %>% filter(`Sample Name` != "H3K9me3 30min E10.5HL")
data <- data %>% filter(`Sample Name` != "H3K9me3 AmpureBeads E10.5HL")


# Set up cell content in variables for analysis
samples <- data$`Sample Name`
targets <- data$`Target Name`
ct_values <- data$CT

# Prepare an empty list to store results for each target
results_list <- list()

# Create an empty data frame to store all summary statistics
summary_df <- data.frame(Target = character(),
                         R_squared = numeric(),
                         Intercept = numeric(),
                         Intercept_p_value = numeric(),
                         Slope = numeric(),
                         Slope_p_value = numeric(),
                         stringsAsFactors = FALSE)

# Loop through each unique target
unique_targets <- unique(targets)
for (target in unique_targets) {
  
  # Filter data for the current target
  target_data <- data %>% filter(`Target Name` == target)
  
  # Extract ct_values and samples for the current target
  ct_values_target <- target_data$CT
  samples_target <- target_data$`Sample Name`
  
  # Perform linear regression
  lm_model <- lm(ct_values_target ~ samples_target)
  
  # Summarize the model
  model_summary <- summary(lm_model)
  
  # Extract the relevant statistics
  r_squared <- model_summary$r.squared
  intercept <- coef(lm_model)[1]
  slope <- coef(lm_model)[2]
  intercept_p_value <- model_summary$coefficients[1, 4]
  slope_p_value <- model_summary$coefficients[2, 4]
  
  # Add this information to the summary data frame
  summary_df <- rbind(summary_df, data.frame(
    Target = target,
    R_squared = r_squared,
    Intercept = intercept,
    Intercept_p_value = intercept_p_value,
    Slope = slope,
    Slope_p_value = slope_p_value
  ))
}

# Sort the data frame by Target
summary_df <- summary_df %>% arrange(Target)

# Define the output Excel file path
output_file <- "~/Documents/postdoc/results/pioneering_activity_project/qPCR/2025-06-20_C&R14+15_Primertrial/regression_summary.xlsx"

# Write the summary data to an Excel file
write_xlsx(summary_df, path = output_file)

# Print the summary dataframe to the console (optional)
print(summary_df)
