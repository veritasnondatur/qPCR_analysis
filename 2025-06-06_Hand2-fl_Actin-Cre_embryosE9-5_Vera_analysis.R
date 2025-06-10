### Analysis of RT-qPCR to determine expression levels relative to control
# author: Vera Laub
# last edited: 2025-06-09
# input: .xlsx of Ct values, set Hand2 #9+10 values that were "undetermined" 
  # manually to 35 (maximum amount of amplifications) to enable analysis of all
  # samples
# output: CSV file containing summary statistics and plots

# Load required packages
library(readxl)
library(dplyr)
library(tidyr)

# ---- READ DATA ----
qpcr_data <-read_excel("Documents/postdoc/results/pioneering_activity_project/qPCR/2025-06-06 Hand2-fl Actin-Cre embryosE9-5/2025-06-06 Hand2-fl Actin-Cre embryosE9-5_analysis.xlsx", 
                       sheet = "results_for_R")
qpcr_data <- qpcr_data %>% filter(`CT` != "Undetermined")

# Check the structure
head(qpcr_data)

# Define output folder
output <- "~/Documents/postdoc/results/pioneering_activity_project/qPCR/2025-06-06 Hand2-fl Actin-Cre embryosE9-5/analysis/"

# Convert CT column to numeric, coercing invalid entries to NA
qpcr_data$CT <- as.numeric(qpcr_data$CT)

####################################################
### DELTA CT WITH GAPDH AS REFERENCE GENE
ref_gene <- "Gapdh #3+4"                  # Housekeeping gene

# Calculate mean CT per sample and gene
mean_ct <- qpcr_data %>%
  group_by(`Sample Name`, `Target Name`) %>%
  summarise(mean_CT = mean(CT, na.rm = TRUE), .groups = "drop")

# Pivot to wide format so each gene becomes a column
ct_wide <- mean_ct %>%
  pivot_wider(names_from = `Target Name`, values_from = mean_CT)

# Calculate ΔCT = CT(ref_gene) - CT(gene)
delta_ct <- ct_wide %>%
  mutate(across(
    .cols = -c(`Sample Name`, all_of(ref_gene)),
    .fns = ~ .data[[ref_gene]] - .,
    .names = "delta_{.col}"
  )) %>%
  select(`Sample Name`, starts_with("delta_")) %>%
  pivot_longer(
    cols = -`Sample Name`,
    names_to = "Target",
    values_to = "delta_CT"
  ) %>%
  mutate(`Target Name` = sub("delta_", "", Target)) %>%
  select(-Target)


# ---- EXPORT RESULTS ----
write.csv(delta_ct, file = file.path(output, "qpcr_dctGapdh_results.csv"), row.names = FALSE)

# Load libraries
library(ggplot2)
library(dplyr)

# ---- Preview the delta_ct table ----
head(delta_ct)

# Define your desired gene and sample order
desired_order <- c("Actb #1+2", "Pbx1 #5+6", "Pbx2 #7+8", "Hand2 #9+10")
desired_sample_order <- c("E9.5 embryo #1", "E9.5 embryo #2", "E9.5 embryo #3",
                          "E9.5 embryo #4", "E9.5 embryo #5", "E9.5 embryo #6",
                          "E9.5 embryo #7", "E9.5 embryo #8", "E9.5 embryo #9",
                          "E8.5 LPM", "E10.5 HL")

# Reorder the 'Target Name' factor and sample names (this affects legend and fill color order)
delta_ct$`Target Name` <- factor(delta_ct$`Target Name`, levels = desired_order)
delta_ct$`Sample Name` <- factor(delta_ct$`Sample Name`, levels = desired_sample_order)

# ---- Create the ΔCT barplot ----
ggplot(delta_ct, aes(x = `Target Name`, y = delta_CT, fill = `Sample Name`)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  labs(
    title = expression("ΔCT per Gene (normalized to Gapdh)"),
    x = "Target Gene",
    y = expression(Delta * "CT value")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  ) +
  scale_fill_brewer(palette = "Set3")

# Save the plot
ggsave(file.path(output, "qpcr_deltaCTGapdh.png"), width = 10, height = 6, dpi = 300)

################################################################################
### DELTA DELTA CT WITH Gapdh AS REFERENCE GENE AND E8.5 LPM AS REFERENCE SAMPLE
# Define reference sample again
ref_sample <- "E8.5 LPM"

# Pivot delta_ct to wide format so each sample becomes a column
delta_ct_wide <- delta_ct %>%
  pivot_wider(names_from = `Sample Name`, values_from = delta_CT)

# Extract the reference delta CT values (as a vector)
# Step 1: Extract delta CT values of the reference sample (E8.5 LPM)
ref_values <- delta_ct_wide %>%
  select(`Target Name`, ref_delta_CT = `E8.5 LPM`)

# Step 2: Join reference delta CTs back to the wide table
delta_delta_ct <- delta_ct_wide %>%
  left_join(ref_values, by = "Target Name")

# Step 3: Reshape long and compute ΔΔCT and fold change
delta_delta_ct_long <- delta_delta_ct %>%
  pivot_longer(
    cols = -c(`Target Name`, ref_delta_CT),
    names_to = "Sample Name",
    values_to = "delta_CT"
  ) %>%
  mutate(
    delta_delta_CT = delta_CT - ref_delta_CT,
    fold_change = 2^(-delta_delta_CT)
  )

# Step 4: Preview and save
head(delta_delta_ct_long)

write.csv(delta_delta_ct_long, file = file.path(output, "qpcr_ddctGapdh_results.csv"), row.names = FALSE)


## Make plot
# Reorder factors to control the order in the plot
delta_delta_ct_long$`Target Name` <- factor(delta_delta_ct_long$`Target Name`,
                                            levels = c("Actb #1+2", "Pbx1 #5+6", "Pbx2 #7+8", "Hand2 #9+10"))

delta_delta_ct_long$`Sample Name` <- factor(delta_delta_ct_long$`Sample Name`,
                                            levels = c("E8.5 LPM", "E9.5 embryo #1", "E9.5 embryo #2", "E9.5 embryo #3",
                                                       "E9.5 embryo #4", "E9.5 embryo #5", "E9.5 embryo #6",
                                                       "E9.5 embryo #7", "E9.5 embryo #8", "E9.5 embryo #9", "E10.5 HL"))

# Plot fold change (log2 scale recommended for visual clarity)
ggplot(delta_delta_ct_long, aes(x = `Target Name`, y = fold_change, fill = `Sample Name`)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black") +
  scale_y_continuous(trans = "log2",  # log2 scale
                     breaks = c(0.25, 0.5, 1, 2, 4, 8),
                     labels = c("0.25", "0.5", "1", "2", "4", "8")) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  labs(
    title = expression("Fold Change (2"^-ΔΔ*"CT) relative to E8.5 LPM"),
    x = "Target Gene",
    y = "Fold Change (log2 scale)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  ) +
  scale_fill_brewer(palette = "Set3")

# Save the plot
ggsave(file.path(output, "qpcr_fold_change_log2_anti-Gapdh+E8.5LPM.png"), width = 10, height = 6, dpi = 300)


####################################################
### DELTA CT WITH ACTB AS REFERENCE GENE
ref_gene <- "Actb #1+2"

# Calculate mean CT per sample and gene
mean_ct <- qpcr_data %>%
  group_by(`Sample Name`, `Target Name`) %>%
  summarise(mean_CT = mean(CT, na.rm = TRUE), .groups = "drop")

# Pivot to wide format so each gene becomes a column
ct_wide <- mean_ct %>%
  pivot_wider(names_from = `Target Name`, values_from = mean_CT)

# Calculate ΔCT = CT(ref_gene) - CT(gene)
delta_ct <- ct_wide %>%
  mutate(across(
    .cols = -c(`Sample Name`, all_of(ref_gene)),
    .fns = ~ .data[[ref_gene]] - .,
    .names = "delta_{.col}"
  )) %>%
  select(`Sample Name`, starts_with("delta_")) %>%
  pivot_longer(
    cols = -`Sample Name`,
    names_to = "Target",
    values_to = "delta_CT"
  ) %>%
  mutate(`Target Name` = sub("delta_", "", Target)) %>%
  select(-Target)

# ---- EXPORT ΔCT RESULTS ----
write.csv(delta_ct, file = file.path(output, "qpcr_dctActb_results.csv"), row.names = FALSE)

# ---- Set factor levels for plotting ----
desired_order <- c("Gapdh #3+4", "Pbx1 #5+6", "Pbx2 #7+8", "Hand2 #9+10")
desired_sample_order <- c("E9.5 embryo #1", "E9.5 embryo #2", "E9.5 embryo #3",
                          "E9.5 embryo #4", "E9.5 embryo #5", "E9.5 embryo #6",
                          "E9.5 embryo #7", "E9.5 embryo #8", "E9.5 embryo #9",
                          "E8.5 LPM", "E10.5 HL")

delta_ct$`Target Name` <- factor(delta_ct$`Target Name`, levels = desired_order)
delta_ct$`Sample Name` <- factor(delta_ct$`Sample Name`, levels = desired_sample_order)

# ---- Plot ΔCT ----
ggplot(delta_ct, aes(x = `Target Name`, y = delta_CT, fill = `Sample Name`)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  labs(
    title = expression("ΔCT per Gene (normalized to Actb)"),
    x = "Target Gene",
    y = expression(Delta * "CT value")
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank()) +
  scale_fill_brewer(palette = "Set3")

ggsave(file.path(output, "qpcr_deltaCT_Actb.png"), width = 10, height = 6, dpi = 300)

################################################################################
### ΔΔCT WITH Actb AS REFERENCE GENE AND E8.5 LPM AS REFERENCE SAMPLE
ref_sample <- "E8.5 LPM"

# Pivot delta_ct to wide format so each sample becomes a column
delta_ct_wide <- delta_ct %>%
  pivot_wider(names_from = `Sample Name`, values_from = delta_CT)

# Extract delta CTs of reference sample (E8.5 LPM)
ref_values <- delta_ct_wide %>%
  select(`Target Name`, ref_delta_CT = `E8.5 LPM`)

# Join reference values back to full data
delta_delta_ct <- delta_ct_wide %>%
  left_join(ref_values, by = "Target Name")

# Reshape and compute ΔΔCT and fold change
delta_delta_ct_long <- delta_delta_ct %>%
  pivot_longer(
    cols = -c(`Target Name`, ref_delta_CT),
    names_to = "Sample Name",
    values_to = "delta_CT"
  ) %>%
  mutate(
    delta_delta_CT = delta_CT - ref_delta_CT,
    fold_change = 2^(-delta_delta_CT)
  )

# Export
write.csv(delta_delta_ct_long, file = file.path(output, "qpcr_ddctActb_results.csv"), row.names = FALSE)

# ---- Plot Fold Change ----
delta_delta_ct_long$`Target Name` <- factor(delta_delta_ct_long$`Target Name`, levels = desired_order)
delta_delta_ct_long$`Sample Name` <- factor(delta_delta_ct_long$`Sample Name`, levels = desired_sample_order)

ggplot(delta_delta_ct_long, aes(x = `Target Name`, y = fold_change, fill = `Sample Name`)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black") +
  scale_y_continuous(trans = "log2", 
                     breaks = c(0.25, 0.5, 1, 2, 4, 8),
                     labels = c("0.25", "0.5", "1", "2", "4", "8")) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  labs(
    title = expression("Fold Change (2"^-ΔΔ*"CT) relative to E8.5 LPM (Actb-normalized)"),
    x = "Target Gene",
    y = "Fold Change (log2 scale)"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank()) +
  scale_fill_brewer(palette = "Set3")

ggsave(file.path(output, "qpcr_fold_change_log2_anti-Actb+E8.5LPM.png"), width = 10, height = 6, dpi = 300)
