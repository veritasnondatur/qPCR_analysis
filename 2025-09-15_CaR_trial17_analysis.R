### Analysis of CUT&RUN-qPCR
# author: Vera Laub
# last edited: 2025-09-15
# input: .xlsx of Ct values (set all "undetermined" values manually to 40 cycles)
# output: plots of enrichment over IgG for C&R experiment

# Load libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)

### LOAD DATA
data <- read_excel(
  "~/Documents/postdoc/results/pioneering_activity_project/CUT&RUN/qPCR/2025-09-12_CaR_17/2025-09-12_cr17_Vera_forR.xlsx", 
  sheet = "results_R"
) %>%
  filter(`Sample Name` != "H2O") %>%
  mutate(
    CT = as.character(CT),
    `Ct Mean` = as.character(`Ct Mean`),
    CT = na_if(CT, "Undetermined"),
    `Ct Mean` = na_if(`Ct Mean`, "Undetermined"),
    CT = replace_na(as.numeric(CT), 40),
    `Ct Mean` = replace_na(as.numeric(`Ct Mean`), 40)
  )

# --- FILTER ROWS CLOSE TO Ct Mean ---
filtered_data <- data %>%
  filter(abs(CT - `Ct Mean`) <= 1)

# --- FUNCTION TO COMPUTE ENRICHMENT PER CONDITION ---
compute_enrichment <- function(cond_name, igg_name, df) {
  # Compute mean IgG Ct per target
  igg_ct <- df %>%
    filter(`Sample Name` == igg_name) %>%
    group_by(`Target Name`) %>%
    summarise(IgG_CT = mean(CT), .groups = "drop")
  
  # Join IgG CT to samples of this condition
  enrichment_data <- df %>%
    filter(grepl(cond_name, `Sample Name`), `Sample Name` != igg_name) %>%
    left_join(igg_ct, by = "Target Name") %>%
    mutate(
      delta_ct = CT - IgG_CT,
      fold_enrichment = 2^(IgG_CT - CT)
    )
  
  return(enrichment_data)
}

# --- FUNCTION TO PLOT ENRICHMENT ---
plot_enrichment <- function(enrichment_data, cond_name, output_dir) {
  sample_name <- unique(enrichment_data$`Sample Name`)
  
  # Linear plot
  p_linear <- ggplot(enrichment_data, aes(x = `Sample Name`, y = fold_enrichment, fill = `Target Name`)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    labs(title = paste("Enrichment Over IgG (", cond_name, ")", sep = ""),
         x = "Sample",
         y = "Fold Enrichment (linear)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(
    filename = paste0(output_dir, "/", cond_name, "_enrichment_linear.png"),
    plot = p_linear,
    width = 10, height = 6, dpi = 300
  )
  
  # Log plot
  p_log <- p_linear +
    scale_y_log10() +
    labs(y = "Fold Enrichment (log10 scale)") +
    annotation_logticks(sides = "l")
  
  ggsave(
    filename = paste0(output_dir, "/", cond_name, "_enrichment_log.png"),
    plot = p_log,
    width = 10, height = 6, dpi = 300
  )
  
  return(enrichment_data)
}

# --- OUTPUT DIRECTORY ---
output_dir <- "~/Documents/postdoc/results/pioneering_activity_project/CUT&RUN/qPCR/2025-09-12_CaR_17/plots"
dir.create(output_dir, showWarnings = FALSE)

# --- CONDITION 1 ---
cond1_data <- compute_enrichment("#1", "IgG #1", filtered_data)
cond1_data <- plot_enrichment(cond1_data, "cond1", output_dir)

# --- CONDITION 2 ---
cond2_data <- compute_enrichment("#2", "IgG #2", filtered_data)
cond2_data <- plot_enrichment(cond2_data, "cond2", output_dir)

# --- CONDITION 3 NO LIBRARY ---
cond3_noLib_data <- compute_enrichment("#3 noLib", "IgG #3 noLib", filtered_data)
cond3_noLib_data <- plot_enrichment(cond3_noLib_data, "cond3_noLib", output_dir)

# --- CONDITION 3 LIBRARY ---
cond3_lib_data <- compute_enrichment("#3 lib", "IgG #3 lib", filtered_data)
cond3_lib_data <- plot_enrichment(cond3_lib_data, "cond3_lib", output_dir)

# --- CONDITION 4 NO LIBRARY ---
cond4_noLib_data <- compute_enrichment("#4 noLib 1:10", "IgG #4 noLib 1:10", filtered_data)
cond4_noLib_data <- plot_enrichment(cond4_noLib_data, "cond4_noLib", output_dir)


# --- COMBINED SUMMARY TABLE ---

# Initialize empty list to store data for all conditions
all_enrichment_list <- list()

# Loop over all conditions (same as before)
conditions <- c("#1", "#2", "#3 noLib", "#3 lib", "#4 noLib 1:10")

for (cond_name in conditions) {
  # Identify IgG sample for this condition
  igg_name <- paste0("IgG ", cond_name)
  
  # Compute mean IgG Ct per target
  igg_ct <- filtered_data %>%
    filter(`Sample Name` == igg_name) %>%
    group_by(`Target Name`) %>%
    summarise(igg_ct = mean(CT), .groups = "drop")
  
  # Join IgG Ct to non-IgG samples of this condition
  enrichment_data <- filtered_data %>%
    filter(grepl(cond_name, `Sample Name`), `Sample Name` != igg_name) %>%
    left_join(igg_ct, by = "Target Name") %>%
    mutate(
      condition = cond_name,
      delta_ct = CT - igg_ct,
      fold_enrichment = 2^(igg_ct - CT)
    )
  
  # Append to list
  all_enrichment_list[[cond_name]] <- enrichment_data
}

# Combine all conditions into one table
combined_summary <- bind_rows(all_enrichment_list)

# Reorder columns for clarity
combined_summary <- combined_summary %>%
  select(`Sample Name`, condition, `Target Name`, CT, igg_ct, delta_ct, fold_enrichment, everything())

# --- SAVE AS CSV ---
write.csv(
  combined_summary,
  file = "~/Documents/postdoc/results/pioneering_activity_project/CUT&RUN/qPCR/2025-09-12_CaR_17/combined_enrichment_summary.csv",
  row.names = FALSE
)

# Print message
cat("Combined summary table saved as CSV.\n")
