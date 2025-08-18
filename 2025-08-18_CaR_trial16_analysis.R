### Analysis of CUT&RUN-qPCR
# author: Vera Laub
# last edited: 2025-08-18
# input: .xlsx of Ct values (set all "undetermined" values manually to 40 cycles)
# output: plots of enrichment over IgG for C&R experiment

# Load libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(writexl)

### PRE-PROCESSING
# Load data
data <- read_excel("~/Documents/postdoc/results/pioneering_activity_project/CUT&RUN/qPCR/2025-08-15_CutAndRun16_Vera/2025-08-15_CutAndRun16_Vera_forR.xlsx", 
                   sheet = "results_R") %>%
  filter(`Sample Name` != "H2O") %>%
  filter(`Sample Name` != "ChIP_Input_MAX_07/19_1:100")

# Optional: filter rows close to Ct mean
filtered_data <- data %>%
  filter(abs(CT - `Ct Mean`) < 1)

### FUNCTION TO RUN ANALYSIS PER CONDITION
analyze_condition <- function(cond_name, data) {
  # Get IgG name for this condition
  igg_name <- paste0("IgG_", cond_name)
  
  # Compute mean IgG Ct per target
  igg_ct <- data %>%
    filter(`Sample Name` == igg_name) %>%
    group_by(`Target Name`) %>%
    summarise(igg_ct = mean(CT), .groups = "drop")
  
  # Join IgG Ct to non-IgG samples of this condition
  enrichment_data <- data %>%
    filter(grepl(cond_name, `Sample Name`), `Sample Name` != igg_name) %>%
    left_join(igg_ct, by = "Target Name") %>%
    mutate(
      delta_ct = CT - igg_ct,
      fold_enrichment = 2^(igg_ct - CT)
    )
  
  # Plot
  enrichment_plot <- ggplot(enrichment_data, aes(x = `Sample Name`, y = fold_enrichment, fill = `Target Name`)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    annotate("text", x = Inf, y = 1, label = "FE = 1", vjust = -2, hjust = 1.1, color = "red", size = 3.5) +
    labs(title = paste("Enrichment Over IgG (", cond_name, ")", sep = ""),
         x = "Sample",
         y = "Fold Enrichment (vs IgG)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save plot
  ggsave(
    filename = paste0("2025-08-18_C&R16_qPCR_", cond_name, "_enrichment_plot.png"),
    plot = enrichment_plot,
    path = "~/Documents/postdoc/results/pioneering_activity_project/CUT&RUN/qPCR/2025-08-15_CutAndRun16_Vera/",
    width = 10, height = 6, dpi = 300
  )
}

### RUN ANALYSIS FOR ALL CONDITIONS
conditions <- c("cond1", "cond2", "cond3", "cond3_library", "cond3_NOlibrary")

for (cond in conditions) {
  analyze_condition(cond, filtered_data)
}


### Code per condition, including linear and log scale plots

# --- CONDITION 1 ONLY ---

# Define condition name and IgG sample name
cond_name <- "cond1"
igg_name  <- "IgG_cond1"

# Get mean IgG Ct per target
igg_ct <- filtered_data %>%
  filter(`Sample Name` == igg_name) %>%
  group_by(`Target Name`) %>%
  summarise(igg_ct = mean(CT), .groups = "drop")

# Join IgG Ct to condition 1 samples (excluding IgG itself)
enrichment_data <- filtered_data %>%
  filter(grepl(cond_name, `Sample Name`), `Sample Name` != igg_name) %>%
  left_join(igg_ct, by = "Target Name") %>%
  mutate(
    delta_ct = CT - igg_ct,
    fold_enrichment = 2^(igg_ct - CT)
  )

# --- PLOT 1: LINEAR SCALE ---
enrichment_plot <- ggplot(enrichment_data, aes(x = `Sample Name`, y = fold_enrichment, fill = `Target Name`)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(title = paste("Enrichment Over IgG (", cond_name, ")", sep = ""),
       x = "Sample",
       y = "Fold Enrichment (vs IgG, linear)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = "2025-08-18_C&R16_qPCR_cond1_enrichment_plot_linear.png",
  plot = enrichment_plot,
  path = "~/Documents/postdoc/results/pioneering_activity_project/CUT&RUN/qPCR/2025-08-15_CutAndRun16_Vera/",
  width = 10, height = 6, dpi = 300
)

# --- PLOT 2: LOG SCALE ---
enrichment_plot_log <- enrichment_plot +
  scale_y_log10() +
  labs(y = "Fold Enrichment (vs IgG, log10 scale)") +
  annotation_logticks(sides = "l")

ggsave(
  filename = "2025-08-18_C&R16_qPCR_cond1_enrichment_plot_log.png",
  plot = enrichment_plot_log,
  path = "~/Documents/postdoc/results/pioneering_activity_project/CUT&RUN/qPCR/2025-08-15_CutAndRun16_Vera/",
  width = 10, height = 6, dpi = 300
)


# --- CONDITION 2 ONLY ---

# Define condition name and IgG sample name
cond_name <- "cond2"
igg_name  <- "IgG_cond2"

# Get mean IgG Ct per target
igg_ct <- filtered_data %>%
  filter(`Sample Name` == igg_name) %>%
  group_by(`Target Name`) %>%
  summarise(igg_ct = mean(CT), .groups = "drop")

# Join IgG Ct to condition 2 samples (excluding IgG itself)
enrichment_data <- filtered_data %>%
  filter(grepl(cond_name, `Sample Name`), `Sample Name` != igg_name) %>%
  left_join(igg_ct, by = "Target Name") %>%
  mutate(
    delta_ct = CT - igg_ct,
    fold_enrichment = 2^(igg_ct - CT)
  )

# --- PLOT 1: LINEAR SCALE ---
enrichment_plot <- ggplot(enrichment_data, aes(x = `Sample Name`, y = fold_enrichment, fill = `Target Name`)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(title = paste("Enrichment Over IgG (", cond_name, ")", sep = ""),
       x = "Sample",
       y = "Fold Enrichment (vs IgG, linear)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = "2025-08-18_C&R16_qPCR_cond2_enrichment_plot_linear.png",
  plot = enrichment_plot,
  path = "~/Documents/postdoc/results/pioneering_activity_project/CUT&RUN/qPCR/2025-08-15_CutAndRun16_Vera/",
  width = 10, height = 6, dpi = 300
)

# --- PLOT 2: LOG SCALE ---
enrichment_plot_log <- enrichment_plot +
  scale_y_log10() +
  labs(y = "Fold Enrichment (vs IgG, log10 scale)") +
  annotation_logticks(sides = "l")

ggsave(
  filename = "2025-08-18_C&R16_qPCR_cond2_enrichment_plot_log.png",
  plot = enrichment_plot_log,
  path = "~/Documents/postdoc/results/pioneering_activity_project/CUT&RUN/qPCR/2025-08-15_CutAndRun16_Vera/",
  width = 10, height = 6, dpi = 300
)


# --- CONDITION 3 (NOlibrary) ONLY ---

# Define condition name and IgG sample name
cond_name <- "cond3_NOlibrary"
igg_name  <- "IgG_cond3_NOlibrary"

# Get mean IgG Ct per target
igg_ct <- filtered_data %>%
  filter(`Sample Name` == igg_name) %>%
  group_by(`Target Name`) %>%
  summarise(igg_ct = mean(CT), .groups = "drop")

# Join IgG Ct to condition 3_NOlibrary samples (excluding IgG itself)
enrichment_data <- filtered_data %>%
  filter(grepl(cond_name, `Sample Name`), `Sample Name` != igg_name) %>%
  left_join(igg_ct, by = "Target Name") %>%
  mutate(
    delta_ct = CT - igg_ct,
    fold_enrichment = 2^(igg_ct - CT)
  )

# --- PLOT 1: LINEAR SCALE ---
enrichment_plot <- ggplot(enrichment_data, aes(x = `Sample Name`, y = fold_enrichment, fill = `Target Name`)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(title = paste("Enrichment Over IgG (", cond_name, ")", sep = ""),
       x = "Sample",
       y = "Fold Enrichment (vs IgG, linear)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = "2025-08-18_C&R16_qPCR_cond3_NOlibrary_enrichment_plot_linear.png",
  plot = enrichment_plot,
  path = "~/Documents/postdoc/results/pioneering_activity_project/CUT&RUN/qPCR/2025-08-15_CutAndRun16_Vera/",
  width = 10, height = 6, dpi = 300
)

# --- PLOT 2: LOG SCALE ---
enrichment_plot_log <- enrichment_plot +
  scale_y_log10() +
  labs(y = "Fold Enrichment (vs IgG, log10 scale)") +
  annotation_logticks(sides = "l")

ggsave(
  filename = "2025-08-18_C&R16_qPCR_cond3_NOlibrary_enrichment_plot_log.png",
  plot = enrichment_plot_log,
  path = "~/Documents/postdoc/results/pioneering_activity_project/CUT&RUN/qPCR/2025-08-15_CutAndRun16_Vera/",
  width = 10, height = 6, dpi = 300
)
