
# ------------------------------------------------------------------------------
# Script: 03_stable_genes_protein_vs_rna_ratios.qmd
# Author: Luke
#
# Overview:
#   - Loads the RNA stable genes list and processed proteomics data.
#   - Computes differences between protein and RNA log₂ fold-changes.
#   - Generates density and variance plots to compare RNA vs. protein data.
#   - Filters candidate genes exceeding a set threshold and creates a 
#     formatted table for downstream analysis.
#
# Inputs:
#   - Stable genes list (RNA): translation_on_demand/data/mm11/processed/stable_genes_list.rds
#   - Processed proteomics data: translation_on_demand/data/mm11/processed/processed_prot_data.rds
#
# Outputs:
#   - Passed genes list (RDS): translation_on_demand/data/mm11/processed/passed_genes_list.rds
#   - Formatted CSV table: translation_on_demand/data/mm11/processed/passed_genes_table.csv
# ------------------------------------------------------------------------------
# -------------------------------
# 03_stable_genes_protein_vs_rna_ratio
# -------------------------------
# Set working directory and update library paths
setwd('/cellfile/cellnet/mESC_differentiation')
new_library_path <- "/cellfile/cellnet/mESC_differentiation/Rlibs/Rlibs_433_Fabian"
.libPaths(c(new_library_path, .libPaths()))

# Load necessary libraries
library(tidyverse)

# -------------------------------
# Load RNA Stable Genes and Protein Data
# -------------------------------
# Load stable genes list (RDS file saved earlier)
stable_genes_list <- readRDS("translation_on_demand/data/mm11/processed/stable_genes_list.rds")

# Load processed protein data (RDS file saved earlier)
prot_data <- readRDS("translation_on_demand/data/mm11/processed/processed_prot_data.rds")

# Preview the structure of the stable genes list
cat("Stable genes list names (contrast comparisons):\n")
print(names(stable_genes_list))
#glimpse(stable_genes_list[[1]])  # preview one contrast

# Preview the processed protein data
#glimpse(prot_data)

# -------------------------------
# Calculating Protein - RNA Ratios (Log2FC differences)
# -------------------------------
# Define mapping: RNA contrast -> corresponding protein contrast column.
contrast_mapping <- c(
  "1h_vs_0h"   = "FC_0m_1h",   # RNA 1h_vs_0h corresponds to protein FC_0m_1h
  "6h_vs_1h"   = "FC_1h_6h",
  "12h_vs_6h"  = "FC_6h_12h",
  "24h_vs_12h" = "FC_12h_24h",
  "36h_vs_24h" = "FC_24h_36h",
  "48h_vs_36h" = "FC_36h_48h",
  "72h_vs_48h" = "FC_48h_72h"
)

logFC_diff_list <- list()

for (contrast_name in names(contrast_mapping)) {
  # Extract RNA stable gene data for this contrast.
  rna_df <- stable_genes_list[[contrast_name]]
  
  # If mgi_symbol is not a column (it's stored as rownames), convert rownames to a column.
  if (!"mgi_symbol" %in% colnames(rna_df)) {
    rna_df <- rna_df %>% rownames_to_column("mgi_symbol")
  }
  
  # Get the corresponding protein contrast column.
  prot_col <- contrast_mapping[contrast_name]
  
  # Extract the relevant protein data: GeneID and the contrast column.
  prot_df <- prot_data %>%
    dplyr::select(GeneID, !!sym(prot_col))
  
  # Join RNA and protein data on gene identifiers.
  merged_df <- inner_join(rna_df, prot_df, by = c("mgi_symbol" = "GeneID"))
  
  # Compute the difference: (Protein logFC) - (RNA log2FoldChange) and flag if diff > 1.
  merged_df <- merged_df %>%
    mutate(logFC_diff = !!sym(prot_col) - log2FoldChange,
           exceeds_threshold = if_else(logFC_diff > 1, TRUE, FALSE))
  
  # Save the result.
  logFC_diff_list[[contrast_name]] <- merged_df
  
  # Report the number of matched genes.
  cat("Contrast:", contrast_name, "-> Matched genes:", nrow(merged_df), "\n")
}

# Optionally, inspect one of the contrasts.
#glimpse(logFC_diff_list[["1h_vs_0h"]])

# -------------------------------
# Plotting Ratio Density Curves
# -------------------------------
# Combine matched stable gene data from all contrasts into one data frame.
diff_df <- bind_rows(lapply(names(logFC_diff_list), function(comp) {
  df <- logFC_diff_list[[comp]]
  df$Comparison <- comp  # add a column with the contrast name
  df
}))

# Reformat comparison labels and set factor levels.
diff_df <- diff_df %>%
  mutate(Comparison = gsub("_vs_", " vs ", Comparison)) %>%
  mutate(Comparison = factor(Comparison, levels = c("1h vs 0h", 
                                                    "6h vs 1h", 
                                                    "12h vs 6h", 
                                                    "24h vs 12h", 
                                                    "36h vs 24h", 
                                                    "48h vs 36h", 
                                                    "72h vs 48h")))

# Remove rows with non-finite logFC_diff values.
diff_df <- diff_df %>% filter(is.finite(logFC_diff))

# Create a density plot of the logFC_diff values.
p_diff <- ggplot(diff_df, aes(x = logFC_diff, color = Comparison)) +
  geom_density(linewidth = 1.5, na.rm = TRUE) +
  labs(title = "Density of log₂FC Differences for Matched Stable Genes",
       x = "Difference (Proteomics - RNA log₂FC)",
       y = "Density",
       color = "Timepoint Comparison") +
  theme_classic(base_size = 24) +
  theme(axis.text.x = element_text(size = 22, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 22),
        axis.title.x = element_text(size = 26),
        axis.title.y = element_text(size = 26),
        plot.title = element_text(size = 28, face = "bold"),
        legend.title = element_text(size = 26),
        legend.text = element_text(size = 24),
        legend.position = c(0.8, 0.8),
        legend.position.inside = TRUE) +
  guides(color = guide_legend(reverse = TRUE, ncol = 1))

print(p_diff)

# -------------------------------
# Filtering ToD Candidates
# -------------------------------
# For each contrast, extract genes that pass the filter (exceeds_threshold == TRUE)
# and order them in descending order of logFC_diff.
passed_genes_list <- lapply(names(logFC_diff_list), function(comp) {
  df <- logFC_diff_list[[comp]]
  df %>% 
    filter(exceeds_threshold == TRUE) %>% 
    mutate(logFC_diff = as.numeric(logFC_diff)) %>% 
    arrange(desc(logFC_diff))
})
names(passed_genes_list) <- names(logFC_diff_list)

# Report the number of passing genes for each contrast.
for (comp in names(passed_genes_list)) {
  cat("Contrast:", comp, "-> Passing genes:", nrow(passed_genes_list[[comp]]), "\n")
}

# Optionally, inspect one contrast.
#glimpse(passed_genes_list[["1h_vs_0h"]])

# Save the passed genes list to an RDS file
saveRDS(passed_genes_list, "translation_on_demand/data/mm11/processed/passed_genes_list.rds")

# -------------------------------
# Density Plot of log₂FC Distributions
# -------------------------------
library(tidyverse)

# Define mapping: RNA contrast -> corresponding protein column.
prot_map <- c(
  "1h_vs_0h"   = "FC_0m_1h",
  "6h_vs_1h"   = "FC_1h_6h",
  "12h_vs_6h"  = "FC_6h_12h",
  "24h_vs_12h" = "FC_12h_24h",
  "36h_vs_24h" = "FC_24h_36h",
  "48h_vs_36h" = "FC_36h_48h",
  "72h_vs_48h" = "FC_48h_72h"
)
contrast_names <- names(prot_map)

# Build a combined data frame with RNA and Protein log₂FC values (only for matched genes).
plot_list <- list()
for (contrast in contrast_names) {
  
  # RNA: from full gene list for this contrast.
  rna_df <- all_genes_list[[contrast]]
  if (!"mgi_symbol" %in% colnames(rna_df)) {
    rna_df <- rna_df %>% rownames_to_column("mgi_symbol")
  }
  
  # Determine matched genes between RNA and Protein.
  matched_genes <- intersect(rna_df$mgi_symbol, prot_data$GeneID)
  
  # Filter RNA data for matched genes.
  rna_df <- rna_df %>% filter(mgi_symbol %in% matched_genes)
  rna_plot <- tibble(
    logFC = rna_df$log2FoldChange,
    DataType = "RNA",
    Contrast = contrast
  )
  
  # Protein: from processed protein data; filter for matched genes.
  prot_df <- prot_data %>% filter(GeneID %in% matched_genes)
  prot_plot <- tibble(
    logFC = prot_df[[prot_map[contrast]]],
    DataType = "Protein",
    Contrast = contrast
  )
  
  # Combine RNA and Protein values for this contrast.
  plot_list[[contrast]] <- bind_rows(rna_plot, prot_plot)
}
plot_df <- bind_rows(plot_list)

# Reformat contrast labels and order them.
plot_df <- plot_df %>%
  mutate(Contrast = gsub("_vs_", " vs ", Contrast),
         Contrast = factor(Contrast, levels = c("1h vs 0h", 
                                                  "6h vs 1h", 
                                                  "12h vs 6h",
                                                  "24h vs 12h", 
                                                  "36h vs 24h", 
                                                  "48h vs 36h", 
                                                  "72h vs 48h")))

# Create the density plot.
p_logFC <- ggplot(plot_df, aes(x = logFC, color = DataType)) +
  geom_density(linewidth = 1, na.rm = TRUE) +
  facet_wrap(~ Contrast, scales = "free") +
  labs(title = "Density of log₂FC Distributions by Contrast (Matched Genes Only)",
       x = "log₂ Fold Change",
       y = "Density",
       color = "Data Type") +
  theme_classic(base_size = 24) +
  theme(
    strip.background = element_blank(),  # remove black boxes around facet titles
    strip.text = element_text(size = 22),
    axis.text.x = element_text(size = 22, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 22),
    axis.title.x = element_text(size = 26),
    axis.title.y = element_text(size = 26),
    plot.title = element_text(size = 28, face = "bold"),
    legend.title = element_text(size = 26),
    legend.text = element_text(size = 24),
    legend.position = c(0.95, 0.05),  # place legend in the bottom-right (inside plot)
    legend.justification = c("right", "bottom")
  ) +
  guides(color = guide_legend(reverse = TRUE, ncol = 1))
print(p_logFC)

#### Barplot of Global Variance for Each Contrast ####

# Ensure each RNA full gene data frame has an explicit mgi_symbol column.
all_genes_list <- lapply(all_genes_list, function(df) {
  if (!"mgi_symbol" %in% colnames(df)) {
    df <- df %>% rownames_to_column("mgi_symbol")
  }
  return(df)
})

# Compute variance for each contrast using only the matched genes.
variance_list <- lapply(contrast_names, function(comp) {
  # Extract RNA gene symbols from the RNA data for this contrast.
  rna_symbols <- all_genes_list[[comp]] %>% pull(mgi_symbol) %>% unique()
  # Extract protein gene IDs from the processed protein data.
  prot_symbols <- prot_data$GeneID %>% unique()
  
  # Matched genes: intersection of RNA symbols and protein GeneIDs.
  matched_genes <- intersect(rna_symbols, prot_symbols)
  
  # Subset RNA values for the matched genes.
  rna_values <- all_genes_list[[comp]] %>% 
    filter(mgi_symbol %in% matched_genes) %>% 
    pull(log2FoldChange)
  
  # Subset Protein values for the matched genes using the appropriate contrast column.
  prot_values <- prot_data %>% 
    filter(GeneID %in% matched_genes) %>% 
    pull(!!sym(prot_map[comp]))
  
  tibble(
    Contrast = comp,
    DataType = c("RNA", "Protein"),
    Variance = c(var(rna_values, na.rm = TRUE),
                 var(prot_values, na.rm = TRUE))
  )
})
var_df <- bind_rows(variance_list)

# Reformat contrast labels.
var_df <- var_df %>%
  mutate(Contrast = gsub("_vs_", " vs ", Contrast),
         Contrast = factor(Contrast, levels = c("1h vs 0h", 
                                                  "6h vs 1h", 
                                                  "12h vs 6h",
                                                  "24h vs 12h", 
                                                  "36h vs 24h", 
                                                  "48h vs 36h", 
                                                  "72h vs 48h")))

# Create the grouped barplot.
p_variance <- ggplot(var_df, aes(x = Contrast, y = Variance, fill = DataType)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  labs(title = "Global Variance of log₂FC by Contrast (Matched Genes Only)",
       x = "Contrast",
       y = "Variance",
       fill = "Data Type") +
  theme_classic(base_size = 24) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 22),
        axis.text.y = element_text(size = 22),
        axis.title.x = element_text(size = 26),
        axis.title.y = element_text(size = 26),
        plot.title = element_text(size = 28, face = "bold"),
        legend.title = element_text(size = 26),
        legend.text = element_text(size = 24))
print(p_variance)

