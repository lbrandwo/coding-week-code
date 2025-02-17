
# ------------------------------------------------------------------------------
# Script: 02_DESeq2_stable_genes.qmd
# Author: Luke
#
# Overview:
#   Performs differential expression analysis using DESeq2 on matched RNA data.
#   - Prepares count data and sample metadata.
#   - Runs DESeq2 to identify stable genes (with minimal log₂ fold-change) 
#     across various contrasts.
#   - Extracts log₂ fold changes for all genes for each contrast.
#
# Inputs:
#   - Matched RNA data: translation_on_demand/data/mm11/processed/processed_rna_data_matched.rds
#
# Outputs:
#   - Stable genes list: translation_on_demand/data/mm11/processed/stable_genes_list.rds
#   - All genes log₂FC list: translation_on_demand/data/mm11/processed/all_genes_logFC_list.rds
# ------------------------------------------------------------------------------
# -------------------------------
# 02_DESeq2_stable_genes
# -------------------------------
# Set working directory and update library paths
setwd('/cellfile/cellnet/mESC_differentiation')
new_library_path <- "/cellfile/cellnet/mESC_differentiation/Rlibs/Rlibs_433_Fabian"
.libPaths(c(new_library_path, .libPaths()))

# Load necessary libraries
library(tidyverse)
library(DESeq2)

# ------------------------------
# Load Processed RNA Data
# ------------------------------
rna_data <- readRDS("translation_on_demand/data/mm11/processed/processed_rna_data_matched.rds")
# The rna_data contains:
#   - mgi_symbol (gene name)
#   - Count columns for each timepoint (e.g. "X0h_rep1", "X0h_rep2", "X1h_rep1", "X1h_rep2", etc.)

# ------------------------------
# Prepare Count Data for DESeq2
# ------------------------------
# Extract count columns (begin with "X") and convert to a matrix.
count_data <- rna_data %>% 
  dplyr::select(starts_with("X")) %>% 
  as.matrix()

# Use mgi_symbol as the row names.
rownames(count_data) <- rna_data$mgi_symbol

# ------------------------------
# Create Sample Metadata (colData)
# ------------------------------
sample_names <- colnames(count_data)
colData <- tibble(sample = sample_names) %>%
  mutate(
    # Remove the leading "X" and extract the timepoint (the part before the underscore)
    timepoint = str_extract(sample, "(?<=X)[^_]+") %>% str_trim(),
    # Extract replicate info (e.g. "rep1" or "rep2")
    replicate = str_extract(sample, "(?i)(?<=_rep)\\d+")
  ) %>%
  mutate(
    timepoint = factor(timepoint, levels = c("0h", "1h", "6h", "12h", "24h", "36h", "48h", "72h")),
    replicate = factor(replicate)
  ) %>%
  column_to_rownames("sample")

# Check that the number of samples matches:
stopifnot(ncol(count_data) == nrow(colData))

# ------------------------------
# Create DESeq2 Dataset and Run DESeq2
# ------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = count_data,    # genes x samples count matrix
  colData   = colData,        # sample metadata
  design    = ~ timepoint     # experimental design
)
dds <- DESeq(dds)             # run DESeq2 analysis

# ------------------------------
# Define Contrasts and Extract Stable Genes
# ------------------------------
contrasts <- list(
  "1h_vs_0h"   = c("timepoint", "1h", "0h"),
  "6h_vs_1h"   = c("timepoint", "6h", "1h"),
  "12h_vs_6h"  = c("timepoint", "12h", "6h"),
  "24h_vs_12h" = c("timepoint", "24h", "12h"),
  "36h_vs_24h" = c("timepoint", "36h", "24h"),
  "48h_vs_36h" = c("timepoint", "48h", "36h"),
  "72h_vs_48h" = c("timepoint", "72h", "48h")
)

stable_genes_list <- list()
for (comp_name in names(contrasts)) {
  contrast <- contrasts[[comp_name]]
  
  # Get DESeq2 results with an lfcThreshold of 1
  res_stable <- results(dds, contrast = contrast,
                        lfcThreshold = 1,
                        altHypothesis = "greaterAbs")
  res_df <- as.data.frame(res_stable)
  
  # If padj is missing, compute it manually.
  if (!"padj" %in% colnames(res_df)) {
    res_df$padj <- p.adjust(res_df$pvalue, method = "BH")
  }
  
  # Annotate results with mgi_symbol
  res_df$mgi_symbol <- rowData(dds)$mgi_symbol[match(rownames(res_df), rownames(dds))]
  
  # Filter genes with 95% CI entirely within [-1, +1] and non-significant change.
  stable <- res_df %>%
    filter(padj > 0.05,  # using padj if available
           (log2FoldChange - 1.96 * lfcSE) >= -1.5,
           (log2FoldChange + 1.96 * lfcSE) <= 1.5)
  
  stable_genes_list[[comp_name]] <- stable
  
  cat("Comparison:", comp_name, "\n")
  cat("Number of stable genes (95% CI within ±1):", nrow(stable), "\n\n")
}

# Save the stable genes list
saveRDS(stable_genes_list, "translation_on_demand/data/mm11/processed/stable_genes_list.rds")

# ------------------------------
# Define Contrasts and Extract All Genes Log₂FCs
# ------------------------------
contrasts <- list(
  "1h_vs_0h"   = c("timepoint", "1h", "0h"),
  "6h_vs_1h"   = c("timepoint", "6h", "1h"),
  "12h_vs_6h"  = c("timepoint", "12h", "6h"),
  "24h_vs_12h" = c("timepoint", "24h", "12h"),
  "36h_vs_24h" = c("timepoint", "36h", "24h"),
  "48h_vs_36h" = c("timepoint", "48h", "36h"),
  "72h_vs_48h" = c("timepoint", "72h", "48h")
)

all_genes_list <- list()
for (comp_name in names(contrasts)) {
  contrast <- contrasts[[comp_name]]
  
  # Get DESeq2 results with an lfcThreshold of 1 (for all genes)
  res_all <- results(dds, contrast = contrast,
                     lfcThreshold = 1,
                     altHypothesis = "greaterAbs")
  res_df_all <- as.data.frame(res_all)
  
  # Annotate with mgi_symbol; note that mgi_symbol is stored as rownames.
  res_df_all$mgi_symbol <- rowData(dds)$mgi_symbol[match(rownames(res_df_all), rownames(dds))]
  
  all_genes_list[[comp_name]] <- res_df_all
  
  cat("Comparison:", comp_name, "\n")
  cat("Number of genes (all):", nrow(res_df_all), "\n\n")
}

# Save the complete list of gene log₂FCs for each contrast for later plotting.
saveRDS(all_genes_list, "translation_on_demand/data/mm11/processed/all_genes_logFC_list.rds")
