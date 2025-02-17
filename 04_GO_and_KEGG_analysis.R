
# ------------------------------------------------------------------------------
# Script: 04_GO_and_KEGG_analysis.qmd
# Author: Luke
#
# Overview:
#   - Performs functional enrichment analyses using GO and KEGG on RNA data.
#   - Loads gene lists to define the background and candidate genes.
#   - Runs GO enrichment (BP, CC, MF) and creates dot plots for each ontology.
#   - Performs KEGG pathway enrichment analysis and generates corresponding plots.
#
# Inputs:
#   - All genes logFC list (RNA): translation_on_demand/data/mm11/processed/all_genes_logFC_list.rds
#   - Passed genes list (stable genes RNA): translation_on_demand/data/mm11/processed/passed_genes_list.rds
#
# Outputs:
#   - GO enrichment plots (BP, CC, MF) saved as PNG files.
#   - KEGG enrichment plots displayed for each contrast.
# ------------------------------------------------------------------------------
# -------------------------------
# 04_GO_and_KEGG_analysis
# -------------------------------
# Set working directory and update library paths
setwd('/cellfile/cellnet/mESC_differentiation')
new_library_path <- "/cellfile/cellnet/mESC_differentiation/Rlibs/Rlibs_433_Fabian"
.libPaths(c(new_library_path, .libPaths()))

# Load necessary libraries
library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)

# Load the full gene logFC list (for all contrasts) from RNA
all_genes_list <- readRDS("translation_on_demand/data/mm11/processed/all_genes_logFC_list.rds")
# Load the passed gene logFC list (for stable genes) from RNA
passed_genes_list <- readRDS("translation_on_demand/data/mm11/processed/passed_genes_list.rds")

# ------------------------------
# GO Analysis Plotting: 1 plot for each type (MF, CC, BP): Define Background Gene List
# ------------------------------
# Combine gene symbols from all RNA contrasts
background_genes <- bind_rows(lapply(all_genes_list, function(df) {
  # Ensure gene symbols are in a column
  if (!"mgi_symbol" %in% colnames(df)) df <- df %>% rownames_to_column("mgi_symbol")
  df
})) %>% pull(mgi_symbol) %>% unique()

# Convert background gene symbols to Entrez IDs
bg_entrez <- bitr(background_genes, fromType = "SYMBOL", 
                  toType = "ENTREZID", OrgDb = org.Mm.eg.db)
cat("Number of background (common) genes (converted):", nrow(bg_entrez), "\n\n")

# ------------------------------
# Run GO Enrichment Analysis for Each Ontology Across All Contrasts
# ------------------------------
comparison_names <- c("1h_vs_0h", "6h_vs_1h", "12h_vs_6h", 
                     "24h_vs_12h", "36h_vs_24h", "48h_vs_36h", "72h_vs_48h")
ontologies <- c("BP", "CC", "MF")
results_list <- list()

for (ont in ontologies) {
  comp_results <- lapply(comparison_names, function(comp) {
    if (!comp %in% names(passed_genes_list)) {
      cat("Comparison", comp, "not found in passed_genes_list.\n")
      return(NULL)
    }
    
    # Get candidate gene symbols from the passed genes list
    candidate_genes <- passed_genes_list[[comp]] %>% 
      { if("mgi_symbol" %in% colnames(.)) . else rownames_to_column(., "mgi_symbol") } %>%
      pull(mgi_symbol) %>% unique()
    
    if (length(candidate_genes) == 0) {
      cat("No candidate genes found for", comp, "\n")
      return(NULL)
    }
    
    # Convert candidate genes to Entrez IDs
    cand_entrez <- bitr(candidate_genes, fromType = "SYMBOL", 
                        toType = "ENTREZID", OrgDb = org.Mm.eg.db)
    
    cat("-------------------------------------------------\n")
    cat("Comparison:", comp, "\n")
    cat("Number of candidate genes:", length(candidate_genes), "\n")
    cat("Candidate genes (converted):", nrow(cand_entrez), "\n\n")
    
    # Run GO enrichment with stricter thresholds
    ego <- enrichGO(gene         = cand_entrez$ENTREZID,
                    universe     = bg_entrez$ENTREZID,
                    OrgDb        = org.Mm.eg.db,
                    ont          = ont,
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.1,
                    qvalueCutoff = 0.1,
                    readable     = TRUE)
    
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      df <- as.data.frame(ego)
      df$Contrast <- comp
      return(df)
    } else {
      return(NULL)
    }
  })
  
  # Combine the results for the current ontology (remove NULLs)
  comp_results <- bind_rows(comp_results)
  results_list[[ont]] <- comp_results
}

# ------------------------------
# Create Combined Dot Plots for Each Ontology Using ggplot2
# ------------------------------
plots <- list()
for (ont in ontologies) {
  df <- results_list[[ont]]
  if (!is.null(df) && nrow(df) > 0) {
    # Reformat contrast labels and sort by p.adjust within each Contrast
    df <- df %>% 
      mutate(Contrast = gsub("_vs_", " vs ", Contrast),
             Contrast = factor(Contrast, levels = c("1h vs 0h", "6h vs 1h", 
                                                  "12h vs 6h", "24h vs 12h", 
                                                  "36h vs 24h", "48h vs 36h", "72h vs 48h"))) %>%
      group_by(Contrast) %>%
      slice_max(order_by = -p.adjust, n = 5) %>%  # Changed to order by p.adjust
      arrange(Contrast, p.adjust) %>%  # Sort by p.adjust within each contrast
      mutate(Description = factor(Description, levels = rev(unique(Description)))) %>%  # Reverse to put most significant at top
      ungroup()
    
    p <- ggplot(df, aes(x = Count, y = reorder(Description, -p.adjust))) +  # Reorder by significance
      geom_point(aes(size = Count, color = p.adjust)) +
      scale_color_gradient(low = "red", high = "blue", name = "Adjusted\np-value") +
      scale_size_continuous(name = "Gene\ncount", labels = function(x) sprintf("%.0f", x)) +
      facet_wrap(~ Contrast, scales = "free_y", ncol = 2) +  # Changed to 2 columns
      theme_bw() +
      ggtitle(paste("GO Enrichment -", ont)) +
      theme(
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8, lineheight = 0.8),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        legend.position = "right",
        panel.spacing = unit(1, "lines"),
        plot.margin = margin(10, 10, 10, 10)
      ) +
      labs(x = "Gene Count", y = "GO Term") +
      scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 25))
    
    plots[[ont]] <- p
  }
}

# Example: Print the BP, MF, and CC plots.
print(plots[["BP"]])
print(plots[["MF"]])
print(plots[["CC"]])

# Save the BP plot
ggsave("translation_on_demand/scripts/m11_scripts/plots/GO_Enrichment_BP.png", plot = plots[["BP"]], width = 12, height = 10, dpi = 300)

# Save the CC plot
ggsave("translation_on_demand/scripts/m11_scripts/plots/GO_Enrichment_CC.png", plot = plots[["CC"]], width = 12, height = 10, dpi = 300)

# Save the MF plot
ggsave("translation_on_demand/scripts/m11_scripts/plots/GO_Enrichment_MF.png", plot = plots[["MF"]], width = 12, height = 10, dpi = 300)

# ------------------------------
# KEGG Analysis
# ------------------------------
# --- Step 0. Define the Background Gene List ---
# Combine all gene symbols across contrasts; assume gene symbols are in the 'mgi_symbol' column
# Convert rownames to a column ("gene_symbol") for each contrast and combine them.
background_genes <- bind_rows(lapply(all_genes_list, function(df) {
  df %>% rownames_to_column("gene_symbol")
})) %>% 
  pull(gene_symbol) %>% 
  unique()

# Convert background gene symbols to Entrez IDs.
bg_entrez <- bitr(background_genes, fromType = "SYMBOL", 
                  toType = "ENTREZID", OrgDb = org.Mm.eg.db)
cat("Number of background (common) genes (converted):", nrow(bg_entrez), "\n\n")

# --- Step 1: Run KEGG Enrichment Analysis for Each Comparison ---
# Define the comparisons to run (using those contrasts with clear matching).
comparison_names <- c("1h_vs_0h", "6h_vs_1h", "12h_vs_6h", 
                      "24h_vs_12h", "36h_vs_24h", "48h_vs_36h", "72h_vs_48h")

for (comp in comparison_names) {
  
  if (!comp %in% names(passed_genes_list)) {
    cat("Comparison", comp, "not found in passed_genes_list.\n")
    next
  }
  
  # Extract candidate gene symbols from the passed genes list.
  # If the gene symbols are stored as rownames, convert them into a column.
  candidate_genes <- passed_genes_list[[comp]] %>% 
    { if ("mgi_symbol" %in% colnames(.)) . else rownames_to_column(., "mgi_symbol") } %>% 
    pull(mgi_symbol) %>% unique()
  
  if (length(candidate_genes) == 0) {
    cat("No candidate genes found for", comp, "\n")
    next
  }
  
  # Convert candidate gene symbols to Entrez IDs.
  cand_entrez <- bitr(candidate_genes, fromType = "SYMBOL", 
                      toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  
  cat("-------------------------------------------------\n")
  cat("Comparison:", comp, "\n")
  cat("Number of candidate genes:", length(candidate_genes), "\n")
  cat("Candidate genes (converted):", nrow(cand_entrez), "\n\n")
  
  # Run KEGG enrichment analysis using enrichKEGG.
  kegg_res <- enrichKEGG(gene         = cand_entrez$ENTREZID,
                         universe     = bg_entrez$ENTREZID,
                         organism     = "mmu",       # mouse KEGG organism code
                         pvalueCutoff = 0.1,         # lenient threshold for exploration
                         qvalueCutoff = 0.1)
  
  # Remove unwanted text from pathway descriptions.
  if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res)) > 0) {
    kegg_res@result$Description <- gsub(" - Mus musculus.*", "", kegg_res@result$Description)
  }
  
  # Create the dot plot.
  if (nrow(as.data.frame(kegg_res)) > 0) {
    p_kegg <- dotplot(kegg_res, showCategory = 20) +
      ggtitle(paste("KEGG Pathway Enrichment (", comp, ")", sep = ""))
  } else {
    p_kegg <- ggplot() +
      geom_blank() +
      ggtitle(paste("No enriched KEGG pathways for", comp))
  }
  
  # Save the plot to a global variable (e.g. p_1h_vs_0h_KEGG).
  assign(paste0("p_", comp, "_KEGG"), p_kegg, envir = .GlobalEnv)
  
  cat("KEGG enrichment analysis completed for comparison:", comp, "\n\n")
}

# --- Step 2: Print Each KEGG Enrichment Plot ---
print(p_1h_vs_0h_KEGG)
print(p_6h_vs_1h_KEGG)
print(p_12h_vs_6h_KEGG)
print(p_24h_vs_12h_KEGG)
print(p_36h_vs_24h_KEGG)
print(p_48h_vs_36h_KEGG)
print(p_72h_vs_48h_KEGG)

