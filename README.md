# Multi-Omics Analysis Pipeline for mESC Differentiation

This repository contains a multi-step analysis pipeline for integrating RNA-seq and proteomics data from a study on mouse embryonic stem cell (mESC) differentiation. The workflow includes data preprocessing, differential expression analysis, protein versus RNA expression comparisons, and functional enrichment analysis using Gene Ontology (GO) and KEGG pathways.

Files are .qmd at this exploratory stage where things are changed regularly for testing. Viewing and plotting will be removed once exploratory analysis is concluded.

## Project Structure

- **01_preprocess_raw_data.qmd**  
  Preprocesses raw RNA-seq and proteomics data by filtering, normalizing, and saving processed datasets.  
  **Inputs:**  
  - RNA-seq counts CSV file (`translation_on_demand/counts/mm11_counts_mgi.csv`)  
  - Proteomics Excel file (`/cellfile/cellnet/mESC_differentiation/translation_on_demand/data/mm9/raw/proteomic_data_raw.xlsx`)  
  **Outputs:**  
  - Processed RNA data (RDS)  
  - Processed proteomics data (RDS)  
  - Matched RNA-protein gene list (RDS)

- **02_DESeq2_stable_genes.qmd**  
  Runs differential expression analysis using DESeq2 on the processed RNA data to identify stable genes across timepoints.  
  **Inputs:**  
  - Processed and matched RNA data  
  **Outputs:**  
  - DESeq2 results for stable genes (RDS)  
  - All genes log₂ fold-change lists (RDS)

- **03_stable_genes_protein_vs_rna_ratios.qmd**  
  Compares protein and RNA log₂ fold-changes for stable genes, generates density and variance plots, and filters candidate genes with notable differences.  
  **Inputs:**  
  - Stable genes list from RNA  
  - Processed proteomics data  
  **Outputs:**  
  - Density and variance plots  
  - Passed (candidate) genes list (RDS)  
  - Formatted gene table (CSV)

- **04_GO_and_KEGG_analysis.qmd**  
  Performs functional enrichment analyses using GO and KEGG on candidate genes to identify enriched biological processes and pathways.  
  **Inputs:**  
  - All genes log₂ fold-change list (RNA)  
  - Passed genes list (RNA)  
  **Outputs:**  
  - GO enrichment dot plots (BP, CC, MF; PNG files)  
  - KEGG enrichment plots for each contrast

## Setup

Before running the scripts, ensure you have the following:

- **Software Requirements:**
  - R (version ≥ 4.0) and RStudio or a Quarto-compatible editor.
  - Required R packages: `tidyverse`, `DESeq2`, `clusterProfiler`, `org.Mm.eg.db`, `enrichplot`, etc.

- **Data Files:**
  - RNA-seq data: `translation_on_demand/counts/mm11_counts_mgi.csv`
  - Proteomics data: `/cellfile/cellnet/mESC_differentiation/translation_on_demand/data/mm9/raw/proteomic_data_raw.xlsx`

- **Configuration:**
  - Update the working directory and library paths in each script to match your local setup.

## Running the Pipeline

1. **Preprocess Data:**  
   Run `01_preprocess_raw_data.qmd` to filter and process the raw RNA-seq and proteomics datasets.

2. **Differential Expression Analysis:**  
   Execute `02_DESeq2_stable_genes.qmd` to run DESeq2 and identify stable genes in the RNA data.

3. **Protein vs RNA Comparison:**  
   Run `03_stable_genes_protein_vs_rna_ratios.qmd` to compute differences between protein and RNA log₂ fold-changes and to generate comparison plots.

4. **Functional Enrichment:**  
   Finally, execute `04_GO_and_KEGG_analysis.qmd` to perform GO and KEGG enrichment analyses and visualize the results.

## Troubleshooting

- **GitHub Authentication Issues:**  
  If you have trouble pushing changes to GitHub (e.g., "username or password is wrong"), make sure you are using a Personal Access Token (PAT) in place of your GitHub password for HTTPS operations, or consider setting up SSH keys.

- **Data Path Errors:**  
  Verify that all file paths in the scripts are correct and accessible from your working environment.

## License

This project is licensed under the [MIT License](LICENSE).

## Contact

For questions or further support, please contact Luke at lbrandw1@uni-koeln.de.
