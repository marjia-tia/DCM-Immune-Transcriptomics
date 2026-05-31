# ============================================================
# Script 00: Package Installation
# Run ONCE before any other script
# ============================================================

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# --- Scripts 01-02: Data and DE Analysis ---
BiocManager::install(c(
  "GEOquery", "DESeq2", "limma", "sva", "edgeR",
  "EnhancedVolcano", "org.Hs.eg.db", "Biobase",
  "BiocGenerics", "Biostrings", "genefilter"
), ask = FALSE)

install.packages(c("tidyverse", "pheatmap", "ggplot2", "ggrepel"))

# --- Script 03: Enrichment Analysis ---
BiocManager::install(c(
  "GO.db", "clusterProfiler", "enrichplot"
), ask = FALSE)

# --- Script 04: WGCNA ---
BiocManager::install(c(
  "WGCNA", "impute", "preprocessCore", "dynamicTreeCut"
), ask = FALSE)

# --- Script 05: Regulatory Network ---
BiocManager::install(c("multiMiR"), ask = FALSE)

# dorothea and viper needs special installation
BiocManager::install(c("dorothea", "viper"), ask = FALSE)

BiocManager::install("decoupleR", ask = FALSE)

# --- Script 06: Validation ---
install.packages(c("VennDiagram", "ggvenn"))


### verification of all installation 
required_packages <- c(
  "GEOquery", "DESeq2", "limma", "sva", "edgeR",
  "EnhancedVolcano", "org.Hs.eg.db", "Biobase",
  "BiocGenerics", "Biostrings", "genefilter",
  "tidyverse", "pheatmap", "ggplot2", "ggrepel",
  "GO.db", "clusterProfiler", "enrichplot", "ReactomePA",
  "WGCNA", "impute", "preprocessCore", "dynamicTreeCut",
  "multiMiR", "igraph", "ggraph", "VennDiagram", "ggvenn"
)

missing <- required_packages[!required_packages %in% 
                               installed.packages()[,"Package"]]

if (length(missing) == 0) {
  message("All packages installed successfully.")
} else {
  message("Missing packages: ", paste(missing, collapse = ", "))
}











