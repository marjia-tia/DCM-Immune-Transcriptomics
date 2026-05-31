# ============================================================
# Script 02: Differential Expression Analysis
# ============================================================


v           <- readRDS("data/processed/voom_object.rds")
sample_info <- readRDS("data/processed/sample_info.rds")

library(limma)
library(edgeR)
library(EnhancedVolcano)
library(pheatmap)
library(tidyverse)

dir_figures <- "results/figures/"
dir_tables  <- "results/tables/"

# --- Fit linear model ---
design <- model.matrix(~ group, data = sample_info)
fit    <- lmFit(v, design)
fit    <- eBayes(fit)

# --- Extract results ---
deg_results <- topTable(fit, 
                        coef     = "groupDCM",
                        number   = Inf,
                        sort.by  = "P")

head(deg_results)
nrow(deg_results)
colnames(deg_results)

# Filter significant DEGs
degs <- deg_results %>%
  filter(adj.P.Val < 0.05, abs(logFC) > 1.5)

nrow(degs)
table(degs$logFC > 0)

# Volcano plot
EnhancedVolcano(deg_results,
                lab     = rownames(deg_results),
                x       = "logFC",
                y       = "adj.P.Val",
                title   = "DCM vs Non-Failing: Differential Expression",
                pCutoff = 0.05,
                FCcutoff = 1.5,
                pointSize = 2.0,
                labSize = 3.0)

ggsave(paste0(dir_figures, "volcano_plot.png"), width = 10, height = 8)

# Save DEG table
write.csv(degs, file = paste0(dir_tables, "DEG_results.csv"), 
          row.names = TRUE)

# Save objects for Script 03
saveRDS(degs, "data/processed/degs.rds")
saveRDS(deg_results, "data/processed/deg_results_full.rds")



# Add gene symbols to DEG results
library(org.Hs.eg.db)

gene_symbols <- mapIds(org.Hs.eg.db,
                       keys      = rownames(degs),
                       column    = "SYMBOL",
                       keytype   = "ENSEMBL",
                       multiVals = "first")

degs$gene_symbol <- gene_symbols
head(degs[, c("gene_symbol", "logFC", "adj.P.Val")])

# Add symbols to full results for volcano
deg_results$gene_symbol <- mapIds(org.Hs.eg.db,
                                  keys     = rownames(deg_results),
                                  column   = "SYMBOL", 
                                  keytype  = "ENSEMBL",
                                  multiVals = "first")

EnhancedVolcano(deg_results,
                lab      = deg_results$gene_symbol,
                x        = "logFC",
                y        = "adj.P.Val",
                title    = "DCM vs Non-Failing: Differential Expression",
                pCutoff  = 0.05,
                FCcutoff = 1.5,
                pointSize = 2.0,
                labSize  = 3.5,
                drawConnectors = TRUE)

ggsave(paste0(dir_figures, "volcano_plot_labeled.png"), 
       width = 12, height = 9)

# Heatmap of top 50 DEGs
top50 <- head(degs[order(degs$adj.P.Val), ], 50)
mat   <- v$E[rownames(top50), ]
rownames(mat) <- top50$gene_symbol

# Annotation for samples
annotation_col <- data.frame(
  Group = sample_info$group,
  row.names = colnames(mat)
)

pheatmap(mat,
         annotation_col  = annotation_col,
         scale           = "row",
         show_colnames   = FALSE,
         fontsize_row    = 7,
         main            = "Top 50 DEGs: DCM vs Non-Failing",
         filename        = paste0(dir_figures, "heatmap_top50.png"))






