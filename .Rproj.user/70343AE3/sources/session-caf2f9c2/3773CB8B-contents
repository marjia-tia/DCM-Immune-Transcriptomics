# ============================================================
# Script 03: Functional Enrichment Analysis
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(tidyverse)

dir_figures <- "results/figures/"
dir_tables  <- "results/tables/"

# Load DEGs from Script 02
degs        <- readRDS("data/processed/degs.rds")
deg_results <- readRDS("data/processed/deg_results_full.rds")


## --- Step 1: Convert Ensembl IDs to Entrez IDs ---
entrez_ids <- mapIds(org.Hs.eg.db,
                     keys      = rownames(degs),
                     column    = "ENTREZID",
                     keytype   = "ENSEMBL",
                     multiVals = "first")

entrez_ids <- entrez_ids[!is.na(entrez_ids)]

## --- Step 2: GO Enrichment Analysis ---
go_bp <- enrichGO(gene         = entrez_ids,
                       OrgDb        = org.Hs.eg.db,
                       ont          = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       readable      = TRUE)

# GO Molecular Function
go_mf <- enrichGO(gene          = entrez_ids,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "MF",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  readable      = TRUE)

# GO Cellular Component  
go_cc <- enrichGO(gene          = entrez_ids,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "CC",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  readable      = TRUE)

# Plot _ BP, MF, CC
dotplot(go_results, showCategory = 20, 
        title = "GO Biological Process - DCM DEGs")

ggsave(paste0(dir_figures, "GO_BP_dotplot.png"), 
       width = 10, height = 8)

dotplot(go_mf, showCategory = 20,
        title = "GO Molecular Function - DCM DEGs")
ggsave(paste0(dir_figures, "GO_MF_dotplot.png"), width = 10, height = 8)

dotplot(go_cc, showCategory = 20,
        title = "GO Cellular Component - DCM DEGs")
ggsave(paste0(dir_figures, "GO_CC_dotplot.png"), width = 10, height = 8)

# Save all results
write.csv(as.data.frame(go_bp), 
          paste0(dir_tables, "GO_BP_results.csv"))
write.csv(as.data.frame(go_mf), 
          paste0(dir_tables, "GO_MF_results.csv"))
write.csv(as.data.frame(go_cc), 
          paste0(dir_tables, "GO_CC_results.csv"))

saveRDS(go_bp, "data/processed/go_bp.rds")
saveRDS(go_mf, "data/processed/go_mf.rds")
saveRDS(go_cc, "data/processed/go_cc.rds")



