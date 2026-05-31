# ============================================================
# Script 04: WGCNA Co-expression Network Analysis
# ============================================================

library(WGCNA)
library(tidyverse)
library(pheatmap)

options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# Load data
v           <- readRDS("data/processed/voom_object.rds")
sample_info <- readRDS("data/processed/sample_info.rds")

# Prepare expression matrix
expr_mat <- t(v$E)
dim(expr_mat)


## Step 1: Choose soft thresholding power
powers <- c(1:20)

sft <- pickSoftThreshold(expr_mat,
                         powerVector  = powers,
                         verbose      = 5,
                         networkType  = "signed")

# Plot scale-free topology fit
par(mfrow = c(1,2))
plot(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab  = "Soft Threshold (power)",
     ylab  = "Scale Free Topology Model Fit (R^2)",
     main  = "Scale independence")
abline(h = 0.80, col = "red")

plot(sft$fitIndices[,1],
     sft$fitIndices[,5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     main = "Mean connectivity")


## Step 2: Build network and identify modules
net <- blockwiseModules(expr_mat,
                        power             = 9,
                        networkType       = "signed",
                        TOMType           = "signed",
                        minModuleSize     = 30,
                        mergeCutHeight    = 0.25,
                        numericLabels     = FALSE,
                        verbose           = 3)

table(net$colors)


## Step 3: Module-trait correlation
# Create trait matrix
trait <- as.data.frame(as.numeric(sample_info$group == "DCM"))
rownames(trait) <- rownames(expr_mat)
colnames(trait) <- "DCM_status"

# Calculate module eigengenes
MEs <- moduleEigengenes(expr_mat, net$colors)$eigengenes
MEs <- orderMEs(MEs)

# Correlate with trait
moduleTraitCor  <- cor(MEs, trait, use = "p")
moduleTraitPval <- corPvalueStudent(moduleTraitCor, nrow(expr_mat))

# Plot heatmap
png(paste0("results/figures/module_trait_heatmap.png"), 
    width = 800, height = 1200)

par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix     = moduleTraitCor,
               xLabels    = "DCM Status",
               yLabels    = names(MEs),
               colorLabels = FALSE,
               colors     = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text   = 0.5,
               main       = "Module-Trait Correlation")

dev.off()

# Find top positively and negatively correlated modules
module_cor_df <- data.frame(
  module     = rownames(moduleTraitCor),
  correlation = moduleTraitCor[,1],
  pvalue     = moduleTraitPval[,1]
)

# Sort by correlation
module_cor_df <- module_cor_df[order(module_cor_df$correlation, 
                                     decreasing = TRUE), ]

# Top positive
head(module_cor_df, 5)

# Top negative
tail(module_cor_df, 5)


# Step 4: Extract hub genes from blue module
blue_genes <- names(net$colors[net$colors == "blue"])

# Calculate module membership
kME <- signedKME(expr_mat, MEs)

# Get blue module membership scores
blue_kME <- kME[blue_genes, "kMEblue"]
names(blue_kME) <- blue_genes
blue_kME <- sort(blue_kME, decreasing = TRUE)

# Top 20 hub genes
hub_genes <- names(head(blue_kME, 20))
hub_genes

# Convert hub genes to symbols
library(org.Hs.eg.db)

hub_symbols <- mapIds(org.Hs.eg.db,
                      keys      = hub_genes,
                      column    = "SYMBOL",
                      keytype   = "ENSEMBL",
                      multiVals = "first")

hub_df <- data.frame(
  ensembl    = hub_genes,
  symbol     = hub_symbols,
  kME        = blue_kME[hub_genes]
)

print(hub_df)

# Save
write.csv(hub_df, "results/tables/hub_genes_blue_module.csv",
          row.names = FALSE)

saveRDS(hub_df, "data/processed/hub_genes.rds")
saveRDS(net, "data/processed/wgcna_net.rds")
saveRDS(MEs, "data/processed/module_eigengenes.rds")









