# ============================================================
# Script 01: Data Acquisition and Preprocessing
# ============================================================

# --- 1. Load Libraries ---
library(GEOquery)
library(tidyverse)
library(DESeq2)
library(limma)
library(sva)        # for ComBat batch correction
library(pheatmap)
library(ggplot2)

# --- 2. Set Directories ---
dir_raw       <- "data/raw/"
dir_processed <- "data/processed/"
dir_figures   <- "results/figures/"

# --- 3. Download GSE116250 (Primary Discovery Dataset) ---
gse116250 <- getGEO("GSE116250", GSEMatrix = TRUE, getGPL = FALSE)

# GSE116250 contains multiple expression sets - extract the first one
eset <- gse116250[[1]]

# Look at what we downloaded
dim(eset)           # rows = genes, columns = samples
head(pData(eset))   # sample metadata

counts_raw <- read.delim(
  "data/raw/GSE116250/GSE116250_rpkm.txt.gz",
  header = TRUE,
  row.names = NULL
)

dim(counts_raw)
head(counts_raw[, 1:5])

# Check column names to see all sample IDs
colnames(counts_raw)

# Check how many duplicated gene names exist
sum(duplicated(counts_raw$Gene))

# Step A: Remove duplicate genes by keeping the one with highest mean expression
counts_raw$mean_expr <- rowMeans(counts_raw[, 3:ncol(counts_raw)])
counts_clean <- counts_raw %>%
  group_by(Gene) %>%
  slice_max(mean_expr, n = 1, with_ties = FALSE) %>%
  ungroup()

# Step B: Set gene names as row names
counts_clean <- as.data.frame(counts_clean)
rownames(counts_clean) <- counts_clean$Gene
counts_clean <- counts_clean %>% 
  select(-Gene, -Common_name, -mean_expr)

# Step C: Keep only DCM and NF samples
dcm_cols <- grep("^DCM", colnames(counts_clean), value = TRUE)
nf_cols  <- grep("^NF", colnames(counts_clean), value = TRUE)
counts_filtered <- counts_clean[, c(nf_cols, dcm_cols)]

dim(counts_filtered)

# Keep only DCM and NF samples, remove ICM
table(gsub("[0-9]+", "", colnames(counts_filtered)))

# Step D: Filter low-expression genes
# Keep genes with RPKM > 1 in at least 50% of samples
keep <- rowSums(counts_filtered > 1) >= (ncol(counts_filtered) * 0.5)
counts_filtered <- counts_filtered[keep, ]

dim(counts_filtered)

# Step E: Create sample metadata
sample_info <- data.frame(
  sample = colnames(counts_filtered),
  group  = gsub("[0-9]+", "", colnames(counts_filtered)),
  row.names = colnames(counts_filtered)
)

sample_info$group <- factor(sample_info$group, levels = c("NF", "DCM"))
table(sample_info$group)

# Step F: Create DGEList and normalize with voom
library(edgeR)

dge <- DGEList(counts = counts_filtered, 
               group = sample_info$group)
dge <- calcNormFactors(dge, method = "TMM")

# Apply voom transformation
design <- model.matrix(~ group, data = sample_info)
v <- voom(dge, design, plot = TRUE)

# Step G: PCA plot
pca <- prcomp(t(v$E), scale. = TRUE)

pca_df <- data.frame(
  PC1   = pca$x[,1],
  PC2   = pca$x[,2],
  group = sample_info$group
)

ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "PCA: DCM vs Non-Failing Samples",
       x = paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 1), "%)"))

ggsave(paste0(dir_figures, "PCA_plot.png"), width = 7, height = 5)

# Step H: Save processed data
saveRDS(v, file = "data/processed/voom_object.rds")
saveRDS(sample_info, file = "data/processed/sample_info.rds")
saveRDS(counts_filtered, file = "data/processed/counts_filtered.rds")



