# ============================================================
# Script 06: Cross-Validation Analysis
# ============================================================

library(GEOquery)
library(limma)
library(edgeR)
library(tidyverse)
library(VennDiagram)
library(org.Hs.eg.db)
library(pheatmap)

dir_figures <- "results/figures/"
dir_tables  <- "results/tables/"

# Load hub genes from Script 04
hub_df <- readRDS("data/processed/hub_genes.rds")
degs   <- readRDS("data/processed/degs.rds")

# Step 1: Download validation dataset GSE141910
gse141910 <- getGEO("GSE141910", GSEMatrix = TRUE, getGPL = FALSE)
eset_val   <- gse141910[[1]]

# Check metadata
dim(eset_val)
head(pData(eset_val)[, c("title", "characteristics_ch1")])

# Step 2: Download supplementary files
getGEOSuppFiles("GSE141910", makeDirectory = TRUE,
                baseDir = "data/raw/")

# Step 3: Extract the tar archive
untar("data/raw/GSE141910/GSE141910_RAW.tar",
      exdir = "data/raw/GSE141910/")

# List extracted files
files <- list.files("data/raw/GSE141910/", 
                    pattern = "*.gz",
                    full.names = TRUE)
length(files)
head(files)

# Step 4: Get sample metadata and filter DCM + NF only
metadata_val <- pData(eset_val)[, c("title", "characteristics_ch1",
                                    "characteristics_ch1.1")]

head(metadata_val, 10)
table(metadata_val$characteristics_ch1)

# Step 5: Filter to DCM and Non-Failing only
metadata_val$group <- ifelse(
  grepl("Non-Failing", metadata_val$characteristics_ch1), "NF",
  ifelse(grepl("Dilated cardiomyopathy", 
               metadata_val$characteristics_ch1), "DCM", NA)
)

metadata_val <- metadata_val[!is.na(metadata_val$group), ]
table(metadata_val$group)

# Match files to filtered samples
sample_ids <- metadata_val$title
keep_files <- files[grepl(paste(sample_ids, collapse="|"), files)]
length(keep_files)

# Step 6: Load and merge count files
# Read first file to check structure
test <- read.csv(keep_files[1], header = TRUE)
head(test)
dim(test)

# Step 7: Load all files and merge
count_list <- lapply(keep_files, function(f) {
  df <- read.csv(f, header = TRUE)
  colnames(df) <- c("gene", basename(f))
  df
})

# Merge all by gene column
counts_val <- Reduce(function(x, y) 
  merge(x, y, by = "gene", all = FALSE), count_list)

# Set gene names as rownames
rownames(counts_val) <- counts_val$gene
counts_val <- counts_val[, -1]

dim(counts_val)

# Step 8: Create sample metadata with correct group assignment
# Matching filenames back to metadata using patient IDs
file_ids <- gsub(".*GSM[0-9]+_(.+)\\.csv\\.gz", "\\1",
                 basename(keep_files))

sample_info_val <- data.frame(
  sample  = colnames(counts_val),
  file_id = file_ids,
  group   = metadata_val$group[match(file_ids, metadata_val$title)],
  row.names = colnames(counts_val)
)

sample_info_val$group <- factor(sample_info_val$group,
                                levels = c("NF", "DCM"))
table(sample_info_val$group)

# Step 9: Differential expression analysis on validation dataset
# Using limma directly since data is already log-transformed
design_val <- model.matrix(~ group, data = sample_info_val)
fit_val    <- lmFit(as.matrix(counts_val), design_val)
fit_val    <- eBayes(fit_val)

deg_val <- topTable(fit_val,
                    coef    = "groupDCM",
                    number  = Inf,
                    sort.by = "P")

degs_val <- deg_val %>%
  filter(adj.P.Val < 0.05, abs(logFC) > 1.5)

nrow(degs_val)

# Step 10: Validate hub genes using relaxed threshold
# Hub genes identified by WGCNA connectivity may show smaller but consistent fold changes
# directional consistency is key
hub_in_val <- deg_val[rownames(deg_val) %in% hub_df$ensembl, ]
hub_in_val$symbol <- hub_df$symbol[match(rownames(hub_in_val),
                                         hub_df$ensembl)]

# Show validation results for all hub genes
hub_in_val[, c("symbol", "logFC", "P.Value", "adj.P.Val")]

# Save validated hub genes — those significant in validation
validated_hubs <- hub_in_val[hub_in_val$adj.P.Val < 0.05, ]
write.csv(validated_hubs,
          paste0(dir_tables, "validated_hub_genes.csv"))

# Step 11: DEG overlap between primary and validation datasets
primary_degs    <- rownames(degs)
validation_degs <- rownames(degs_val)
overlap         <- intersect(primary_degs, validation_degs)

cat("Primary DEGs:", length(primary_degs), "\n")
cat("Validation DEGs:", length(validation_degs), "\n")
cat("Overlapping DEGs:", length(overlap), "\n")
cat("Overlap percentage:", 
    round(length(overlap)/length(primary_degs)*100, 1), "%\n")

# Step 12: Venn diagram of DEG overlap
venn.diagram(
  x = list(Primary    = primary_degs,
           Validation = validation_degs),
  category.names = c("Primary\n(GSE116250)",
                     "Validation\n(GSE141910)"),
  filename = paste0(dir_figures, "venn_diagram.png"),
  output   = TRUE,
  fill     = c("#3498DB", "#E74C3C"),
  alpha    = 0.5,
  cex      = 1.5,
  cat.cex  = 1.2,
  main     = "DEG Overlap: Primary vs Validation"
)

#Save all validation results
write.csv(as.data.frame(deg_val),
          paste0(dir_tables, "validation_DEG_results.csv"))
saveRDS(deg_val, "data/processed/validation_deg_results.rds")





### Pathway-level validation
# Get entrez IDs for validation DEGs
val_entrez <- mapIds(org.Hs.eg.db,
                     keys      = rownames(degs_val),
                     column    = "ENTREZID",
                     keytype   = "ENSEMBL",
                     multiVals = "first")

val_entrez <- val_entrez[!is.na(val_entrez)]

# GO enrichment on validation DEGs
library(clusterProfiler)
go_val <- enrichGO(gene          = val_entrez,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   readable      = TRUE)

# Compare top terms with primary enrichment
head(go_val@result[, c("Description", "p.adjust")], 15)

# Save pathway validation results
write.csv(go_val@result,
          paste0(dir_tables, "validation_GO_BP_results.csv"))

dotplot(go_val, showCategory = 15,
        title = "GO BP Enrichment - Validation Dataset (GSE141910)")

ggsave(paste0(dir_figures, "validation_GO_dotplot.png"),
       width = 10, height = 8)

saveRDS(go_val, "data/processed/validation_go_results.rds")


