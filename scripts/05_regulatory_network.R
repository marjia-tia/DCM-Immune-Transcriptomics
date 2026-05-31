# ============================================================
# Script 05: Regulatory Network Analysis 
#              (miRNA-mRNA Interactions + Drug Repurposing)
# ============================================================

library(multiMiR)
library(igraph)
library(ggraph)
library(tidyverse)
library(org.Hs.eg.db)

dir_figures <- "results/figures/"
dir_tables  <- "results/tables/"

# Load hub genes from Script 04
hub_df <- readRDS("data/processed/hub_genes.rds")
degs   <- readRDS("data/processed/degs.rds")

head(hub_df)


## Step 1: Query miRNA-mRNA interactions for hub genes
# Remove NA symbols first
hub_symbols <- hub_df$symbol[!is.na(hub_df$symbol)]
hub_symbols <- hub_symbols[hub_symbols != "<NA>"]

# Query multiMiR
mirna_results <- get_multimir(
  org        = "hsa",
  target     = hub_symbols,
  table      = "validated",
  summary    = TRUE
)

# Extract results
mirna_df <- mirna_results@data
dim(mirna_df)
head(mirna_df)


## Step 2: Build miRNA-mRNA network
# Get top miRNAs by frequency
mirna_freq <- mirna_df %>%
  group_by(mature_mirna_id) %>%
  summarise(n_targets = n_distinct(target_symbol)) %>%
  arrange(desc(n_targets))

# Keep top 20 miRNAs
top_mirnas <- head(mirna_freq$mature_mirna_id, 20)

# Filter interactions to top miRNAs
net_df <- mirna_df %>%
  filter(mature_mirna_id %in% top_mirnas) %>%
  dplyr::select(mature_mirna_id, target_symbol) %>%
  distinct()

# Build igraph network
g <- graph_from_data_frame(net_df, directed = TRUE)

# Add node type
V(g)$type <- ifelse(V(g)$name %in% top_mirnas, "miRNA", "Hub Gene")

# Plot
ggraph(g, layout = "fr") +
  geom_edge_arc(alpha = 0.3, color = "grey50") +
  geom_node_point(aes(color = type), size = 4) +
  geom_node_text(aes(label = name), size = 2.5, repel = TRUE) +
  scale_color_manual(values = c("miRNA" = "#E74C3C", 
                                "Hub Gene" = "#3498DB")) +
  theme_void() +
  labs(title = "miRNA-Hub Gene Regulatory Network",
       color = "Node Type")

ggsave(paste0(dir_figures, "miRNA_network.png"), 
       width = 14, height = 10)


## Step 3: Drug repurposing via DGIdb API
# Save miRNA results
write.csv(mirna_df, paste0(dir_tables, "miRNA_interactions.csv"),
          row.names = FALSE)

# Top miRNAs by target count
write.csv(mirna_freq, paste0(dir_tables, "miRNA_frequency.csv"),
          row.names = FALSE)

# Drug repurposing - query DGIdb for hub genes
hub_symbols_clean <- hub_df$symbol[!is.na(hub_df$symbol) & 
                                     hub_df$symbol != "<NA>"]

# Drug Repurposing — DGIdb
# Results obtained manually from https://www.dgidb.org
# Query genes: all hub genes from blue module
# Results downloaded and saved as CSV

drug_results <- read.csv("results/tables/drug_repurposing_DGIdb.csv", 
                         sep = "\t", header = TRUE)
head(drug_results)
colnames(drug_results)


saveRDS(mirna_df, "data/processed/mirna_interactions.rds")
saveRDS(drug_results, "data/processed/drug_results.rds")

  

  
  



