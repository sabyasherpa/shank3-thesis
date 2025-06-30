# This script is to analyse GSE201621 data with DESeq2. 

#install required packages and load them
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("ggplot2")
BiocManager::install("tidyr")
BiocManager::install("dplyr")
BiocManager::install("stringr")
BiocManager::install("ggfortify")
BiocManager::install("ggrepel")
BiocManager::install("pheatmap")
BiocManager::install("orthogene")
BiocManager::install("bioMart")
BiocManager::install("remotes")
BiocManager::install("vitkl/orthologsBioMART")

library(DESeq2)
library(ggplot2)
library(ggfortify)
library(tidyr)
library(dplyr)
library(stringr)
library(ggrepel)
library(pheatmap)
library(orthogene)
library(gprofiler2)
library(biomaRt)
library(orthologsBioMART)

#load processed count data and create some overview plots.
data <- read.table("GSE201621_processed_data_DESeq_count.txt", header = TRUE, sep = "\t", check.names = FALSE)

#Some raw count data visualisation with boxplot and PCA
# Add a small pseudocount (like 1) to avoid log(0)
boxplot(log(data[ , -1] + 1))
PCA <- prcomp(t(log(data[ , -1] + 1)), scale. = TRUE)
autoplot(PCA, label = TRUE, label.size = 3)


# Load metadata file
metadata <- read.csv("SraRunTable (1).csv", header = TRUE, check.names = FALSE)

# Reverse metadata order (because it's flipped vs. count data)
metadata <- metadata[rev(seq_len(nrow(metadata))), ]

# assign the column names from the count matrix in the correct order
metadata$sample_name <- colnames(data)[-1]

#this is analysis of prefrontal cortex

# Subset metadata for prefrontal cortex only
metadata_ctx <- metadata[metadata$tissue == "prefrontal cortex", ]

# Subset the columns: keep gene IDs + cortex sample names
data_ctx <- data[ , c("ensembl_gene_id", metadata_ctx$sample_name)]

# Set gene IDs as row names
rownames(data_ctx) <- data_ctx$ensembl_gene_id
data_ctx <- data_ctx[ , -1]  # remove gene ID column

# Set sample names as rownames in metadata
rownames(metadata_ctx) <- metadata_ctx$sample_name

#make all data round numbers
data_ctx <- round(data_ctx)

# Create DESeq2 object
dds_ctx <- DESeqDataSetFromMatrix(
  countData = data_ctx,
  colData = metadata_ctx,
  design = ~ genotype  # only comparing HT vs WT
)
#level the deseq2
metadata_ctx$genotype <- factor(metadata_ctx$genotype, levels = c("wildtype", "heterozygous"))

# Run DESeq2
dds_ctx <- DESeq(dds_ctx)

# Extract results: heterozygous vs wildtype
res_ctx <- results(dds_ctx, contrast = c("genotype", "heterozygous", "wildtype"))
#Sort by adjusted p-value (most significant genes first)
res_ctx <- res_ctx[order(res_ctx$padj), ]  # sort by adjusted p-value

# View the top differentially expressed genes
head(res_ctx)

# Save all DESeq2 results
# Filter significant DEGs (padj < 0.05)
SigDEG <- subset(res_ctx, padj < 0.05 & !is.na(padj))
write.csv(as.data.frame(SigDEG), "DESeq2_HT_vs_WT_Cortex_SigDEG.csv")

relaxed_DEGs_ctx_test <- subset(res_ctx, 
                           padj < 0.05 & abs(log2FoldChange) > 0.5 & !is.na(padj) & !is.na(log2FoldChange))

# Filter DEGs using relaxed thresholds
relaxed_DEGs_ctx <- subset(res_ctx, 
                       pvalue < 0.05 & abs(log2FoldChange) > 0.5 & !is.na(pvalue) & !is.na(log2FoldChange))

# Save to CSV
write.csv(as.data.frame(relaxed_DEGs_ctx), "cortex_HT_vs_WT_Relaxed_DEGs.csv")

# Count how many
nrow(relaxed_DEGs_ctx)
# Upregulated
up_DEGs_ctx <- subset(relaxed_DEGs_ctx, log2FoldChange > 0.5)
nrow(up_DEGs_ctx)

# Downregulated
down_DEGs_ctx <- subset(relaxed_DEGs_ctx, log2FoldChange < -0.5)
nrow(down_DEGs_ctx)

#get the full list of genes degs of cortex for pathway overlay analysis
#full list dds_ctx
full_list_ctx <- results(dds_ctx)

# Convert to data frame and add gene IDs as a column
full_list_ctx_df <- as.data.frame(full_list_ctx)
full_list_ctx_df$Mouse_Ensembl <- rownames(full_list_ctx_df)

# Replace NA values with blank (optional but helps Excel)
full_list_ctx_df[is.na(full_list_ctx_df)] <- ""

# Export as CSV for BioMart mapping (safe for Excel)
write.csv(full_list_ctx_df, "full_list_clean_for_mapping.csv", row.names = FALSE, na = "")

# Run full mapping:
mapped <- findOrthologsMmHs(from_filters = "ensembl_gene_id",
                            from_values = full_list_ctx_df$Mouse_Ensembl,
                            to_attributes = "hgnc_symbol")
# Merge with full DESeq2 output
final_merge <- merge(full_list_df, mapped, by.x = "Mouse_Ensembl", by.y = "mouse_ensembl_gene_id")

# Export the final full mapped file
write.csv(final_merge, "full_DESeq2_mapped_cortex.csv", row.names = FALSE)
write.csv(final_merge, "full_DESeq2_mapped_cortex_clean.csv", 
          row.names = FALSE, na = "", quote = FALSE)

#this is analysis of the hippocampus

#  Subset metadata: hippocampus samples, exclude KI
metadata_hpc <- metadata[
  metadata$tissue == "hippocampus" &
    metadata$genotype %in% c("wildtype", "heterozygous"),
]

# Set factor levels (WT = control)
metadata_hpc$genotype <- factor(metadata_hpc$genotype, levels = c("wildtype", "heterozygous"))

#  Subset count matrix for these samples
data_hpc <- data[ , c("ensembl_gene_id", metadata_hpc$sample_name)]
data_hpc <- as.data.frame(data_hpc)


# Set gene names as rownames, and remove gene ID column
rownames(data_hpc) <- data_hpc$ensembl_gene_id
data_hpc <- data_hpc[ , -1]

#  Match row names in metadata
rownames(metadata_hpc) <- metadata_hpc$sample_name

#make all data round numbers
data_ctx <- round(data_ctx)

#  Run DESeq2
dds_hpc <- DESeqDataSetFromMatrix(countData = round(data_hpc), colData = metadata_hpc, design = ~ genotype)
dds_hpc <- DESeq(dds_hpc)

#  Get results (HT vs WT)
res_hpc <- results(dds_hpc, contrast = c("genotype", "heterozygous", "wildtype"))
res_hpc <- res_hpc[order(res_hpc$padj), ]

#  Save all results
write.csv(as.data.frame(res_hpc), "DESeq2_HT_vs_WT_Hippocampus_FULL.csv")

# Filter DEGs using relaxed thresholds
relaxed_DEGs_hpc <- subset(res_hpc, 
                       pvalue < 0.05 & abs(log2FoldChange) > 0.5 & !is.na(pvalue) & !is.na(log2FoldChange))

# Save to CSV
write.csv(as.data.frame(relaxed_DEGs_hpc), "hippocampus_HT_vs_WT_Relaxed_DEGs.csv")

# Count how many
nrow(relaxed_DEGs_hpc)

# Upregulated
up_DEGs_hpc <- subset(relaxed_DEGs_hpc, log2FoldChange > 0.5)
nrow(up_DEGs_hpc)

# Downregulated
down_DEGs_hpc <- subset(relaxed_DEGs_hpc, log2FoldChange < -0.5)
nrow(down_DEGs_hpc)

#get the full list of genes degs of hippocampus for pathway overlay analysis
full_list_hpc <- results(dds_hpc)
full_list_hpc_df <- as.data.frame(full_list_hpc)
full_list_hpc_df$Mouse_Ensembl <- rownames(full_list_hpc_df)

# Perform full ortholog mapping using orthologsBioMART
mapped_hpc <- findOrthologsMmHs(
  from_filters = "ensembl_gene_id",
  from_values = full_list_hpc_df$Mouse_Ensembl,
  to_attributes = "hgnc_symbol"
)

# Merge the DESeq2 result with the orthology mapping
final_merge_hpc <- merge(full_list_hpc_df, mapped_hpc, by.x = "Mouse_Ensembl", by.y = "mouse_ensembl_gene_id")

# Clean up HGNC symbols (remove quotes problem)
final_merge_hpc$hgnc_symbol <- toupper(final_merge_hpc$hgnc_symbol)

# Export PathVisio-safe CSV
write.csv(final_merge_hpc, "full_DESeq2_mapped_hippocampus.csv", row.names = FALSE, na = "", quote = FALSE)

# this is analysis of striatum

# Subset metadata for striatum samples, exclude KI
metadata_str <- metadata[
  metadata$tissue == "striatum" & metadata$genotype %in% c("wildtype", "heterozygous"),
]

#  Set genotype factor levels (WT = reference)
metadata_str$genotype <- factor(metadata_str$genotype, levels = c("wildtype", "heterozygous"))

# Subset count matrix to matching striatum samples
data_str <- data[ , c("ensembl_gene_id", metadata_str$sample_name)]
data_str <- as.data.frame(data_str)

# Set gene IDs as rownames and remove gene ID column
rownames(data_str) <- data_str$ensembl_gene_id
data_str <- data_str[ , -1]
rownames(metadata_str) <- metadata_str$sample_name

#  Run DESeq2
dds_str <- DESeqDataSetFromMatrix(countData = round(data_str), colData = metadata_str, design = ~ genotype)
dds_str <- DESeq(dds_str)

#  Get results
res_str <- results(dds_str, contrast = c("genotype", "heterozygous", "wildtype"))
res_str <- res_str[order(res_str$padj), ]

 Save results
write.csv(as.data.frame(res_str), "DESeq2_HT_vs_WT_Striatum_FULL.csv")

# Filter DEGs using relaxed thresholds
relaxed_DEGs <- subset(res_str, 
                       pvalue < 0.05 & abs(log2FoldChange) > 0.5 & !is.na(pvalue) & !is.na(log2FoldChange))

# Save to CSV
write.csv(as.data.frame(relaxed_DEGs), "Striatum_HT_vs_WT_Relaxed_DEGs.csv")

# Count how many
nrow(relaxed_DEGs)

# Total DEGs
relaxed_DEGs <- subset(res_str, pvalue < 0.05 & abs(log2FoldChange) > 0.5 & !is.na(pvalue))

# Upregulated
up_DEGs <- subset(relaxed_DEGs, log2FoldChange > 0.5)
n_up <- nrow(up_DEGs)
n_up

# Downregulated
down_DEGs <- subset(relaxed_DEGs, log2FoldChange < -0.5)
n_down <- nrow(down_DEGs)
n_down
# Print results
cat("Striatum - Relaxed threshold:\n")
cat("Upregulated genes:", n_up, "\n")
cat("Downregulated genes:", n_down, "\n")
cat("Total DEGs:", n_up + n_down, "\n")
nrow(up_DEGs) + nrow(down_DEGs) == nrow(relaxed_DEGs)

#get the full list of genes degs of striatum for pathway overlay analysis
full_list_str <- results(dds_str)
full_list_str_df <- as.data.frame(full_list_str)
full_list_str_df$Mouse_Ensembl <- rownames(full_list_str_df)

# Perform full ortholog mapping
mapped_str <- findOrthologsMmHs(
  from_filters = "ensembl_gene_id",
  from_values = full_list_str_df$Mouse_Ensembl,
  to_attributes = "hgnc_symbol"
)

# Merge the DESeq2 result with the orthology mapping
final_merge_str <- merge(full_list_str_df, mapped_str, by.x = "Mouse_Ensembl", by.y = "mouse_ensembl_gene_id")

# Clean up HGNC symbols (remove quotes problem)
final_merge_str$hgnc_symbol <- toupper(final_merge_str$hgnc_symbol)

# Export PathVisio-safe CSV
write.csv(final_merge_str, "full_DESeq2_mapped_striatum.csv", row.names = FALSE, na = "", quote = FALSE)



#orthology mapping of cortex and go enrichment
# Load orthologs correctly with proper column names
orthologs <- read.csv("mart_export cortex degs.txt", header = TRUE, stringsAsFactors = FALSE)
colnames(orthologs) <- c("Gene_stable_ID", "Human_gene_stable_ID", "Human_gene_name")

# Add Gene ID as a column to your DEGs
relaxed_DEGs_ctx$Gene_stable_ID <- rownames(relaxed_DEGs_ctx)


# Merge DEGs with orthologs using mouse Gene ID
mapped_DEGs_ctx <- merge(relaxed_DEGs_ctx, orthologs, by = "Gene_stable_ID")

# Keep only DEGs with a human gene symbol (i.e., successful mapping)
mapped_DEGs_ctx <- mapped_DEGs_ctx[!is.na(mapped_DEGs_ctx$Human_gene_name) & mapped_DEGs_ctx$Human_gene_name != "", ]
# Merge using only common genes (inner join behavior)
mapped_DEGs_ctx <- merge(
  x = relaxed_DEGs_ctx,
  y = orthologs,
  by = "Gene_stable_ID"
)



# cortex mapped DEGs
gost_ctx <- gost(mapped_DEGs_ctx$Human_gene_name, organism = "hsapiens", significant = TRUE)
nrow(gost_ctx$result)


# See top terms
head(gost_ctx$result[, c("term_name", "intersection_size", "p_value")])



# Extract relevant info
go_df <- gost_ctx$result
go_df <- go_df[order(go_df$p_value), ][1:10, ]  # top 10 terms
go_df[, c("term_name", "intersection_size", "p_value")]

# Plot
ggplot(go_df, aes(x = reorder(term_name, -log10(p_value)), 
                  y = -log10(p_value), 
                  size = intersection_size)) +
  geom_point(aes(color = -log10(p_value))) +
  coord_flip() +
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  labs(
    title = " ",
    x = "GO Biological Process",
    y = "-log10(p-value)",
    size = "Gene Count"
  ) +
  theme_minimal()





#orthology mapping of hippocampus and go enrichment
# Load orthologs correctly with proper column names
orthologs <- read.csv("mart_export hippocampus.txt", header = TRUE, stringsAsFactors = FALSE)
colnames(orthologs) <- c("Gene_stable_ID", "Human_gene_stable_ID", "Human_gene_name")

# Add Gene ID as a column to your DEGs
relaxed_DEGs_hpc$Gene_stable_ID <- rownames(relaxed_DEGs_hpc)
# Merge DEGs with orthologs using mouse Gene ID
mapped_DEGs_hpc <- merge(relaxed_DEGs_hpc, orthologs, by = "Gene_stable_ID")

# Keep only DEGs with a human gene symbol (i.e., successful mapping)
mapped_DEGs_hpc <- mapped_DEGs_hpc[!is.na(mapped_DEGs_hpc$Human_gene_name) & mapped_DEGs_hpc$Human_gene_name != "", ]

# Example: cortex mapped DEGs
gost_hpc <- gost(mapped_DEGs_hpc$Human_gene_name, organism = "hsapiens", significant = TRUE)
nrow(gost_hpc$result)

# See top terms
head(gost_hpc$result[, c("term_name", "intersection_size", "p_value")])


library(ggplot2)

# Extract relevant info
go_df <- gost_hpc$result
go_df <- go_df[order(go_df$p_value), ][1:10, ]  # top 10 terms

# Plot
ggplot(go_df, aes(x = reorder(term_name, -log10(p_value)), 
                  y = -log10(p_value), 
                  size = intersection_size)) +
  geom_point(aes(color = -log10(p_value))) +
  coord_flip() +
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  labs(
    title = " ",
    x = "GO Biological Process",
    y = "-log10(p-value)",
    size = "Gene Count"
  ) +
  theme_minimal()



#orthology mapping of striatum and go enrichment
# Load orthologs correctly with proper column names
orthologs <- read.csv("mart_export striatum.txt", header = TRUE, stringsAsFactors = FALSE)
colnames(orthologs) <- c("Gene_stable_ID", "Human_gene_stable_ID", "Human_gene_name")

# Add Gene ID as a column to your DEGs
relaxed_DEGs$Gene_stable_ID <- rownames(relaxed_DEGs)
relaxed_DEGs_str_df <- as.data.frame(relaxed_DEGs)
relaxed_DEGs_str_df$Gene_stable_ID <- rownames(relaxed_DEGs_str_df)

mapped_DEGs_str <- merge(relaxed_DEGs_str_df, orthologs, by = "Gene_stable_ID")
# Merge DEGs with orthologs using mouse Gene ID

head(mapped_DEGs_str)

# Keep only DEGs with a human gene symbol (i.e., successful mapping)
mapped_DEGs_str <- mapped_DEGs_str[!is.na(mapped_DEGs_str$Human_gene_name) & mapped_DEGs_str$Human_gene_name != "", ]
gene_list_str <- unique(na.omit(mapped_DEGs_str$Human_gene_name))
length(gene_list_str)
head(gene_list_str)
gost_str <- gost(gene_list_str, organism = "hsapiens", significant = TRUE)
gost_str <- gost(gene_list_str, organism = "hsapiens", significant = FALSE)
nrow(gost_str$result)  # Should be > 0

# Example: cortex mapped DEGs
gost_str <- gost(mapped_DEGs_str$Human_gene_name, organism = "hsapiens", significant = TRUE)
nrow(gost_str)
# See top terms
head(gost_str$result[, c("term_name", "intersection_size", "p_value")])


# Extract relevant info
go_df <- gost_str$result
go_df <- go_df[order(go_df$p_value), ][1:10, ]  # top 10 terms


# Plot
ggplot(go_df, aes(x = reorder(term_name, -log10(p_value)), 
                  y = -log10(p_value), 
                  size = intersection_size)) +
  geom_point(aes(color = -log10(p_value))) +
  coord_flip() +
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  labs(
    title = "Striatum ",
    x = "GO Biological Process",
    y = "-log10(p-value)",
    size = "Gene Count"
  ) +
  theme_minimal()


#pathvisio pathway genes go enrichment
# Read PathVisio gene list (1 column with gene symbols)
pathvisio_genes <- read.csv("shank3 gene list.txt", header = TRUE, stringsAsFactors = FALSE)

# Make sure it's a clean character vector
gene_list <- unique(pathvisio_genes[[1]])  # or: pathvisio_genes$Symbol, depending on column name
# Remove everything after the tab character
gene_list <- sub("\t.*", "", gene_list)

# Check now
head(gene_list)

length(gene_list)



gost_result <- gost(
  query = gene_list,
  organism = "hsapiens",
  sources = "GO:BP",         # Or add KEGG, Reactome etc.
  significant = TRUE
)


# Take top 10 terms
go_df <- gost_result$result
go_df <- go_df[order(go_df$p_value), ][1:10, ]

n_enriched_terms <- nrow(go_df)
print(n_enriched_terms)
go_df[ , c("term_name", "intersection_size")]

# Plot pathvisio genes
ggplot(go_df, aes(x = reorder(term_name, -log10(p_value)), 
                  y = -log10(p_value), 
                  size = intersection_size)) +
  geom_point(aes(color = -log10(p_value))) +
  coord_flip() +
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  labs(
    title = "PathVisio Genes",
    x = "GO Biological Process",
    y = "-log10(p-value)",
    size = "Gene Count"
  ) +
  theme_minimal()







#volcano plots for all the three tissue

plot_volcano_top10 <- function(res_df, title) {
  res_df <- as.data.frame(res_df)
  
  # Assign direction
  res_df$direction <- "Not Sig"
  res_df$direction[res_df$pvalue < 0.05 & res_df$log2FoldChange > 0.5] <- "Upregulated"
  res_df$direction[res_df$pvalue < 0.05 & res_df$log2FoldChange < -0.5] <- "Downregulated"
  
  # Remove missing or empty HGNC symbols
  res_df <- res_df[res_df$hgnc_symbol != "" & !is.na(res_df$hgnc_symbol), ]
  
  # Select top 5 by lowest p-value
  top5 <- res_df[order(res_df$pvalue), ][1:5, ]
  
  # Add labels
  res_df$label <- ifelse(res_df$hgnc_symbol %in% top5$hgnc_symbol, res_df$hgnc_symbol, "")
  
  # Plot
  ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
    geom_point(aes(color = direction), alpha = 0.6) +
    geom_text_repel(aes(label = label), size = 3.5, max.overlaps = Inf) +
    scale_color_manual(values = c(
      "Upregulated" = "red",
      "Downregulated" = "blue",
      "Not Sig" = "gray"
    )) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    labs(
      title = title,
      x = "Log2 Fold Change",
      y = "-log10(p-value)"
    ) +
    xlim(-23, 10) +
    ylim(0, 20) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
}

res_hippocampus_annotated <- read.csv("full_DESeq2_mapped_hippocampus.csv", stringsAsFactors = FALSE)
res_striatum_annotated <- read.csv("full_DESeq2_mapped_striatum.csv", stringsAsFactors = FALSE)
res_cortex_annotated <- read.csv("full_DESeq2_mapped_cortex.csv", stringsAsFactors = FALSE)

# Remove empty HGNC symbols
res_cortex_cleaned <- res_cortex_annotated %>%
  filter(hgnc_symbol != "" & !is.na(hgnc_symbol))

# Keep only one row per gene symbol — lowest p-value
res_cortex_cleaned <- res_cortex_cleaned %>%
  group_by(hgnc_symbol) %>%
  slice_min(order_by = pvalue, n = 1) %>%
  ungroup()

#plot the volcano plots with top 5 genes
plot_volcano_top10(res_cortex_cleaned, "Volcano Plot: HT vs WT (Cortex)")
plot_volcano_top10(res_hippocampus_annotated, "")
plot_volcano_top10(res_striatum_annotated, "")
