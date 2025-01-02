# Differential Gene Expression Analysis Script for Alzheimer's Dataset (GSE125583)

# Loading the raw data
# This file contains the raw RNA-Seq counts from the GSE125583 dataset.
data <- read.delim("GSE125583_raw_counts_GRCh38.p13_NCBI.tsv", header = TRUE, sep = "\t")
View(data)  # Checking the structure of the dataset

# Importing necessary libraries
library("dplyr")

# Selecting relevant columns
# Keeping only the GeneID and the counts for the samples being analyzed
cts <- select(data, GeneID, GSM3577568, GSM3577569, GSM3577570, GSM3577571, GSM3577572, 
              GSM3577580, GSM3577583, GSM3577584, GSM3577587, GSM3577589)
View(cts)

# Saving the filtered count table for reference
write.csv(cts, file = "proj_data_count_table_1.csv", row.names = FALSE)

# Preparing the count data for DESeq2 analysis
cts1 <- cts[, -1]  # Removing the GeneID column to isolate count data
View(cts1)

# Creating metadata for DESeq2
# This metadata specifies whether each column corresponds to a sample or control
colData <- data.frame(
  condition = c("sample", "sample", "sample", "sample", "sample", 
                "control", "control", "control", "control", "control"),
  row.names = colnames(cts1)
)
write.csv(colData, file = "colData_1.csv", row.names = FALSE)
View(colData)

# Setting row names of count data to match Gene IDs
rownames(cts1) <- cts$GeneID

# Verifying that the metadata matches the count data
all(rownames(colData) == colnames(cts1))

# Loading DESeq2 for differential expression analysis
library("DESeq2")

# Creating a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = cts1,
                              colData = colData,
                              design = ~ condition)
dds

# Specifying the reference level for the condition factor
dds$condition <- factor(dds$condition, levels = c("sample", "control"))

# Running the DESeq pipeline for differential expression analysis
dds <- DESeq(dds)
result <- results(dds)
result  # Displaying the results

# Generating results for a specific contrast
res <- results(dds, contrast = c("condition", "sample", "control"))
res

# Creating an MA plot to visualize the results
plotMA(res, ylim = c(-5, 5))

# Plotting counts for the most significant gene
plotCounts(dds, gene = which.min(res$padj), intgroup = "condition")

# Extracting significantly upregulated and downregulated genes
resSigind <- res[which(res$padj < 0.05 & res$log2FoldChange > 2), ]
resSigrep <- res[which(res$padj < 0.05 & res$log2FoldChange < -2), ]
resSig <- rbind(resSigind, resSigrep)
resSig

# Saving the upregulated and downregulated genes separately
upregulated <- data.frame(rownames(resSigind))
downregulated <- data.frame(rownames(resSigrep))
View(upregulated)

# Annotating genes using a reference annotation file
anootated_file <- read.delim("Human.GRCh38.p13.annot.tsv", header = TRUE, sep = "\t")
View(anootated_file)

# Matching upregulated genes to the annotation file
gene_match_rows = match(rownames(resSigind), anootated_file[, 1])
upreg_annot <- anootated_file[gene_match_rows, c(1, 2, 3)]
View(upreg_annot)

# Matching downregulated genes to the annotation file
gene_match_rows_downreg = match(rownames(resSigrep), anootated_file[, 1])
down_reg_annot <- anootated_file[gene_match_rows_downreg, c(1, 2, 3)]
View(down_reg_annot)

# Saving the annotated results
write.csv(upreg_annot, file = "upregulated_genes_1.csv", row.names = FALSE)
write.csv(down_reg_annot, file = "downregulated_genes_1.csv", row.names = FALSE)

# Performing gene enrichment analysis
library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)

# Preparing gene list for enrichment analysis
original_gene_list <- res$log2FoldChange
names(original_gene_list) <- rownames(res)
original_gene_list <- na.omit(original_gene_list)
