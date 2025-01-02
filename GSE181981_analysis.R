# Loading and preparing alcohol-related dataset

# Viewing the raw alcohol dataset to check its structure
View(alcohol_raw_data)

# Removing the first column (index ) to isolate count data
alcohol_raw_data_1 <- alcohol_raw_data[, -1]
View(alcohol_raw_data_1)

# Setting Gene IDs as row names for the dataset
rownames(alcohol_raw_data_1) <- alcohol_raw_data$GeneID

# Creating metadata for samples
# Metadata includes sample IDs and descriptions indicating whether a sample is a control or from subjects with Alcohol Use Disorder
metadata <- data.frame(
  sample_id = c(
    "GSM5515131", "GSM5515132", "GSM5515133", "GSM5515134", "GSM5515135",
    "GSM5515136", "GSM5515137", "GSM5515138", "GSM5515139", "GSM5515140",
    "GSM5515141", "GSM5515142", "GSM5515143", "GSM5515144", "GSM5515145",
    "GSM5515146", "GSM5515147", "GSM5515148", "GSM5515149", "GSM5515150",
    "GSM5515151", "GSM5515152", "GSM5515153", "GSM5515154"
  ),
  description = c(
    "Subject_1_Control", "Subject_2_Control", "Subject_3_Control", "Subject_4_AlcoholUseDisorder", "Subject_5_AlcoholUseDisorder",
    "Subject_6_AlcoholUseDisorder", "Subject_7_AlcoholUseDisorder", "Subject_8_Control", "Subject_9_AlcoholUseDisorder", "Subject_10_Control",
    "Subject_11_Control", "Subject_12_AlcoholUseDisorder", "Subject_13_Control", "Subject_14_AlcoholUseDisorder", "Subject_15_Control",
    "Subject_16_AlcoholUseDisorder", "Subject_17_Control", "Subject_18_AlcoholUseDisorder", "Subject_19_Control", "Subject_20_Control",
    "Subject_21_AlcoholUseDisorder", "Subject_22_Control", "Subject_23_AlcoholUseDisorder", "Subject_24_AlcoholUseDisorder"
  )
)

# Assigning conditions (sample or control) based on the description in metadata
metadata$condition <- ifelse(grepl("Control", metadata$description), "control", "sample")

# Creating colData for DESeq2
# colData specifies the condition (sample/control) for each sample
colData_1 <- data.frame(
  condition = metadata$condition,
  row.names = metadata$sample_id
)

# Viewing the resulting colData structure
View(colData_1)

# Ensuring that the row names of colData match the column names of the count data
all(rownames(colData_1) == colnames(alcohol_raw_data_1))

# Performing differential gene expression analysis using DESeq2
library(DESeq2)

# Creating a DESeqDataSet object
dds_alc <- DESeqDataSetFromMatrix(countData = alcohol_raw_data_1,
                                  colData = colData_1,
                                  design = ~ condition)
dds_alc

# Setting the reference level for conditions
dds_alc$condition <- factor(dds_alc$condition, levels = c("sample", "control"))

# Running the DESeq analysis pipeline
dds_alc <- DESeq(dds_alc)
result_alc <- results(dds_alc)
result_alc  # Viewing the overall results

# Generating specific contrasts for results
res_alc <- results(dds_alc, contrast = c("condition", "sample", "control"))
res_alc

# Creating an MA plot to visualize differential expression
plotMA(res_alc, ylim = c(-5, 5))

# Plotting counts for the most significant gene
plotCounts(dds_alc, gene = which.min(res_alc$padj), intgroup = "condition")

# Extracting upregulated and downregulated genes based on thresholds
res_alc_up <- res_alc[which(res_alc$padj < 0.05 & res_alc$log2FoldChange > 2), ]
res_alc_down <- res_alc[which(res_alc$padj < 0.05 & res_alc$log2FoldChange < -2), ]
res_up_down_total <- rbind(res_alc_up, res_alc_down)
res_up_down_total

# Saving upregulated and downregulated genes to data frames
upregulated_genes_alz <- data.frame(rownames(res_alc_up))
downregulated_genes_alz <- data.frame(rownames(res_alc_down))
View(upregulated_genes_alz)
View(downregulated_genes_alz)

# Annotating genes using a reference annotation file
anootated_file <- read.delim("Human.GRCh38.p13.annot.tsv", header = TRUE, sep = "\t")
View(anootated_file)

# Annotating upregulated genes
gene_match_rows_alz_upregulated <- match(rownames(res_alc_down), anootated_file[, 1])
upreg_alz_annot <- anootated_file[gene_match_rows_alz_upregulated, c(1, 2, 3)]
View(data.frame(gene_match_rows_alz_upregulated))
View(upreg_alz_annot)

# Annotating downregulated genes
gene_match_rows_downreg <- match(rownames(res_alc_up), anootated_file[, 1])
downreg_alz_annot <- anootated_file[gene_match_rows_downreg, c(1, 2, 3)]
View(downreg_alz_annot)

# Saving annotated results to CSV files
write.csv(upreg_alz_annot, file = "upregulated_genes_alz.csv", row.names = FALSE)
write.csv(downreg_alz_annot, file = "downregulated_genes_alz.csv", row.names = FALSE)

# Extracting gene description for a specific gene
library(dplyr)
gene_description <- anootated_file %>%
  filter(GeneID == 999) %>%
  pull(Description, Symbol)
gene_description
