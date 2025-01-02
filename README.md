# Differential Gene Expression Analysis: Alzheimer's Disease and Alcohol Use Disorder

## Project Overview
This repository contains R scripts and documentation for the analysis of differential gene expression (DEG) in Alzheimer's Disease (AD) and Alcohol Use Disorder (AUD) datasets. The project aims to identify common genes implicated in these conditions using publicly available RNA-Seq datasets from the Gene Expression Omnibus (GEO).

## Objectives
1. Analyze RNA-Seq data from:
   - GSE125583: Alzheimer's disease dataset.
   - GSE181981: Alcohol Use Disorder dataset.
2. Perform differential expression analysis to identify DEGs for each condition.
3. Investigate potential overlapping DEGs between the two datasets.

## Methodology
### 1. Data Collection
The datasets were downloaded from GEO:
- **[GSE125583]**: RNA-Seq data from the fusiform gyrus of 5 AD and 5 control post-mortem. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125583
- **GSE181981**: RNA-Seq data from the amygdala of 12 AUD and 12 control sampleshttps://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181981

### 2. Preprocessing and Differential Expression Analysis
#### Alzheimer's Dataset (GSE125583):
1. Preprocessed raw count matrix by:
   - Filtering samples (5 AD, 5 controls).
   - Setting `GeneID` as row names.
   - Creating metadata (`colData`) for experimental conditions.
2. Performed DEG analysis using the `DESeq2` package in R with the following thresholds:
   - |log2FoldChange| > 2
   - Adjusted p-value < 0.05
3. Annotated DEGs using a reference annotation file.

#### Alcohol Dataset (GSE181981):
1. Followed the same preprocessing and analysis steps as the AD dataset.

#### Identification of Common Genes:
- Identified overlapping DEGs between the two datasets using R.

### 3. Results
#### Alzheimer's Dataset:
- **Upregulated Genes:** 83
- **Downregulated Genes:** 242

#### Alcohol Dataset:
- **Upregulated Genes:** 2
- **Downregulated Genes:** 3

#### Common Genes:
- No common genes were identified between the two datasets.

## Key Findings
1. **Alzheimer's Dataset:**
   - The most highly expressed gene was AFAP1L2 (GeneID: 84632).
2. **Alcohol Dataset:**
   - The most highly expressed gene was CDH1 (GeneID: 999).
3. **No Overlap:**
   - The lack of overlapping DEGs may be due to differences in brain regions, experimental designs, or biological variation.

## plots 
**MA plot showing differential gene expression between Alzheimer’s and control samples.**

![image](https://github.com/user-attachments/assets/574034a4-25c8-4267-aaf0-59927d74a60f)
 -A total of 83 genes are significantly upregulated (blue dots above the horizontal line), while 242 genes are downregulated (blue dots below the line). The x-axis represents the mean normalized counts (log scale), and the y-axis shows the log2 fold change in expression. This plot highlights key genes potentially contributing to Alzheimer’s pathology or altered due to disease-related processes.

 **Most highly expressed gene in the  Alzheimer’s dataset compared to the normal.**

 ![image](https://github.com/user-attachments/assets/e5d6c9c2-5695-4cfe-ae9b-5fc8f077d76b)
-Plot of normalized counts for the most significantly differentially expressed gene (ID: 84632) identified by DESeq2 (lowest adjusted p-value). **GeneID 84632 (AFAP1L2)** is actin filament associated protein 1 like 2.

**MA plot showing differential gene expression between Alcohol affected  and control samples.**

![image](https://github.com/user-attachments/assets/9186df30-3312-4a5e-a614-4cde11a3c686)
-This MA plot visualizes the differential gene expression analysis for the alcohol dataset. Each point represents a gene, plotted with its mean normalized expression counts (x-axis) against its log fold change (y-axis). Highlighted in blue are 2 upregulated genes and 3 downregulated genes that passed the significance threshold. The majority of genes (in gray) showed no significant differential expression between conditions.

**Most highly expressed gene in the  Alcohol dataset compared to the normal.**

![image](https://github.com/user-attachments/assets/ef31442e-e752-417c-83bf-ed34168ab466)
-This MA plot visualizes the differential gene expression analysis for the alcohol dataset. Each point represents a gene, plotted with its mean normalized expression counts (x-axis) against its log fold change (y-axis).**Gene ID 999 is a (CDH1) cadherin 1.**

## References
1. Joshi, A., Giorgi, F. M., & Sanna, P. P. (2024). Transcriptional Patterns in Stages of Alzheimer's Disease Are Cell-Type Specific and Partially Converge with the Effects of Alcohol Use Disorder in Humans. *ENeuro*, 11(10). https://doi.org/10.1523/ENEURO.0118-24.2024
2. Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15(12), 550. https://doi.org/10.1186/s13059-014-0550-8


