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
- **GSE181981**: RNA-Seq data from the amygdala of 12 AUD and 12 control samples.https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi

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

## Requirements
### Software and Tools:
- R (version >= 4.0)
- RStudio
- GEOquery, DESeq2, ggplot2

### R Packages:
```r
install.packages(c("GEOquery", "DESeq2", "ggplot2"))
```


## References
1. Joshi, A., Giorgi, F. M., & Sanna, P. P. (2024). Transcriptional Patterns in Stages of Alzheimer's Disease Are Cell-Type Specific and Partially Converge with the Effects of Alcohol Use Disorder in Humans. *ENeuro*, 11(10). https://doi.org/10.1523/ENEURO.0118-24.2024
2. Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15(12), 550. https://doi.org/10.1186/s13059-014-0550-8


