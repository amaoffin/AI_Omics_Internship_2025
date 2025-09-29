# AI OMICS INTERNSHIP 2025

Research internship conducted by AI and Biotechnology to gain hands-on experience in R programming for biological data science, including scripting, data wrangling, and statistical analysis of high-dimensional omics datasets.

## class_Ib.R

Basics of R programming and handling datasets; cleaning and factoring patient datasets, as well as summarizing data.

## class_Ic.R

Practice using basic R functions and operations.

## class_II.R

Working with the results of differential gene expression (DGE) analysis, particularly:

  log2FoldChange (logFC): magnitude and direction of change in gene expression; positive values suggest upregulation, negative values suggest downregulation
  
  adjusted p-value (padj): statistical significance of the observed difference, corrected for multiple testing; padj < 0.05 was used
  
In this project, I wrote a function classify_gene() to categorize genes as upregulated, downregulated, or insignificant based on thresholds of log2FC and padj, processed two provided datasets using this function, and generated summary counts of upregulated and downregulated genes.
  
## class_III.R
Full preprocessing workflow for a microarray expression dataset, including quality control, normalization, filtering, and phenotype labeling.

### Quality Control
Performed quality control checks before and after normalization, using array-level diagnostic plots to detect outliers based on metrics such as signal distribution, MA plots, and boxplots.

### Normalization & Filtering

Applied background correction and quantile normalization to the expression dataset. After normalization, filtered out low-intensity probes using an intensity threshold and detection p-values.

### Phenotype Annotation

Used associated phenotype metadata to define the experimental groups. Original labels were mapped to simplified categories:

"normal" for control/non-diseased samples

"diseased" for non-control/diseased samples

These relabeled groups were used for downstream differential expression and classification analysis.
