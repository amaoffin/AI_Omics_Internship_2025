AI OMICS INTERNSHIP 2025
Research internship conducted by AI and Biotechnology to gain hands-on experience in R programming for biological data science, including scripting, data wrangling, and statistical analysis of high-dimensional omics datasets.

class_Ib.R
Basics of R programming and handling datasets; cleaning and factoring patient datasets, as well as summarizing data.

class_Ic.R
Practice using basic R functions and operations.

class_II.R
Working with the results of differential gene expression (DGE) analysis, particularly:
  log2FoldChange (logFC): magnitude and direction of change in gene expression; positive values suggest upregulation, negative values suggest downregulation
  adjusted p-value (padj): statistical significance of the observed difference, corrected for multiple testing; padj < 0.05 was used
In this project, I wrote a function classify_gene() to categorize genes as upregulated, downregulated, or insignificant based on thresholds of log2FC and padj, processed two provided datasets using this function, and generated summary counts of upregulated and downregulated genes.
  
