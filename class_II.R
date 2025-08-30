classify_gene <- function(logFC, padj) {
  ifelse(padj < 0.05 & logFC > 1, "Upregulated",
         ifelse(padj < 0.05 & logFC < -1, "Downregulated", "Not_Significant"))
}
input_dir <- "raw_data" 
output_dir <- "results"
result_list <- list()

if(!dir.exists(output_dir)){
  dir.create(output_dir)
}

files_to_process <- c("DEGs_Data_1.csv", "DEGs_Data_2.csv") 

result_list <- list()

for (file_names in files_to_process) {
  cat("\nProcessing:", file_names, "\n")
  
  input_file_path <- file.path(input_dir, file_names)
  
  data <- read.csv(input_file_path, header = TRUE)
  cat("File imported. Checking for missing values...\n")
  summary(data)
  
  if("padj" %in% names(data)){
    missing_count <- sum(is.na(data$padj))
    
    cat("Missing values in 'padj':", missing_count, "\n")
    data$padj[is.na(data$padj)] <- 1
  }
  data$status <- classify_gene(data$logFC, data$padj)
  result_list[[file_names]] <- data 
  
  output_file_path <- file.path(output_dir, paste0("Processed_", file_names))
  write.csv(data, output_file_path, row.names = FALSE)
  cat("Results saved to:", output_file_path, "\n")
}
results_1 <- result_list[[1]] 
results_2 <- result_list[[2]]

cat("\nSummary for DEGs_Data_1.csv:\n")
print(table(results_1$status))

cat("\nSummary for DEGs_Data_2.csv:\n")
print(table(results_2$status))
