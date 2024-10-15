#!/usr/bin/env Rscript

# DESCRIPTION:
# This R script performs the following tasks:
# 1. Normalizes count data using the TMM method and calculates log2 counts per million (CPM) plus 1.
# 2. Performs differential expression analysis using the limma package, comparing bacterial, viral, and control groups.
# 3. Extracts and writes the top differentially expressed genes to files.
# 4. Sorts and writes results to separate TSV and XLSX files based on log-fold change and p-value.

# USAGE:
# Rscript 05_runLIMMA_DE.R

# Load required packages
library(limma)
library(edgeR)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(openxlsx)

# Main function to normalize counts
normalize_counts <- function(input_file, output_file) {
  countTable <- read.delim(input_file, stringsAsFactors = FALSE, check.names = FALSE, header = TRUE)
  counts <- as.matrix(countTable[, -1])
  dgeLst <- DGEList(counts = counts)
  dgeLst$genes <- countTable$miRNA
  dgeTMM <- calcNormFactors(dgeLst, method = "TMM")
  pseudo_TC <- log2(cpm(dgeTMM) + 1)
  rownames(pseudo_TC) <- countTable$miRNA
  write.table(pseudo_TC, file = output_file, sep = "\t", quote = FALSE, row.names = TRUE)
}

# Main function to perform differential expression analysis for specific group comparisons
perform_differential_analysis <- function(count_table_file, output_file, count_columns, group_definitions, contrast_name) {
  countTable <- read.delim(count_table_file, stringsAsFactors = FALSE, check.names = FALSE, header = TRUE)
  
  # Combine the counts for the samples as specified in count_columns
  counts <- cbind(countTable[, count_columns])
  
  # Define group labels
  groups <- as.factor(group_definitions)
  
  # Create a design matrix for the linear model
  design <- model.matrix(~0 + groups)
  colnames(design) <- levels(groups)
  
  # Fit a linear model using the already normalized data (logCPM)
  fit <- lmFit(counts, design)
  fit$genes <- countTable$miRNA
  
  # Define the contrast matrix and fit the linear model with contrasts
  cont.matrix <- makeContrasts(contrast_name, levels = design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  
  # Calculate moderated t-statistics using an empirical Bayes approach
  eBayesFit <- eBayes(fit2, trend = TRUE)
  
  # Extract the top differentially expressed genes and write them to a file
  allGenes <- topTable(eBayesFit, number = Inf)
  write.table(allGenes, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
}

# Main function to sort results by logFC and p-value
sort_by_fc_pval <- function(input_file, output_up_file, output_down_file) {
  inDF <- read.delim(input_file, stringsAsFactors = FALSE, check.names = FALSE)
  up <- subset(inDF, logFC > 0)
  down <- subset(inDF, logFC < 0)
  ordDFup <- up[order(-up$logFC, up$adj.P.Val), ]
  ordDFdown <- down[order(down$logFC, down$adj.P.Val), ]
  
  # Write sorted data frames to TSV format
  write.table(ordDFup, file = paste0(output_up_file, ".tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  write.table(ordDFdown, file = paste0(output_down_file, ".tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  # Write sorted data frames to XLSX format
  wb <- createWorkbook()
  addWorksheet(wb, sheetName = "Upregulated")
  addWorksheet(wb, sheetName = "Downregulated")
  writeData(wb, sheet = "Upregulated", x = ordDFup)
  writeData(wb, sheet = "Downregulated", x = ordDFdown)
  saveWorkbook(wb, paste0(output_up_file, ".xlsx"), overwrite = TRUE)
  saveWorkbook(wb, paste0(output_down_file, ".xlsx"), overwrite = TRUE)
}

# Main function to run all steps
main <- function() {
  # Step 1: Normalize counts
  normalize_counts('../01_featureCounts/countDataFilt10Counts.txt', '../01_featureCounts/normalizedTMMlog2Cpmplus1.csv')
  
  # Step 2: Differential analysis with defined sample groups
  
  # Comparison 1: Bacterial (B) vs Control (C)
  perform_differential_analysis(
    'countDataFilt10Counts_DE.txt',
    '03_1_DEGs_B_vs_C.csv',
    2:8, 16:23,  # First 7 columns for B, last 8 for C
    c(rep("B", 7), rep("C", 8)),
    'B-C'
  )
  
  # Comparison 2: Viral (V) vs Control (C)
  perform_differential_analysis(
    'countDataFilt10Counts_DE.txt',
    '03_2_DEGs_V_vs_C.csv',
    9:15, 16:23,  # Columns 9-15 for V, last 8 for C
    c(rep("V", 7), rep("C", 8)),
    'V-C'
  )
  
  # Comparison 3: Bacterial (B) vs Viral (V)
  perform_differential_analysis(
    'countDataFilt10Counts_DE.txt',
    '03_3_DEGs_B_vs_V.csv',
    2:8, 9:15,  # First 7 columns for B, columns 9-15 for V
    c(rep("B", 7), rep("V", 7)),
    'B-V'
  )
  
  # Step 3: Sort results by logFC and p-value
  sort_by_fc_pval('03_1_DEGs_B_vs_C.csv', '03_1_DEGs_B_vs_C_UP', '03_1_DEGs_B_vs_C_DOWN')
  sort_by_fc_pval('03_2_DEGs_V_vs_C.csv', '03_2_DEGs_V_vs_C_UP', '03_2_DEGs_V_vs_C_DOWN')
  sort_by_fc_pval('03_3_DEGs_B_vs_V.csv', '03_3_DEGs_B_vs_V_UP', '03_3_DEGs_B_vs_V_DOWN')
}

# Run main function
main()
