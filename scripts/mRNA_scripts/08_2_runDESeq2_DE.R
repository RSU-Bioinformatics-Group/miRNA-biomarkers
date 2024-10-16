#!/usr/bin/env Rscript

# DESCRIPTION:
# This R script performs differential expression analysis using DESeq2.
# It reads in an external annotation file to define sample groups and 
# runs comparisons between Bacterial vs Control, Viral vs Control, and Bacterial vs Viral.

# USAGE:
# Rscript 08_2_run_DESeq2_DE.R <countData.txt> <annotationFile.txt> <outputFilePrefix>

# Load required libraries
library("DESeq2")
library("IHW")

# Function to load the count data and metadata (sample group annotations)
load_data <- function(count_table_file, annotation_file) {
  # Step 1: Read in the count data
  countTable <- read.delim(count_table_file, stringsAsFactors = FALSE, check.names = FALSE, header = TRUE)
  
  # Step 2: Read in the annotation file
  annotations <- read.delim(annotation_file, stringsAsFactors = FALSE, check.names = FALSE, header = TRUE)
  
  # Step 3: Ensure gene IDs are unique (if necessary)
  rownames(countTable) <- make.unique(countTable$Gene_ID)
  
  return(list(countTable = countTable, annotations = annotations))
}

# Function to run DESeq2 analysis and save the results
run_deseq2_analysis <- function(counts, meta, condition_levels, output_file) {
  # Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~ condition)
  
  # Set factor levels
  dds$condition <- factor(dds$condition, levels = condition_levels)
  
  # Run DESeq2 analysis
  dds <- DESeq(dds)
  
  # Extract results and order by adjusted p-value (padj)
  res <- results(dds)
  resOrdered <- res[order(res$padj),]
  
  # Write the ordered results to a file
  write.table(resOrdered, file = output_file, sep = "\t", quote = FALSE, row.names = TRUE)
}

# Main function to handle the DE comparisons
main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (length(argv) != 3) {
    stop("USAGE: Rscript 08_2_run_DESeq2_DE.R <countData.txt> <annotationFile.txt> <outputFilePrefix>", call. = FALSE)
  }
  
  countDataFile <- argv[1]
  annotationFile <- argv[2]
  outputFilePrefix <- argv[3]
  
  # Load the count data and annotation file
  data <- load_data(countDataFile, annotationFile)
  countTable <- data$countTable
  annotations <- data$annotations
  
  # Ensure the counts matrix contains only the relevant sample columns
  counts <- countTable[, annotations$SampleID]
  rownames(counts) <- rownames(countTable)  # Set rownames to the gene IDs
  
  # Comparison 1: Bacterial (B) vs Control (C)
  meta_BC <- data.frame(condition = factor(annotations$Group[annotations$Group %in% c('Bacterial', 'Control')],
                                           levels = c('Bacterial', 'Control')))
  counts_BC <- counts[, annotations$Group %in% c('Bacterial', 'Control')]
  
  output_file_BC <- paste0(outputFilePrefix, "_DEGsUsingDESeq2_B_vs_C.csv")
  run_deseq2_analysis(counts_BC, meta_BC, c('Control', 'Bacterial'), output_file_BC)
  
  # Comparison 2: Viral (V) vs Control (C)
  meta_VC <- data.frame(condition = factor(annotations$Group[annotations$Group %in% c('Viral', 'Control')],
                                           levels = c('Viral', 'Control')))
  counts_VC <- counts[, annotations$Group %in% c('Viral', 'Control')]
  
  output_file_VC <- paste0(outputFilePrefix, "_DEGsUsingDESeq2_V_vs_C.csv")
  run_deseq2_analysis(counts_VC, meta_VC, c('Control', 'Viral'), output_file_VC)
  
  # Comparison 3: Bacterial (B) vs Viral (V)
  meta_BV <- data.frame(condition = factor(annotations$Group[annotations$Group %in% c('Bacterial', 'Viral')],
                                           levels = c('Bacterial', 'Viral')))
  counts_BV <- counts[, annotations$Group %in% c('Bacterial', 'Viral')]
  
  output_file_BV <- paste0(outputFilePrefix, "_DEGsUsingDESeq2_B_vs_V.csv")
  run_deseq2_analysis(counts_BV, meta_BV, c('Viral', 'Bacterial'), output_file_BV)
  
  cat("DESeq2 analysis complete. Results saved to files.\n")
}

# Run the main function
main()