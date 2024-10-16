#!/usr/bin/env Rscript

# DESCRIPTION:
# This R script performs the following tasks:
# 1. Normalizes count data using the TMM method and calculates log2 counts per million (CPM) plus 1.
# 2. Performs differential expression analysis using the limma package, comparing bacterial, viral, and control groups based on an external annotation file.
# 3. Writes the results of each differential expression comparison (B vs C, V vs C, and B vs V) into one .csv file for each comparison.

# USAGE:
# Rscript 08_1_runLIMMA_DE.R <countData.txt> <annotationFile.txt> <outputDirectory>

# Load required packages
library(limma)
library(edgeR)

# Function to normalize counts using TMM method
normalize_counts <- function(input_file, output_file) {
  countTable <- read.delim(input_file, stringsAsFactors = FALSE, check.names = FALSE, header = TRUE)
  counts <- as.matrix(countTable[, -1])  # Exclude the first column (gene IDs)
  dgeLst <- DGEList(counts = counts)
  dgeLst$genes <- countTable$miRNA
  dgeTMM <- calcNormFactors(dgeLst, method = "TMM")
  pseudo_TC <- log2(cpm(dgeTMM) + 1)
  rownames(pseudo_TC) <- countTable$miRNA
  write.table(pseudo_TC, file = output_file, sep = "\t", quote = FALSE, row.names = TRUE)
}

# Function to perform differential expression analysis
perform_differential_analysis <- function(count_table_file, annotation_file, output_file, contrast_name) {
  countTable <- read.delim(count_table_file, stringsAsFactors = FALSE, check.names = FALSE, header = TRUE)
  annotation <- read.delim(annotation_file, stringsAsFactors = FALSE, header = TRUE)
  
  # Filter the annotation file to keep only samples present in the count table
  sampleNames <- colnames(countTable)[-1]  # Exclude gene ID column
  annotation <- annotation[annotation$sample %in% sampleNames, ]
  
  # Reorder the count table based on the order of samples in the annotation file
  countTable <- countTable[, c("miRNA", annotation$sample)]
  
  # Define groups using the 'condition' column from the annotation file
  groups <- as.factor(annotation$condition)
  
  # Create a design matrix for the linear model
  design <- model.matrix(~0 + groups)
  colnames(design) <- levels(groups)
  
  # Fit a linear model using the already normalized data (logCPM)
  counts <- as.matrix(countTable[, -1])  # Exclude gene IDs
  fit <- lmFit(counts, design)
  fit$genes <- countTable$miRNA
  
  # Define the contrast matrix and fit the linear model with contrasts
  cont.matrix <- makeContrasts(contrast_name, levels = design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  
  # Calculate moderated t-statistics using empirical Bayes
  eBayesFit <- eBayes(fit2, trend = TRUE)
  
  # Extract the top differentially expressed genes and write them to a single file
  allGenes <- topTable(eBayesFit, number = Inf)
  write.table(allGenes, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
}

# Main function to run all steps
main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (length(argv) != 3) {
    stop("USAGE: Rscript 08_runLIMMA_DE.R <countData.txt> <annotationFile.txt> <outputDirectory>", call. = FALSE)
  }
  
  countDataFile <- argv[1]
  annotationFile <- argv[2]
  outputDir <- argv[3]
  
  # Step 1: Normalize counts
  normalizedCountsFile <- file.path(outputDir, 'normalizedTMMlog2Cpmplus1_LIMMA.csv')
  normalize_counts(countDataFile, normalizedCountsFile)
  
  # Step 2: Differential analysis for Bacterial vs Control (B vs C)
  perform_differential_analysis(
    countDataFile,
    annotationFile,
    file.path(outputDir, '_DEGsUsingLIMMA_B_vs_C.csv'),
    'groupsB - groupsC'  # Contrast for Bacterial vs Control
  )
  
  # Step 3: Differential analysis for Viral vs Control (V vs C)
  perform_differential_analysis(
    countDataFile,
    annotationFile,
    file.path(outputDir, '_DEGsUsingLIMMA_V_vs_C.csv'),
    'groupsV - groupsC'  # Contrast for Viral vs Control
  )
  
  # Step 4: Differential analysis for Bacterial vs Viral (B vs V)
  perform_differential_analysis(
    countDataFile,
    annotationFile,
    file.path(outputDir, '_DEGsUsingLIMMA_B_vs_V.csv'),
    'groupsB - groupsV'  # Contrast for Bacterial vs Viral
  )
  
  cat("LIMMA analysis complete. Results saved as CSV files in the specified output directory.\n")
}

# Run the main function
main()
