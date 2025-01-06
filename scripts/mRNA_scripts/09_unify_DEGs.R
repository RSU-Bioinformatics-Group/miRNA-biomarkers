#!/usr/bin/env Rscript

# DESCRIPTION:
# This script compares DEGs between LIMMA and DESeq2 for the same comparison group, performs Fisher's meta-analysis to combine p-values, and calculates adjusted p-values.

# USAGE:
# Rscript 09_unify_DEGs.R

# Load necessary libraries
library(dplyr)
library(metap)

# Function to compare DEGs between LIMMA and DESeq2
compare_DEGs <- function(limma_file, deseq2_file) {
  # Load the LIMMA and DESeq2 files
  dat1 <- read.csv(limma_file, stringsAsFactors = FALSE, sep = "\t")
  dat2 <- read.csv(deseq2_file, stringsAsFactors = FALSE, sep = "\t")

  # Ensure IDs (gene symbols) are available in both data sets
  dat1$ID <- as.character(dat1$ID)
  dat2$ID <- as.character(row.names(dat2))

  # Adjust column names to reflect the source (LIMMA or DESeq2)
  colnames(dat1)[-1] <- paste0("LIMMA_", colnames(dat1)[-1])
  colnames(dat2)[-1] <- paste0("DESeq2_", colnames(dat2)[-1])

  # Merge the two datasets by the gene ID
  merged_data <- merge(dat1, dat2, by = "ID")

  return(merged_data)
}

# Function to perform Fisher's meta-analysis
perform_meta_analysis <- function(merged_data, output_file) {
  # Calculate the mean log fold change between LIMMA and DESeq2
  meanFC <- rowMeans(cbind(merged_data$LIMMA_logFC, merged_data$DESeq2_log2FoldChange), na.rm = TRUE)

  # Perform Fisher's method to combine p-values
  FishersP <- apply(rbind(merged_data$LIMMA_P.Value, merged_data$DESeq2_pvalue), 2, metap::sumlog)
  FishersPP <- sapply(FishersP, function(x) x$p)

  # Adjust p-values using Benjamini-Hochberg correction
  BH <- p.adjust(FishersPP, method = "BH")

  # Create a data frame for the results
  outList <- data.frame(
    Gene_Symbol = merged_data$ID,
    meanLogFC = meanFC,
    Fishers_Meta_PValue = FishersPP,
    adjPValueBH = BH
  )

  # Write the output to a file
  write.table(outList, output_file, sep = "\t", quote = FALSE, row.names = FALSE)

  # Print completion message
  cat(paste0("Meta-analysis completed. Results saved to: ", output_file, "\n"))
}

# Main function to process multiple comparison groups
main <- function() {
  # Define the file paths for each comparison group
  comparisons <- list(
    list(limma_file = "DEGs_LIMMA_B_vs_C.csv", 
         deseq2_file = "*_DEGsUsingDESeq2_B_vs_C.csv", 
         output_meta = "meta_B_vs_C.csv", 
         comparison = "Bacterial vs Control"),
    list(limma_file = "DEGs_LIMMA_V_vs_C.csv", 
         deseq2_file = "*_DEGsUsingDESeq2_V_vs_C.csv", 
         output_meta = "meta_V_vs_C.csv", 
         comparison = "Viral vs Control"),
    list(limma_file = "DEGs_LIMMA_B_vs_V.csv", 
         deseq2_file = "*_DEGsUsingDESeq2_B_vs_V.csv", 
         output_meta = "meta_B_vs_V.csv", 
         comparison = "Bacterial vs Viral")
  )

  # Loop through each comparison and run the analysis
  for (comp in comparisons) {
    # Find the actual DESeq2 file by replacing the wildcard (*) with a prefix
    deseq2_file <- list.files(pattern = gsub("\\*", "", comp$deseq2_file))

    if (length(deseq2_file) == 0) {
      cat(paste0("No file matching pattern: ", comp$deseq2_file, "\n"))
      next
    }

    # Compare LIMMA and DESeq2 results and get merged data
    merged_data <- compare_DEGs(comp$limma_file, deseq2_file[1])

    # Perform Fisher's meta-analysis on the merged data
    perform_meta_analysis(merged_data, comp$output_meta)
  }
}

# Run the main function
main()
