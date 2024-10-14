#!/usr/bin/env Rscript

# This script performs feature counting on RNA-seq data using Rsubread's featureCounts function.
# It takes a BAM file and a GTF gene annotation file as input and generates two output files:
# 1. A status file that includes the read counts.
# 2. A counts file with Gene IDs and their corresponding counts.
# 
# USAGE:
# Rscript your_script_name.R <bam_file> <gene_annotation_file> <output_file_status> <output_file_counts>
##########################################################################################################

# Load the Rsubread library 
library(Rsubread)

# Main function that processes command line arguments
main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  # Check if the correct number of arguments are provided
  if (length(argv) != 4) {
    stop("Usage: Rscript your_script_name.R <bam_file> <gene_annotation_file> <output_file_status> <output_file_counts>")
  }
  
  # Assign command line arguments to variables
  bamFileName <- argv[1]                       # Input BAM file
  geneAnnotationFile <- argv[2]                # Input gene annotation GTF file
  outFileName1 <- argv[3]                      # Output file for status
  outFileName2 <- argv[4]                      # Output file for counts
  
  # Perform feature counting using Rsubread's featureCounts function
  fCounts <- featureCounts(
    files = bamFileName,
    isPairedEnd = TRUE,
    annot.ext = geneAnnotationFile,          # Use gene annotation file from command line
    isGTFAnnotationFile = TRUE
  )
  
  # Rename columns for the status output
  colnames(fCounts$stat) <- c("Status", "ReadCount")
  
  # Write the status information to the output file
  write.table(fCounts$stat, outFileName1, sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Prepare the counts data frame for output
  colnames(fCounts$counts) <- sub("^([^.]*).*", "\\1", basename(bamFileName))  # Rename columns based on BAM file name
  GeneID <- rownames(fCounts$counts)  # Extract Gene IDs
  fCountsDF <- as.data.frame(cbind(GeneID, fCounts$counts))  # Combine Gene IDs with counts
  
  # Write the counts data frame to the second output file
  write.table(fCountsDF, outFileName2, sep = "\t", quote = FALSE, row.names = FALSE)
}

# Call the main function
main()
