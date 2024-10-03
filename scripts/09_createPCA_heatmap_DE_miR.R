#!/usr/bin/env Rscript

# DESCRIPTION:
# This R script processes miRNA expression data by generating a heatmap and performing Principal Component Analysis (PCA).
# It reads an input file containing miRNA expression levels, scales the data, and creates a heatmap visualizing
# the Z-scores of selected miRNAs. The script also conducts PCA to analyze the distribution of samples and 
# generates PCA plots, distinguishing different sample groups. The resulting heatmap and PCA plots are saved to a PDF.

# USAGE:
# Rscript analysis_script.R --input <expression_data_file> --output <output_pdf_file>

library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(dplyr)
library(ggplot2)
library(optparse)
library(tidyverse)
library(ggrepel)
library(readxl)
library(cowplot)
library(ggfortify)
library(factoextra)

# Helper function to save plots to a PDF
save_plots <- function(plots, output_file) {
  pdf(output_file, width = 10, height = 7)
  for (plot in plots) {
    print(plot)
  }
  dev.off()
}

# Function to process input data
process_input_data <- function(input_file) {
  info_expr <- read.delim(input_file, check.names = FALSE)
  
  # Extract relevant feature names
  extracted_feature_names <- c('hsa-miR-136-5p', 'hsa-miR-513c-3p', 'hsa-miR-514a-5p', 'hsa-miR-514a-3p', 'hsa-miR-507')
  selected_rows <- info_expr[rownames(info_expr) %in% extracted_feature_names, ]
  
  # Extract sample groups from column names
  snames_2 <- colnames(info_expr)
  split_names <- strsplit(snames_2, "_")
  sample_group <- sapply(split_names, function(x) x[3])
  sample_group[is.na(sample_group)] <- "C"
  sample_group <- gsub("BS", "B", sample_group)
  
  list(selected_rows = selected_rows, sample_group = sample_group)
}

# Function to generate heatmap
generate_heatmap <- function(processed_data) {
  genes_expressed_matrix <- data.matrix(processed_data$selected_rows)
  gene_expr_matrix <- t(scale(t(genes_expressed_matrix)))
  
  ha = HeatmapAnnotation(Group = processed_data$sample_group, 
                         col = list(Group = c("B" = "darkred", "C" = "darkgreen", "V" = "lightblue")))
  
  ht_list <- Heatmap(gene_expr_matrix, show_row_names = TRUE, cluster_rows = FALSE, 
                     top_annotation = ha, name = "Z-score")
  
  draw(ht_list, merge_legend = TRUE, heatmap_legend_side = "left", annotation_legend_side = "left")
}

# Function to generate PCA plot
generate_pca_plot <- function(processed_data) {
  dat <- t(processed_data$selected_rows)
  dat.pca <- prcomp(dat, center = TRUE, scale = FALSE)
  
  # PCA of individuals
  pca_ind <- fviz_pca_ind(dat.pca, geom.ind = c("point", "text"), col.ind = processed_data$sample_group, 
                          addEllipses = TRUE, ellipse.level = 0.8, legend.title = "Groups") +
    labs(title = "PCA with miRNA markers", x = "PC1 (62.7%)", y = "PC2 (18.1%)") +
    theme(text = element_text(size = 20), axis.title = element_text(size = 20), 
          axis.text = element_text(size = 16), legend.position = 'bottom')
  
  # PCA of variables
  pca_var <- fviz_pca_var(dat.pca, geom.ind = c("point", "text"), col.var = "contrib", 
                          gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")) +
    labs(title = '', x = "PC1 (62.7%)", y = "PC2 (18.1%)") +
    theme(text = element_text(size = 18), axis.title = element_text(size = 20), 
          axis.text = element_text(size = 16), legend.position = 'bottom')
  
  list(pca_ind, pca_var)
}

# Main function to execute the analysis
main <- function() {
  option_list <- list(
    make_option(c("-i", "--input"), type = "character", help = "Input CSV file for heatmap/PCA"),
    make_option(c("-o", "--output"), type = "character", help = "Output PDF file")
  )
  
  opt <- parse_args(OptionParser(option_list = option_list))
  
  # Process input data
  processed_data <- process_input_data(opt$input)
  
  # Generate heatmap and PCA plots
  heatmap <- generate_heatmap(processed_data)
  pca_plots <- generate_pca_plot(processed_data)
  
  # Save all plots into a single PDF
  save_plots(c(heatmap, pca_plots), opt$output)
}

# Execute main function only if the script is run directly
if (interactive() == FALSE) {
  main()
}
