#!/usr/bin/env Rscript

# DESCRIPTION:
# This R script generates volcano plots for differential expression analysis results.
# It reads upregulated and downregulated gene lists for different comparisons, labels the significantly
# upregulated and downregulated genes, and highlights the top genes in each comparison. 
# The script creates volcano plots for each comparison (B vs C, V vs C, and B vs V) and combines
# them into a single PDF, with a common legend included at the bottom.

# USAGE:
# Rscript 08_createVolcanoPlots.R --output <output_pdf_file> --input_bc_up <B_vs_C_UP.xlsx> --input_bc_down <B_vs_C_DOWN.xlsx> --input_vc_up <V_vs_C_UP.xlsx> --input_vc_down <V_vs_C_DOWN.xlsx> --input_bv_up <B_vs_V_UP.xlsx> --input_bv_down <B_vs_V_DOWN.xlsx>

# Required Libraries
library(tidyverse)
library(RColorBrewer)
library(ggrepel)
library(cowplot)
library(optparse)
library(readxl)

# Function to read, merge and annotate gene list data
read_and_prepare_data <- function(up_file, down_file) {
  gene_list_up <- read_excel(up_file)
  gene_list_down <- read_excel(down_file)
  gene_list <- rbind(gene_list_up, gene_list_down)
  
  # Annotate differential expression
  gene_list$diffexpressed <- "NO"
  gene_list$diffexpressed[gene_list$logFC > 1 & gene_list$P.Value < 0.05] <- "UP"
  gene_list$diffexpressed[gene_list$logFC < -1 & gene_list$P.Value < 0.05] <- "DOWN"
  
  # Filter for significant genes
  filtered_genes <- gene_list[gene_list$P.Value < 0.05, ]
  filtered_genes <- filtered_genes[order(filtered_genes$logFC), ]
  
  # Identify top 3 upregulated and downregulated genes
  top_downregulated <- head(filtered_genes, 3)
  top_upregulated <- tail(filtered_genes, 3)
  top_genes <- rbind(top_downregulated, top_upregulated)
  
  # Label top genes
  gene_list$delabel <- ifelse(gene_list$ID %in% top_genes$ID, gene_list$ID, NA)
  
  return(gene_list)
}

# Function to generate a volcano plot
generate_volcano_plot <- function(gene_data, title) {
  ggplot(data = gene_data, aes(x = logFC, y = -log10(P.Value), col = diffexpressed, label = delabel)) +
    geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
    geom_point(size = 4) +
    scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"),
                       labels = c("Downregulated", "Not significant", "Upregulated")) +
    coord_cartesian(ylim = c(0, 6.5), xlim = c(-5, 5)) +
    labs(color = '', x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) +
    scale_x_continuous(breaks = seq(-10, 10, 1)) +
    ggtitle(title) +
    theme_minimal(base_size = 24) +
    guides(color = FALSE) +
    theme(plot.title = element_text(face = 'bold', hjust = 0.5, size = 32)) +
    geom_text_repel(max.overlaps = Inf, size = 8, color = 'black')
}

# Function to generate the common legend
generate_legend <- function(gene_data) {
  ggplot(data = gene_data, aes(x = logFC, y = -log10(P.Value), col = diffexpressed)) +
    theme(legend.position = 'bottom') +  
    geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
    geom_point(size = 2) +
    scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"),
                       labels = c("Downregulated", "Not significant", "Upregulated")) +
    coord_cartesian(ylim = c(0, 7), xlim = c(-5, 5)) +
    labs(color = '', x = expression("log"[2]*"FC"), y = '') +
    scale_x_continuous(breaks = seq(-10, 10, 2))
}

# Main function to execute the script
main <- function() {
  # Argument parsing
  option_list <- list(
    make_option("--output", type = "character", default = "volcano_plots.pdf", help = "Output PDF file [default: %default]"),
    make_option("--input_bc_up", type = "character", help = "B vs C upregulated genes file"),
    make_option("--input_bc_down", type = "character", help = "B vs C downregulated genes file"),
    make_option("--input_vc_up", type = "character", help = "V vs C upregulated genes file"),
    make_option("--input_vc_down", type = "character", help = "V vs C downregulated genes file"),
    make_option("--input_bv_up", type = "character", help = "B vs V upregulated genes file"),
    make_option("--input_bv_down", type = "character", help = "B vs V downregulated genes file")
  )
  
  opt_parser <- OptionParser(option_list = option_list)
  opt <- parse_args(opt_parser)
  
  # Read and prepare data for each comparison
  gene_list_BC <- read_and_prepare_data(opt$input_bc_up, opt$input_bc_down)
  gene_list_VC <- read_and_prepare_data(opt$input_vc_up, opt$input_vc_down)
  gene_list_BV <- read_and_prepare_data(opt$input_bv_up, opt$input_bv_down)
  
  # Generate volcano plots for each comparison
  plot_BC <- generate_volcano_plot(gene_list_BC, "B vs C")
  plot_VC <- generate_volcano_plot(gene_list_VC, "V vs C")
  plot_BV <- generate_volcano_plot(gene_list_BV, "B vs V")
  
  # Generate the common legend
  legend <- get_legend(generate_legend(gene_list_VC))
  
  # Combine the plots and save to PDF
  combined_plot <- plot_grid(plot_BC, plot_VC, plot_BV, nrow = 1, labels = c('A', 'B', 'C'), label_size = 28)
  final_plot <- plot_grid(combined_plot, legend, ncol = 1, rel_heights = c(1, .1))
  
  # Save output to PDF
  ggsave(opt$output, final_plot, width = 16, height = 8)
}

# Execute the script
main()
