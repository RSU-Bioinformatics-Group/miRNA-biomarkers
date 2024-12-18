#!/usr/bin/env python3

# DESCRIPTION:
# This Python script performs feature selection on gene expression data using SelectKBest (ANOVA F-test).
# It reads normalized gene expression data and lists of upregulated and downregulated genes.
# The script then selects common genes, applies feature selection, and outputs the selected features
# (those with a p-value < 0.05) to a CSV file.

# USAGE:
# python 07_runANOVA_DEG.py --normalized <normalized_data_file> --up <upregulated_genes_file> --down <downregulated_genes_file> --output <output_csv_file>

import pandas as pd
import numpy as np
from sklearn.feature_selection import SelectKBest, f_classif
import argparse

# Function to load data
def load_data(normalized_file, up_file, down_file):
    raw = pd.read_csv(normalized_file, index_col=0, sep='\t')
    DEG_UP = pd.read_excel(up_file)
    DEG_DOWN = pd.read_excel(down_file)
    
    DEG = pd.concat([DEG_UP, DEG_DOWN])
    DEG = DEG.sort_values(by='logFC', ascending=False)
    DEG.reset_index(inplace=True)
    
    return raw, DEG

# Function to select common genes
def select_common_genes(raw, DEG):
    common_genes = raw.index.intersection(DEG['ID'])
    miRNAs_selected = raw.loc[common_genes].transpose()
    
    return miRNAs_selected

# Function to perform feature selection
def feature_selection(miRNAs_selected, y, k=5):
    selector = SelectKBest(f_classif, k=k)
    X_new = selector.fit_transform(miRNAs_selected, y)
    
    return selector

# Function to save selected features
def save_results(selector, miRNAs_selected, output_file):
    miRNA_names = miRNAs_selected.columns.tolist()
    results = pd.DataFrame({
        'Importance': selector.scores_,
        'P-Value': selector.pvalues_
    }, index=miRNA_names)
    
    # Filter features with p-value < 0.05
    results_filtered = results[results['P-Value'] < 0.05]
    results_filtered = results_filtered.sort_values(by='P-Value', ascending=True)
    
    # Save to CSV
    results_filtered.to_csv(output_file)
    print(f"Results saved to {output_file}")

# Main function to orchestrate the workflow
def main():
    # Argument parsing
    parser = argparse.ArgumentParser(description='Feature selection on miRNA data.')
    
    parser.add_argument('--normalized', type=str, required=True, help='Path to the normalized TMM log2 CPM+1 data file')
    parser.add_argument('--up', type=str, required=True, help='Path to the DEGs upregulated gene list (xlsx)')
    parser.add_argument('--down', type=str, required=True, help='Path to the DEGs downregulated gene list (xlsx)')
    parser.add_argument('--output', type=str, required=True, help='Path to the output CSV file for selected features')
    
    args = parser.parse_args()
    
    # Load data
    raw, DEG = load_data(args.normalized, args.up, args.down)
    
    # Select common genes
    miRNAs_selected = select_common_genes(raw, DEG)
    
    # Labels for classification
    y = [1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0]
    
    # Perform feature selection
    selector = feature_selection(miRNAs_selected, y, k=5)
    
    # Save results to CSV
    save_results(selector, miRNAs_selected, args.output)

# Execute the script
if __name__ == "__main__":
    main()
