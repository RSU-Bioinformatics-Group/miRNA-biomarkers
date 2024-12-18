#!/bin/bash

# This script builds a Bowtie index from a reference genome FASTA file.
# Usage:
# bash 02_build_bowtie-index.sh <genome_fasta> <output_name>

# Check if the correct number of arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <genome_fasta> <output_name>"
    exit 1
fi

# Assign command-line arguments to variables
GENOME_FASTA=$1       # Path to reference genome FASTA file
INDEX_NAME=$2         # Name for the Bowtie index

# Build the Bowtie index
# This command creates an index from the reference genome FASTA file.
echo "Building Bowtie index..."
bowtie-build $GENOME_FILE $INDEX_NAME

# Completion message
echo "Bowtie index built successfully!"
