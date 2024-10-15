#!/bin/bash

# This script builds a Bowtie index from a reference genome FASTA file.
# Usage:
# bash 02_build_bowtie-index.sh <genome_fasta> <index_name> <output_directory>

# Check if the correct number of arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <genome_fasta> <index_name> <output_directory>"
    exit 1
fi

# Assign command-line arguments to variables
GENOME_FASTA=$1       # Path to reference genome FASTA file
INDEX_NAME=$2         # Name for the Bowtie index
OUTDIR=$3            # Directory where the Bowtie index files will be saved

# Create output directory if it doesn't exist
mkdir -p $OUTDIR
chmod a+rwx $OUTDIR

# Copy the reference genome to the working directory
# This ensures that Bowtie can access the genome file.
cp $GENOME_FASTA $OUTDIR
GENOME_FILE=$(basename "$GENOME_FASTA")
cd $OUTDIR

# Build the Bowtie index
# This command creates an index from the reference genome FASTA file.
echo "Building Bowtie index..."
bowtie-build $GENOME_FILE $INDEX_NAME

# Completion message
echo "Bowtie index built successfully! Index files are stored in $OUTDIR"
