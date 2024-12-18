#!/bin/bash

# This script generates a STAR genome index from a reference genome FASTA file and a GTF annotation file.
# Usage:
# bash build_STAR-index.sh <genome_fasta> <gtf_file>

# Check if the correct number of arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <genome_fasta> <gtf_file>"
    exit 1
fi

# Assign command-line arguments to variables
GENOME_FASTA=$1       # Path to reference genome FASTA file
GTF_FILE=$2           # Path to GTF annotation file

# Generate the STAR genome index
echo "Generating STAR genome index..."
STAR --runMode genomeGenerate \
     --genomeDir . \
     --genomeFastaFiles $GENOME_FASTA \
     --sjdbGTFfile $GTF_FILE \
     --sjdbOverhang 149 \
     --outFileNamePrefix star_index_

# Set read-write permissions for the output files
chmod a+rw *

# Completion message
echo "STAR genome index generation complete!"
