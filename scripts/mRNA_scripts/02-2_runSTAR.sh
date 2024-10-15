#!/bin/bash

# This script runs STAR for mapping paired-end RNA-seq reads to a reference genome.
# Usage:
# bash 02-2_runSTAR.sh <genome_directory> <forward_reads_directory> <reverse_reads_directory> <output_directory>

# Check if the correct number of arguments are provided
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <genome_directory> <forward_reads_directory> <reverse_reads_directory> <output_directory>"
    exit 1
fi

# Assign command-line arguments to variables
REFPATH=$1          # Path to STAR genome directory
FPATH=$2            # Path to directory containing forward reads (paired-end)
RPATH=$3            # Path to directory containing reverse reads (paired-end)
OUTPATH=$4          # Output directory for STAR results

# Create output directory if it doesn't exist
mkdir -p $OUTPATH
chmod a+rwx $OUTPATH

# Loop through each forward read file and find its corresponding reverse read file
for ffile in $FPATH/*_1.fq.gz; do
    rfile="$RPATH/$(basename "$ffile" _1.fq.gz)_2.fq.gz"  # Matching reverse read file
    fsample="$(basename "$ffile" _1.fq.gz)"              # Sample name from forward file
    rsample="$(basename "$rfile" _2.fq.gz)"              # Sample name from reverse file

    # Check if the forward and reverse samples match
    if [ "$fsample" == "$rsample" ]; then
        echo "Matching pair found: Forward sample: $fsample, Reverse sample: $rsample"

        # Run STAR for paired-end mapping
        STAR --runThreadN 16 \                          # Number of threads to use
             --genomeDir $REFPATH \                     # Path to genome index directory
             --readFilesIn $ffile $rfile \              # Forward and reverse read files
             --readFilesCommand zcat \                  # Decompress .gz files on-the-fly
             --outFileNamePrefix $OUTPATH/$fsample. \   # Output file prefix
             --outSAMtype BAM SortedByCoordinate        # Output sorted BAM file

    else
        echo "Skipping pair: Forward sample: $fsample, Reverse sample: $rsample"
    fi
done

# Completion message
echo "STAR mapping complete! Results are stored in $OUTPATH"
