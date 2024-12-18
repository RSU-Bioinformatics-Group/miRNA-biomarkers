#!/bin/bash

# This script performs the entire miRDeep2 analysis pipeline:
# 1. Maps reads using Bowtie.
# 2. Prepares files for miRDeep2.
# 3. Runs miRDeep2 to analyze miRNAs.
#
# Usage:
# bash 03_runmirDeep2.sh <input_directory> <output_directory> <bowtie_index_path> <reference_genome> <mature_mirna_file> <hairpin_mirna_file>

# Check if the correct number of arguments are provided
if [ "$#" -ne 6 ]; then
    echo "Usage: $0 <input_directory> <output_directory> <bowtie_index_path> <reference_genome> <mature_mirna_file> <hairpin_mirna_file>"
    exit 1
fi

# Assign command-line arguments to variables
INPATH=$1               # Path to input FASTQ files
OUTPATH=$2              # Path to output directory for mapped reads
BOWTIE_INDEX_PATH=$3    # Path to Bowtie index directory
REFERENCE_GENOME=$4     # Path to reference genome FASTA file
MATURE_MIRNA_FILE=$5    # Path to mature miRNA reference file which was previously modified as indicated in README
HAIRPIN_MIRNA_FILE=$6   # Path to hairpin miRNA reference file which was previously modified as indicated in README

# Step 1: Create output directory if it doesn't exist
mkdir -p $OUTPATH
chmod a+rwx $OUTPATH

# Step 2: Run Bowtie mapper on all FASTQ files in the input directory
for file in $INPATH/*.{fq,fastq,fq.gz,fastq.gz}; do
    bname=$(basename $file '.fq')
    echo "Running Bowtie mapper on $file..."
    mapper.pl $file -e -h -i -j -k AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA -l 18 -m -p $BOWTIE_INDEX_PATH -s $OUTPATH/$bname.fa -t $OUTPATH/$bname._vs_hsa.arf -v -o 16
done

# Step 3: Run miRDeep2 on the processed files
files=("$OUTPATH"/*.fa)
for file in "${files[@]}"; do
    echo "Processing file: $file"
    bname=$(basename $file '.fa')
    echo "Running miRDeep2 on $file..."
    miRDeep2.pl $file $REFERENCE_GENOME $OUTPATH/$bname._vs_hsa.arf $MATURE_MIRNA_FILE none $HAIRPIN_MIRNA_FILE -t Human -P -v -z $bname
done

# Completion message
echo "miRDeep2 pipeline complete! Results are stored in $OUTPATH"
