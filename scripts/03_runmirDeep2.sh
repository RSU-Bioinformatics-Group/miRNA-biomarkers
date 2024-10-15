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
MATURE_MIRNA_FILE=$5    # Path to mature miRNA reference file
HAIRPIN_MIRNA_FILE=$6   # Path to hairpin miRNA reference file

# Step 1: Create output directory if it doesn't exist
mkdir -p $OUTPATH
chmod a+rwx $OUTPATH

# Step 2: Run Bowtie mapper on all FASTQ files in the input directory
for file in $INPATH/*.fq ; do
    bname=$(basename $file '.fq')
    echo "Running Bowtie mapper on $file..."
    mapper.pl $file -e -h -i -j -k AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA -l 18 -m -p $BOWTIE_INDEX_PATH -s $OUTPATH/$bname.fa -t $OUTPATH/$bname._vs_hsa.arf -v -o 16
done

# Step 3: Prepare files for miRDeep2 processing
# Setting locale environment variables
export LC_ALL=en_US.UTF-8
export LC_CTYPE=en_US.UTF-8 

# Clean up formatting issues in reference files
echo "Cleaning up reference files..."
perl -plane 's/\s+.+$//' < $REFERENCE_GENOME > ${REFERENCE_GENOME%.fa}_cleaned.fa
perl -plane 's/\s+.+$//' < $MATURE_MIRNA_FILE > ${MATURE_MIRNA_FILE%.fa}_cleaned.fa
perl -plane 's/\s+.+$//' < $HAIRPIN_MIRNA_FILE > ${HAIRPIN_MIRNA_FILE%.fa}_cleaned.fa

# Extract miRNAs from cleaned reference files
echo "Extracting miRNAs from reference files..."
extract_miRNAs.pl ${MATURE_MIRNA_FILE%.fa}_cleaned.fa hsa > ${MATURE_MIRNA_FILE%.fa}_extracted.fa
extract_miRNAs.pl ${HAIRPIN_MIRNA_FILE%.fa}_cleaned.fa hsa > ${HAIRPIN_MIRNA_FILE%.fa}_extracted.fa

# Step 4: Run miRDeep2 on the processed files
files=("$OUTPATH"/*.fa)
for file in "${files[@]}"; do
    echo "Processing file: $file"
    bname=$(basename $file '.fa')
    echo "Running miRDeep2 on $file..."
    miRDeep2.pl $file $REFERENCE_GENOME $OUTPATH/$bname._vs_hsa.arf ${MATURE_MIRNA_FILE%.fa}_extracted.fa none ${HAIRPIN_MIRNA_FILE%.fa}_extracted.fa -t Human -P -v -z $bname
done

# Completion message
echo "miRDeep2 pipeline complete! Results are stored in $OUTPATH"
