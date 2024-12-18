#!/bin/bash

# This script runs STAR for mapping paired-end RNA-seq reads to a reference genome.
# Usage:
# bash 02_runSTAR.sh <genome_directory> <forward_reads_directory> <reverse_reads_directory> <output_directory>

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

# Loop through all possible forward read files
for ffile in $FPATH/*_{1,R1}.{fq,fastq}{,.gz}; do
    # Check if the forward file exists to avoid errors with the glob pattern
    if [ -e "$ffile" ]; then
        if [[ "$ffile" == *"_1.fq" ]] || [[ "$ffile" == *"_1.fastq" ]]; then
            rfile="$RPATH/$(basename "$ffile" _1.fq)_2.fq"
            rfile="${rfile%.fq}.fastq"  
        elif [[ "$ffile" == *"_1.fq.gz" ]] || [[ "$ffile" == *"_1.fastq.gz" ]]; then
            rfile="$RPATH/$(basename "$ffile" _1.fq.gz)_2.fq.gz"
            rfile="${rfile%.fq.gz}.fastq.gz"  
        elif [[ "$ffile" == *"R1.fq" ]] || [[ "$ffile" == *"R1.fastq" ]]; then
            rfile="$RPATH/$(basename "$ffile" R1.fq)R2.fq"
            rfile="${rfile%.fq}.fastq" 
        elif [[ "$ffile" == *"R1.fq.gz" ]] || [[ "$ffile" == *"R1.fastq.gz" ]]; then
            rfile="$RPATH/$(basename "$ffile" R1.fq.gz)R2.fq.gz"
            rfile="${rfile%.fq.gz}.fastq.gz" 
        fi

        # Extract the sample name from the forward and reverse files
        fsample="$(basename "$ffile" | sed -e 's/_1\.fq//' -e 's/_1\.fastq//' -e 's/R1\.fq//' -e 's/R1\.fastq//' -e 's/\.gz//')"
        rsample="$(basename "$rfile" | sed -e 's/_2\.fq//' -e 's/_2\.fastq//' -e 's/R2\.fq//' -e 's/R2\.fastq//' -e 's/\.gz//')"

        # Check if the forward and reverse samples match
        if [ "$fsample" == "$rsample" ]; then
            echo "Matching pair found: Forward sample: $fsample, Reverse sample: $rsample"
            echo "Forward file: $ffile"
            echo "Reverse file: $rfile"
            # Run STAR for paired-end mapping
            STAR --runThreadN 16 \                          
             --genomeDir $REFPATH \                     
             --readFilesIn $ffile $rfile \              
             --readFilesCommand zcat \                  
             --outFileNamePrefix $OUTPATH/$fsample. \   
             --outSAMtype BAM SortedByCoordinate        
        else
            echo "Sample names do not match: Forward sample: $fsample, Reverse sample: $rsample"
        fi
    else
        echo "No matching forward files found in $FPATH"
        break
    fi
done

# Completion message
echo "STAR mapping complete! Results are stored in $OUTPATH"
