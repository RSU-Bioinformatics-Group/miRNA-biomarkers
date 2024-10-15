#!/bin/bash

# This script runs FastQC and MultiQC using Docker containers for quality control analysis of sequencing data.
# It accepts input and output directories as command-line arguments.
# FastQC is run on all .fq files in the input directory, and MultiQC generates a summary report based on the FastQC results.

# Usage:
#   bash 01_runFastQC_MultiQC.sh /path/to/input_directory /path/to/output_directory

# Check if the correct number of arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_directory> <output_directory>"
    exit 1
fi

# Assign the input and output directories from the command-line arguments
INPATH=$1
OUTPATH=$2

# Create output directory if it does not already exist
mkdir -p $OUTPATH
chmod a+rwx $OUTPATH

# Loop through all .fq files in the input directory and run FastQC on each
# -t 4 specifies 4 threads for FastQC execution (adjust as necessary)
for file in $INPATH/*.fq ;
do
    echo "Running FastQC on $file..."
    docker run --rm -v $INPATH:/data -v $OUTPATH:/output your-container-name fastqc /data/$(basename $file) --outdir /output -t 4
done

# After FastQC completes for all files, run MultiQC to aggregate the results
# -f option forces MultiQC to overwrite any existing reports in the output directory
echo "Running MultiQC on FastQC output (forcing report overwrite with -f)..."
docker run --rm -v $OUTPATH:/output your-container-name multiqc /output -o /output -f

# Print a completion message
echo "FastQC and MultiQC analysis complete!"
