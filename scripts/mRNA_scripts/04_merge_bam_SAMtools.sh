#!/bin/bash

# This script uses SAMtools to merge BAM files from multiple subdirectories into single BAM files.
# It accepts the base input directory and output directory as command-line arguments.
# Usage: bash your_script_name.sh <input_directory> <output_directory>
# Example: bash merge_bam_SAMtools.sh path/to/input/directory path/to/output/directory 

# Check if the correct number of arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_directory> <output_directory>"
    exit 1
fi

# Assign command-line arguments to variables
BASE_INPUT_DIR=$1   # Base directory for input subdirectories containing BAM files
OUTPATH=$2          # Output directory for merged BAM files

# Create output directory if it doesn't exist
mkdir -p "$OUTPATH"
chmod a+rwx "$OUTPATH"  # Ensure the output directory has read, write, and execute permissions

# Load SAMtools module
module load bio/samtools/1.10

# Loop through all subdirectories in the base input directory
for subdir in "$BASE_INPUT_DIR"/*; do
    if [ -d "$subdir" ]; then  # Ensure the item is a directory
        sample_name=$(basename "$subdir")  # Extract the directory name for the output file
        
        # Check if there are any BAM files in the subdirectory
        bam_files=("$subdir"/*.bam)
        if [ -n "${bam_files[0]}" ]; then  # Only process if BAM files are found
            echo "Processing directory: $subdir"
            
            # Merge BAM files in the subdirectory using SAMtools
            samtools merge -c -f "$OUTPATH/$sample_name.bam" "$subdir"/*.bam
            
            # Check if the merge was successful
            if [ $? -eq 0 ]; then
                echo "Successfully merged files for $sample_name."
            else
                echo "Error merging files for $sample_name."
            fi
        else
            echo "No BAM files found in $subdir. Skipping."
        fi
    else
        echo "$subdir is not a directory. Skipping."
    fi
done

# Completion message
echo "SAMtools merging process complete! Merged BAM files are stored in $OUTPATH."

