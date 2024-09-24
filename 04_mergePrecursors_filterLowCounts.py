#!/usr/bin/env python3

import os
import sys
import glob
import statistics
import pandas as pd

# This script processes miRNA count data from raw files and generates a unified counts matrix.
# It performs the following tasks:
# 1. Merges read counts from different files for each miRNA.
# 2. Transforms the count data into a matrix format with miRNAs as rows and samples as columns.
# 3. Filters out miRNAs with median counts below a specified threshold (10).
# 
# USAGE: python script_name.py <input_directory> <output_directory> <clinical_ids_file>

def main():
    """
    Main function to process miRNA counts data. Reads input files, generates a count table, 
    transforms it into a matrix format, and filters out low counts.
    """
    if len(sys.argv) != 4:
        print("USAGE: python script_name.py <input_directory> <output_directory> <clinical_ids_file>")
        sys.exit(1)
    
    rCountsDirName = sys.argv[1]  # Directory with input *.csv files
    outputDirName = sys.argv[2]   # Directory for the output file
    clinIDsFile = sys.argv[3]     # File with clinical IDs

    # Generate count table
    countTableFileName = os.path.join(outputDirName, "countTable.txt")
    with open(countTableFileName, "w") as outputFile1:
        print("\t".join(["miRNA", "Sample_Id", "Count"]), file=outputFile1)
        
        countsMap = {}  # Dictionary to store miRNA counts
        clinIDs = getClinIDs(clinIDsFile)  # Load clinical IDs
        allFiles = glob.glob(os.path.join(rCountsDirName, '*.csv'))  # Find all CSV files

        # Process each file
        for i, aFile in enumerate(sorted(allFiles)):
            if os.path.basename(aFile).startswith("miRNAs"):
                sampleName = os.path.basename(aFile).split(".")[0].split("_")[-1]
                mergeCounts(aFile, countsMap, sampleName)  # Merge counts from file

        # Get sample IDs and condense counts
        condMap = condMiRs(countsMap)

        # Write condensed results to output file
        for geneIdsample in sorted(condMap):
            geneId, sample = geneIdsample.split("|")[0], geneIdsample.split("|")[1]
            countLst = condMap[geneIdsample]
            meanCount = statistics.mean(countLst)
            if sample in clinIDs:
                sample2 = clinIDs[sample]
                print("\t".join([str(x) for x in [geneId, sample2, meanCount]]), file=outputFile1)

    # Transform count table into a matrix format
    transform_to_matrix(countTableFileName)

    # Filter out low counts with a fixed threshold of 10
    filter_out_low_counts(countTableFileName, 10.0)

def getClinIDs(inFile):
    """
    Reads clinical IDs from a file and returns a dictionary mapping sample names to clinical IDs.
    
    :param inFile: Path to the file containing sample names and clinical IDs.
    :return: Dictionary where keys are sample names and values are clinical IDs.
    """
    iMap = {}
    with open(inFile) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            sample_name, clin_id = fields[0], fields[1]
            iMap[sample_name] = clin_id
    return iMap

def mergeCounts(inFile, cMap, sName):
    """
    Merges miRNA counts from an input file into a dictionary.
    
    :param inFile: Path to the input file containing miRNA counts.
    :param cMap: Dictionary to store the miRNA counts.
    :param sName: Sample name to be used as the key in the dictionary.
    """
    with open(inFile) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            miRNA, read_count, precursor = fields[0], fields[1], fields[2]
            geneId = miRNA + "|" + precursor
            fCounts = read_count
            fillMap(cMap, geneId, (sName, fCounts))

def getSampleIds(counts):
    """
    Extracts and returns a sorted list of sample IDs from the counts dictionary.
    
    :param counts: Dictionary containing miRNA counts.
    :return: Sorted list of sample IDs.
    """
    ids = set()
    for gene in counts:
        for (sample, reads) in sorted(counts[gene]):
            ids.add(sample)
    return sorted(list(ids))

def condMiRs(inMap):
    """
    Condenses miRNA counts by combining counts for the same miRNA across different samples.
    
    :param inMap: Dictionary containing miRNA counts.
    :return: Dictionary with condensed miRNA counts.
    """
    oMap = {}
    for key in inMap:
        key1 = key.split("|")[0]
        for (el1, el2) in inMap[key]:
            fillMap(oMap, key1 + "|" + el1, float(el2))
    return oMap

def fillMap(aMap, key, val):
    """
    Adds a value to a list associated with a key in a dictionary. Creates a new list if the key is not present.
    
    :param aMap: Dictionary to update.
    :param key: Key to update in the dictionary.
    :param val: Value to add to the list associated with the key.
    """
    aMap.setdefault(key, []).append(val)

def transform_to_matrix(rawCountFileName):
    """
    Transforms the count table into a matrix format where samples are columns and miRNAs are rows.
    
    :param rawCountFileName: Path to the input count table file.
    """
    dat = pd.read_csv(rawCountFileName, sep="\t", dtype=str)
    sorted_dat = dat.sort_values("Sample_Id")
    mat = sorted_dat.pivot(index="miRNA", columns="Sample_Id", values="Count")
    mat.to_csv(rawCountFileName.replace("countTable.txt", "countData.txt"), sep="\t")

def filter_out_low_counts(rawCountFileName, threshold):
    """
    Filters out miRNAs with median counts below a specified threshold.
    
    :param rawCountFileName: Path to the input count table file.
    :param threshold: Minimum median count required for miRNAs to be included.
    """
    with open(rawCountFileName) as rawCounts:
        header = rawCounts.readline().rstrip("\n")
        print(header)
        for line in rawCounts:
            line = line.rstrip("\n")
            fields = line.split("\t")
            geneId, rCounts = fields[0], [float(x) for x in fields[1:]]
            if statistics.median(rCounts) >= threshold:
                print(line)

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        pass
