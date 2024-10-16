#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
DESCRIPTION:
------------
This script combines two distinct functionalities:
1. Filtering out genes with low counts based on the median across samples.
2. Mapping Ensembl gene IDs to their respective gene symbols using a GTF file.

PART 1:
-------
Filter out genes from the count table whose median count is below a specified threshold (default: 10).
It processes a tab-delimited count table file and prints only rows where the median count is above the threshold.

PART 2:
-------
Maps Ensembl gene IDs to gene symbols using an Ensembl GTF annotation file. The output is a new count table where Ensembl IDs are replaced with gene symbols.

OUTPUT:
-------
1. Filtered count table (Part 1).
2. Count table with Ensembl IDs replaced by gene symbols (Part 2).

USAGE:
------
For Part 1:
python 07_maptoEnsemble_filterLowCounts.py filter <countData.txt> > <countDataFilt10Counts.txt>

For Part 2:
python 07_maptoEnsemble_filterLowCounts.py map <countDataFilt10Counts.txt> <Homo_sapiens.GRCh38.109.gtf> > <countDataFilt10CountsGeneSymbols.txt>
"""

import os
import sys
import glob
import statistics
from optparse import OptionParser

# General usage string
usage_string = "USAGE:\n1. For filtering: python %s filter <countData.txt>\n2. For mapping: python %s map <countDataFilt10Counts.txt> <GTFfile.gtf>" % (sys.argv[0], sys.argv[0])


def main():
    # Parse the command line arguments to determine which part of the script to execute
    if len(sys.argv) < 2:
        print(usage_string)
        sys.exit(1)
    
    mode = sys.argv[1]
    
    if mode == "filter":
        filter_low_counts()
    elif mode == "map":
        map_degs_to_gene_symbols()
    else:
        print(usage_string)
        sys.exit(1)


# PART 1: Filter out genes with low counts
def filter_low_counts():
    """
    Filters out genes from the count data file whose median count is below the threshold (default 10).
    """

    try:
        rawCountFileName = sys.argv[2]
    except IndexError:
        print("USAGE: python %s filter <countData.txt>" % sys.argv[0])
        sys.exit(1)

    # Open and process the count data file
    rawCounts = open(rawCountFileName)
    header = rawCounts.readline().rstrip("\n")
    print(header)

    # Filter rows where the median count across samples is >= 10
    for line in rawCounts:
        line = line.rstrip("\n")
        fields = line.split("\t")
        geneId, counts = fields[0], [int(x) for x in fields[1:]]
        if statistics.median(counts) >= 10:
            print(line)

    rawCounts.close()


# PART 2: Map Ensembl gene IDs to gene symbols
def map_degs_to_gene_symbols():
    """
    Replaces Ensembl gene IDs in the count table with gene symbols using a GTF annotation file.
    """
    try:
        geneFileName, ensemblGTFFileName = sys.argv[2], sys.argv[3]
    except ValueError:
        print("USAGE: python %s map <countDataFilt10Counts.txt> <GTFfile.gtf>" % sys.argv[0])
        sys.exit(1)

    # Load the Ensembl GTF to Gene Symbol mapping
    ensembl = get_mapping(ensemblGTFFileName)

    # Open the gene count file and replace Ensembl IDs with gene symbols
    geneFile = open(geneFileName)
    header = geneFile.readline().rstrip("\n").split("\t")
    header[0] = "Gene_ID"
    print("\t".join(header))

    for line in geneFile:
        line = line.rstrip("\n")
        fields = line.split("\t")
        ensemblId, exprvals = fields[0], fields[1:]
        if ensemblId in ensembl:
            print("\t".join([ensembl[ensemblId]] + exprvals))

    geneFile.close()

# Helper function to create a mapping from Ensembl IDs to gene symbols using the GTF file
def get_mapping(fileName):
    """
    Parses the GTF file and returns a dictionary mapping Ensembl gene IDs to gene symbols.
    """
    annoMap = {}
    with open(fileName) as ensembl:
        for line in ensembl:
            line = line.rstrip("\n")
            fields = line.split("\t")
            if fields[2] == "gene":
                annotations = {kv.split(" ")[0]: kv.split(" ")[1].strip('"') for kv in fields[8].split("; ") if " " in kv}
                geneId = annotations.get("gene_id", "")
                geneName = annotations.get("gene_name", "")
                if geneId and geneName:
                    annoMap[geneId] = geneName
    return annoMap

# Helper function to convert a list into a string with a separator
def convert_list_to_string(lst, sep):
    return sep.join(map(str, lst))

# Entry point of the script
if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        pass