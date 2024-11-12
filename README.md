# miRNA-biomarkers
This repository contains scripts that were used for miRNA and mRNA data pre-processing and miRNA diagnostic features prioritization in paediatric patients with bacterial or viral infections. All the analysis are optimized for execution within a High-Performance Computing (HPC) environment using provided Docker images. 
## Description
The repository contains all the scripts used pre-processing and analysing miRNA and mRNA datasets for diagnostic feature discovering. The steps include:
### miRNA workflow:
- Quality control of the raw data using FastQC and combining the outputs using MultiQC.
- miRNA mapping to human reference genome (add version) and miRBase and quantification using miRDeep2.
- Count table processing for differential expression analysis using Python pandas library.
- Normalisation and differential expression analysis using R packages edgeR and limma.
- Feature selection for important diagnostic miRNA using logistic regression and ANOVA implemented via R package glmnet and Python library scikit-learn.
- Image producting using various R libraries mentioned below.
### mRNA workflow:
- Quality control of the raw data using FastQC and combining the outputs using MultiQC.
- Mapping to human reference genome (add version) using STAR.
- Adding transcripts to genome features using R package featurecounts.
- Count table processing for differential expression analysis using Python pandas library.
- Normalisation and differential expression using R packages edgeR, limma and DESeq2.
- Image producting using various R libraries mentioned below.
## Usage: 
### Prerequisites:
Before running the scripts please ensure:
1. You have Docker installed on your machine. Check this link for more information: https://docs.docker.com/engine/install/
2. Prepare the necessary reference files. Run following commands:
- for reference genome:
```
wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz
gunzip Homo_sapiens.GRCh38.113.gtf.gz Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz
```
- for miRNA sequences from miRBase:
```
wget https://www.mirbase.org/download/CURRENT/hairpin.fa
wget https://www.mirbase.org/download/CURRENT/mature.fa
perl -plane 's/\s+.+$//' < mature.fa > mature_2.fa
perl -plane 's/\s+.+$//' < hairpin.fa > hairpin_2.fa
extract_miRNAs.pl mature_2.fa hsa > miRBase_mature_hsa_v22_3.fa
extract_miRNAs.pl mature_2.fa hsa > miRBase_hairpin_hsa_v22_3.fa
```


4.
### Using the workflow:


For reproducibility and convinience the scripts for data preprocessing and differential expression can be run using a Docker container. The Dockerfile for building the image is provided in this repository. To build the image use the following command:
docker build -t mirna-biomarker-pipeline:latest .

<img src="https://github.com/user-attachments/assets/a7d31e53-1c7b-4bcd-a4a3-8f43b4af1031" width="600">
A summary of the workflow. Urine miRNA data analysis workflow marked in yellow, whole blood transcriptome data analysis workflow – in light red, workflow connecting both data types – in orange.
