# miRNA-biomarkers
This repository contains scripts that were used for miRNA and mRNA data pre-processing and miRNA diagnostic features prioritization in paediatric patients with bacterial or viral infections. All the analysis are optimized for execution within a High-Performance Computing (HPC) environment using provided Docker images. 
## Description
<img src="https://github.com/user-attachments/assets/a7d31e53-1c7b-4bcd-a4a3-8f43b4af1031" width="600">
A summary of workflow. Urine miRNA data analysis workflow marked in yellow, whole blood transcriptome data analysis workflow – in light red, workflow connecting both data types – in orange.
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
Before running the workflow please ensure:
1. You have all the necessary dependencies installed.
2.
3.
### Running the workflow:
