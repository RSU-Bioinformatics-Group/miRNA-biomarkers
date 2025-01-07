# miRNA-biomarkers
This repository contains scripts that were used for miRNA and mRNA data pre-processing and miRNA diagnostic features prioritization in paediatric patients with bacterial or viral infections. All the analysis are optimized for execution using provided Docker image. 
## Description
The repository contains all the scripts used pre-processing and analysing miRNA and mRNA datasets for diagnostic feature discovering. The steps include:
### miRNA workflow:
- Quality control of the raw data using FastQC and combining the outputs using MultiQC.
- miRNA mapping to human reference genome (hg38) and miRBase and quantification using miRDeep2.
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
1. You have Docker installed on your machine. Check this link for more information: [docker installation](https://docs.docker.com/engine/install/).
In case you want to use the workflow within an HPC environment, you should have Singularity installed there. Learn more about Singularity here: [intro to Singularity](https://bioinformaticsworkbook.org/Appendix/HPC/Containers/Intro_Singularity.html#gsc.tab=0).
2. You have cloned the Github repository with the necessary scripts. To do this, run (replace `.` with the directory you want to clone the repository to if necessary):
```
git clone https://github.com/tkiselova/miRNA-biomarkers .
```
3. For reproducibility and convinience the scripts for data preprocessing and differential expression can be run using a Docker container. The Dockerfile for building the image is provided in this repository. To build the image use the following command (replace `.` with the directory where the Dockerfile is located if necessary):
```
docker build -t mirna-biomarker-pipeline:latest .
```
In case you want to use Singularity, use this command. It will pull the Docker image and convert it into a Singularity format (*.sif):
```
sungularity pull /path/to/image
```
4. Create the directories for `.fq` input files and reference files in the same directory where the GitHub repository was cloned. Run:
```
mkdir output_data input_data references
chmod 775 /path/to/output_data
chmod 775 /path/to/references
chmod 775 /path/to/input_data
```
5. Prepare the necessary reference files. Run following commands (replace `/path/to/` with appropriate paths to files):
- for reference genome:
```
wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz
gunzip Homo_sapiens.GRCh38.113.gtf.gz Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz
```
- building bowtie index for miRDeep2:
```
docker run -v /path/to/references:/data2 mirna-biomarker-pipeline bash build_bowtie-index.sh /path/to/genome_fasta output_file_name 
```
If you are using Singularity, instead run:
```
singularity exec --bind /path/to/references:/data2 mirna-biomarker-pipeline bash build_bowtie-index.sh /path/to/genome_fasta output_file_name
```
- building STAR index for mRNA alignment:
```
docker run -v /path/to/references:/data2 mirna-biomarker-pipeline bash build_STAR-index.sh /path/to/genome_fasta /path/to/gtf_file
```
or, for Singularity:
```
singularity exec --bind /path/to/references:/data2 mirna-biomarker-pipeline bash build_STAR-index.sh /path/to/genome_fasta /path/to/gtf_file
```
- for miRNA sequences from miRBase. The `extract_miRNAs.pl` script from miRDeep2 prepares the reference files for further use for miRDeep2:
```
wget https://www.mirbase.org/download/CURRENT/hairpin.fa
wget https://www.mirbase.org/download/CURRENT/mature.fa
perl -plane 's/\s+.+$//' < mature.fa > mature_2.fa
perl -plane 's/\s+.+$//' < hairpin.fa > hairpin_2.fa
docker run -v /path/to/mature_2.fa:/data2 mirna-biomarker-pipeline perl extract_miRNAs.pl mature_2.fa hsa > miRBase_mature_hsa_v22_3.fa
docker run -v /path/to/hairpin_2.fa:/data2 mirna-biomarker-pipeline perl extract_miRNAs.pl hairpin_2.fa hsa > miRBase_hairpin_hsa_v22_3.fa
```
If you are using Singularity, instead run:
```
wget https://www.mirbase.org/download/CURRENT/hairpin.fa
wget https://www.mirbase.org/download/CURRENT/mature.fa
perl -plane 's/\s+.+$//' < mature.fa > mature_2.fa
perl -plane 's/\s+.+$//' < hairpin.fa > hairpin_2.fa
singularity exec --bind /path/to/mature_2.fa:/data2 mirna-biomarker-pipeline perl extract_miRNAs.pl mature_2.fa hsa > miRBase_mature_hsa_v22_3.fa
singularity exec --bind /path/to/hairpin_2.fa:/data2 mirna-biomarker-pipeline perl extract_miRNAs.pl hairpin_2.fa hsa > miRBase_hairpin_hsa_v22_3.fa
```
### Running the workflow:
The commands below describe the launching of the scripts (replace `/path/to/` with paths to appropriate directories). For more information on required inputs, see each of the scripts' description.
1. Running FastQC and MultiQC:
```
docker run -v /path/to/parent_dir:/data2 mirna-biomarker-pipeline bash 01_runFastQC_MultiQC.sh /path/to/input_directory /path/to/output_directory
```
#### miRNA workflow:
2. Running miRDeep2 for miRNA alignment and count generation:
```
docker run -v /path/to/parent_dir:/data2 mirna-biomarker-pipeline bash 03_runmirDeep2.sh path/to/input_directory path/to/output_directory path/to/bowtie_index path/to/reference_genome path/to/mature_mirna_file path/to/hairpin_mirna_file
```
3. Merging miRNA precursors and filter out low (median <10) counts:
```
docker run -v /path/to/parent_dir:/data2 mirna-biomarker-pipeline python3 04_mergePrecursors_filterLowCounts.py path/to/input_directory path/to/output_directory path/to/clinical_ids_file

```
4. Performing differential expression analysis using limma:
```
docker run -v /path/to/parent_dir:/data2 mirna-biomarker-pipeline Rscript 04_runLIMMA_DE.R /path/to/counts_table
```
5. Feature seelction with LASSO regularised linear regression:
```
docker run -v /path/to/parent_dir:/data2 mirna-biomarker-pipeline Rscript 05_runLASSO_DEG.R /path/to/upregulated_genes_file /path/to/downregulated_genes_file /path/to/normalized_expression_file
```
6. Feature selection with ANOVA:
```
docker run -v /path/to/parent_dir:/data2 mirna-biomarker-pipeline python3 06_runANOVA_DEG.py --normalized /path/to/normalized_data_file --up /path/to/upregulated_genes_file --down /path/to/downregulated_genes_file --output /path/to/output_csv_file
```
#### mRNA workflow
2. Running STAR for alignment:
```
docker run -v /path/to/parent_dir:/data2 mirna-biomarker-pipeline bash 02_runSTAR.sh path/to/star_index_path /path/to/forward_reads /path/to/reverse_reads /path/to/output_directory
```
3. Running Qualimap with MultiQC for alignment quality control:
```
docker run -v /path/to/parent_dir:/data2 mirna-biomarker-pipeline bash 03_runQualimap_MultiQC.sh /path/to/input_bam /path/to/output_directory /path/to/gtf_file
```
4. Merging files from different sequencing lanes of one sample with Samtools:
```
docker run -v /path/to/parent_dir:/data2 mirna-biomarker-pipeline bash 04_merge_bam_SAMtools.sh path/to/input_directory path/to/output_directory
```
5. Generating counts with featureCounts:
```
docker run -v /path/to/parent_dir:/data2 mirna-biomarker-pipeline Rscript 05_runFeaturecounts.R /path/to/input_bam /path/to/gtf_file /path/to/output_directory/status-file /path/to/output_directory/count-table
```
6. Merging count table:
```
docker run -v /path/to/parent_dir:/data2 mirna-biomarker-pipeline python3 10_renameGenerateCountTable.py /path/to/featureCounts_directory
```
7.1. Filtering out low counts:
```
docker run -v /path/to/parent_dir:/data2 mirna-biomarker-pipeline python3 07_maptoEnsemble_filterLowCounts.py filter /path/to/countData.txt > path/to/countDataFilt10Counts.txt
```
7.2. Mapping gene names to Ensemble:
```
docker run -v /path/to/parent_dir:/data2 mirna-biomarker-pipeline python3 map /path/to/countDataFilt10Counts.txt /path/to/gtf_file > path/to/countDataFilt10CountsGeneSymbols.txt
```
8.1. Performing differential expression analysis using limma:
```
docker run -v /path/to/parent_dir:/data2 mirna-biomarker-pipeline Rscript 08_1_runLIMMA_DE.R path/to/countDataFilt10CountsGeneSymbols.txt path/to/metadata_file /path/to/output_directory
```
8.2. Performing differential expression analysis using DESeq2:
```
docker run -v /path/to/parent_dir:/data2 mirna-biomarker-pipeline Rscript 08_2_run_DESeq2_DE.R path/to/countData.txt path/to/metadata_file path/to/output_directory
```
9. Unifying the results of differential expression analysis:
```
docker run -v /path/to/parent_dir:/data2 mirna-biomarker-pipeline Rscript 09_unify_DEGs.R
```
The scripts used for visualisations and correlation calculation `07_createVolcanoPlots.ipynb`, `08_createPCA_heatmap_DE_miR.ipynb` and `10_doCorrelation.ipynb` can be run in Jupyter Notebook.
## Visual summary of the workflow:
<div align="center">
<img src="https://github.com/user-attachments/assets/5b579158-6ccd-45fb-b7d7-428a95962f9e" width="600" alt="Summary of the workflow">
</div>
