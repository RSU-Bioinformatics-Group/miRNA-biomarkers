# Dependencies of the pipeline:

# - FastQC/0.11.9          # from official page
# - MultiQC/1.24.1         # from PyPI
#      @requires: python3
#      @requires: python3-pip
# - miRDeep2               # from github
#      @requires: bowtie/1.X
#      @requires: Vienna and RNAfold
#      @requires: SQUID
#      @requires: randfold
#      @requires: PDF::API2 perl package 
# - STAR/2.7.10b           # from github       
# - qualimap/2.3           # from bitbucket.org
# - R/4.4.1
#    - Bioconductor packages DESeq2, edgeR, limma # from bioconductor.org
#    - glmnet                                     # from cran.rstudio.com

LABEL maintainer="Tatjana Kiselova <tatjana.kiselova AT rsu.lv>" \
    base_image="bioconductor/bioconductor_docker:devel" \
    version="v0.1.0"   \
    software="mirna-biomarker-pipeline" \
    about.summary="A pipeline used for miRNA and mRNA data preprocessing and miRNA urine biomarker distinguishing." \
    about.home="" \
    about.tags="Transcriptomics"


############INIT################

# Uses official Bioconductor image as base image
FROM bioconductor/bioconductor_docker:devel AS bioconductor-base

# Install additional tools and libraries in the Bioconductor image neseccary for running other tools
RUN apt-get update && apt-get install -y \
    wget gcc g++ make openjdk-8-jdk python3 python3-pip zlib1g-dev \
    bzip2 build-essential software-properties-common \
    perl cpanminus libopenblas-dev liblapack-dev cmake && \
    add-apt-repository ppa:ubuntu-toolchain-r/test && \
    apt-get update && \
    apt-get install -y g++-11 && \
    update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-11 60 && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Install MultiQC
RUN pip install multiqc

# Install qualimap
RUN wget https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v2.3.zip && \
    unzip qualimap_v2.3.zip -d /opt/ && \
    chmod +x /opt/qualimap_v2.3/qualimap && \
    ln -s /opt/qualimap_v2.3/qualimap /usr/local/bin/qualimap && \
    rm qualimap_v2.3.zip

# Install required R packages from Bioconductor and glmnet
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager', repos='http://cran.rstudio.com/')" && \
    R -e "BiocManager::install(c('DESeq2', 'limma', 'edgeR'), dependencies=TRUE)" && \
    R -e "install.packages('glmnet', dependencies=TRUE, repos='http://cran.rstudio.com/')"

# Install STAR
RUN wget https://github.com/alexdobin/STAR/archive/refs/tags/2.7.10b.tar.gz && \
    tar -xzf 2.7.10b.tar.gz && \
    cd STAR-2.7.10b/source && \
    make STAR && \
    cp STAR /usr/local/bin/ && \
    cd ../.. && \
    rm -rf STAR-2.7.10b*

# Install FastQC
ADD http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip /usr/local/bin/
RUN cd /usr/local/bin/ && \
    unzip fastqc_v0.11.8.zip && \
    chmod 755 FastQC/fastqc && \
    ln -s /usr/local/bin/FastQC/fastqc /usr/local/bin/fastqc && \
    rm fastqc_v0.11.8.zip

# Install miRDeep2
RUN git clone https://github.com/rajewsky-lab/mirdeep2.git /opt/mirdeep2
WORKDIR /opt/mirdeep2
RUN perl install.pl && \
    /bin/bash -c "source ~/.bashrc" && \
    perl install.pl
ENV PATH /opt/mirdeep2/bin:$PATH
ENV PERL5LIB /opt/mirdeep2/lib/perl5

# Set environment variables for R
ENV R_HOME=/usr/lib/R

# Set working directory
WORKDIR /data2

# Set default command
CMD ["bash"]
