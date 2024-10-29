# Uses official Bioconductor image as base image
FROM bioconductor/bioconductor_docker:devel AS bioconductor-base

LABEL maintainer="Tatjana Kiselova <tatjana.kiselova AT rsu.lv>" \
    base_image="bioconductor/bioconductor_docker:devel" \
    version="v0.1.0" \
    software="mirna-biomarker-pipeline" \
    about.summary="A pipeline used for miRNA and mRNA data preprocessing and miRNA urine biomarker distinguishing." \
    about.tags="Transcriptomics"

############ INIT ################

# Install system dependencies
RUN apt-get update && apt-get install -y \
    wget \
    gcc \
    g++ \
    make \
    openjdk-8-jdk \
    zlib1g-dev \
    bzip2 \
    build-essential \
    software-properties-common \
    cmake \
    libopenblas-dev \
    liblapack-dev \
    python3 \
    python3-pip \
    samtools \
    openjdk-8-jdk \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    gfortran \
    unzip \
    perl \
    cpanminus \
    bowtie \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install Python libraries
RUN pip3 install pandas multiqc

# Install Qualimap
RUN wget https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v2.3.zip && \
    unzip qualimap_v2.3.zip -d /opt/ && \
    chmod +x /opt/qualimap_v2.3/qualimap && \
    ln -s /opt/qualimap_v2.3/qualimap /usr/local/bin/qualimap && \
    rm qualimap_v2.3.zip

# Install FastQC
RUN wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && \
    unzip fastqc_v0.11.9.zip -d /opt/ && \
    chmod +x /opt/FastQC/fastqc && \
    ln -s /opt/FastQC/fastqc /usr/local/bin/fastqc && \
    rm fastqc_v0.11.9.zip

# Install STAR
RUN wget https://github.com/alexdobin/STAR/archive/refs/tags/2.7.10b.tar.gz && \
    tar -xzf 2.7.10b.tar.gz && \
    cd STAR-2.7.10b/source && \
    make STAR && \
    cp STAR /usr/local/bin/ && \
    cd ../.. && \
    rm -rf STAR-2.7.10b*

# Install miRDeep2
RUN git clone https://github.com/rajewsky-lab/mirdeep2.git /opt/mirdeep2
WORKDIR /opt/mirdeep2
RUN perl install.pl && \
    /bin/bash -c "source ~/.bashrc" && \
    perl install.pl
ENV PATH /opt/mirdeep2/bin:$PATH
ENV PERL5LIB /opt/mirdeep2/lib/perl5

# Install required R packages from Bioconductor and CRAN
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager', repos='http://cran.rstudio.com/')" && \
    R -e "BiocManager::install(c('DESeq2', 'IHW', 'limma', 'edgeR', 'Rsubread'), dependencies=TRUE)" && \
    R -e "install.packages(c('dplyr', 'metap'), dependencies=TRUE, repos='http://cran.rstudio.com/')"

# Set working directory
WORKDIR /data2

# Set environment variables for R
ENV R_HOME=/usr/lib/R

# Set default command
CMD ["bash"]

############ END ################

