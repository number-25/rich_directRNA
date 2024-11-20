# syntax=docker/dockerfile:1
# R base image used to run bambu

FROM r-base:4.4.2

# Image metadata
LABEL maintainer="deancompbio@pm.me"

# Define version as build variable
ARG R_VER="4.4.2"
ARG BAMBU_VER="3.8.0"
ARG BIOCONDUCTOR_VER="3.2.0"

# Update base packages
RUN apt-get update -y

# Install Bioconductor Manager
RUN R_PROFILE_USER=/.Rprofile && \
    R_ENVIRON_USER=/.Renviron && \
    Rscript -e "install.packages(c('BiocManager', 'dplyr', 'tidyr', 'xgboost', 'data.table', 'Rccp'), dependencies=TRUE, repos='http://cran.rstudio.com/')" && \
    Rscript -e "BiocManager::install(c('SummarizedExperiment', 'S4Vectors', 'BSgenome', 'IRanges', 'BiocGenerics', 'BiocParallel', 'GenomeInfoDb', 'GenomicAlignments', 'GenomicFeatures', 'GenomicRanges', 'Rsamtools', 'bambu'))"

RUN apt-get clean all && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
