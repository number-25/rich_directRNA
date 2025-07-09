# syntax=docker/dockerfile:1
# R base image used to run bambu

FROM r-base:4.4.2

# Image metadata
LABEL maintainer="deancompbio@pm.me"

# Define version as build variable
ARG R_VER="4.4.2"
ARG BAMBU_VER="3.8.0"
ARG BIOCONDUCTOR_VER="3.2.0"

# Set some R environment variables
ENV R_PROFILE_USER=/.Rprofile
ENV R_ENVIRON_USER=/.Renviron

# Update base packages
RUN apt-get update

# Install Bioconductor Manager
RUN R -e "install.packages('BiocManager', quietly=TRUE, repos='http://cran.rstudio.com/')" && \
    R -e "BiocManager::install('bambu')" && \
    R -e "library(bambu)"

# Clean shit up
RUN apt-get clean all && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
