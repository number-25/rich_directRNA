#!/bin/sh

# This script will download our target reference sequences

# hg38 genome - analysis set

rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz ./

# chrom sizes, aliase file and chromToUcsc

rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes ./

rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chromAlias.txt ./

wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/chromToUcsc

chmod +x chromToUcsc

# ensemble annotation 112

rsync -a -P rsync://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz | gunzip

# transcripts fasta sequence - gencode 46 matches ensemble 116
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.transcripts.fa.gz | gunzip

