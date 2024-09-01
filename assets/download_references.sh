#!/bin/sh

set -euox pipefail

# This script will download our target reference sequences

if [ -d ../data ] ; then
    cd ../data
else
    mkdir ../data \
    cd ../data
fi

# hg38 genome - analysis set

#rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz ./
#wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz

# chrom sizes, aliase file and chromToUcsc

#rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes ./
#wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes ./

#rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chromAlias.txt ./
#wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chromAlias.txt ./

wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chromAlias.txt ./

wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/chromToUcsc

chmod a+x chromToUcsc

# ensemble annotation 112

#rsync -a -P rsync://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz | gunzip
wget https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz

# transcripts fasta sequence - gencode 46 matches ensemble 112
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.transcripts.fa.gz

# Convert the chromosome naming in the ensemble annotation to match the genome
# (from 1 to chr1)
zcat Homo_sapiens.GRCh38.112.gtf.gz | chromtoUcsc -a hg38.chromAlias.txt > Homo_sapiens.GRCh38.112_aliased.gtf

# Get CAGE-refTSS data from RIKEN
wget https://reftss.riken.jp/datafiles/current/human/refTSS_v4.1_human_coordinate.hg38.bed.txt.gz

# Get Poly-A motifs
wget https://raw.githubusercontent.com/ConesaLab/SQANTI3/master/data/polyA_motifs/mouse_and_human.polyA_motif.txt

# Get Poly-A sites
wget https://polyasite.unibas.ch/download/atlas/2.0/GRCh38.96/atlas.clusters.2.0.GRCh38.96.bed.gz
# is going to need to be converted to have the correct chromosome naming

# Get intropolis files - hg19 to hg38 liftover
wget https://github.com/Magdoll/images_public/raw/master/SQANTI2_support_data/intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified.gz

# Get 30way PhyloP file - BigWig format
#wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP30way/hg38.phyloP30way.bw
