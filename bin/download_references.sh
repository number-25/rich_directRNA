#!/bin/sh

set -euox pipefail

# This script will download our target reference sequences

if [ -d ../data ] ; then
    cd ../data
else
    mkdir ../data \
    cd ../data
fi

### hg38 genome - analysis set

#rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz ./
if [ -f hg38.analysisSet.fa.gz ] ; then
    wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz
fi

### chrom sizes, aliase file and chromToUcsc

#rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes ./
if [ -f hg38.chrom.sizes ] ; then
    wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes
fi

#rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chromAlias.txt ./
if [ -f hg38.chromAlias.txt ] ; then
    wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chromAlias.txt
fi

if [ -f chromToUcsc ] ; then
    wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/chromToUcsc
    chmod a+x chromToUcsc
fi

### ensemble annotation 112

#rsync -a -P rsync://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz | gunzip
if [ -f Homo_sapiens.GRCh38.112.gtf.gz ] ; then
wget https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz
fi

#### transcripts fasta sequence - gencode 46 matches ensemble 112
if [ -f gencode.v46.transcripts.fa.gz ] ; then
    wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.transcripts.fa.gz
fi

### Convert the chromosome naming in the ensemble annotation to match the genome
### (from 1 to chr1)
zcat Homo_sapiens.GRCh38.112.gtf.gz | chromToUcsc -a hg38.chromAlias.txt > Homo_sapiens.GRCh38.112_aliased.gtf

### Get CAGE-refTSS data from RIKEN
if [ -f refTSS_v4.1_human_coordinate.hg38.bed.txt.gz ] ; then
    wget https://reftss.riken.jp/datafiles/current/human/refTSS_v4.1_human_coordinate.hg38.bed.txt.gz
fi

### Get Poly-A motifs
if [ -f mouse_and_human.polyA_motif.txt ] ; then
    wget https://raw.githubusercontent.com/ConesaLab/SQANTI3/master/data/polyA_motifs/mouse_and_human.polyA_motif.txt
fi

### Get Poly-A sites
if [ -f atlas.clusters.2.0.GRCh38.96.bed.gz ] ; then
    wget https://polyasite.unibas.ch/download/atlas/2.0/GRCh38.96/atlas.clusters.2.0.GRCh38.96.bed.gz
fi
# is going to need to be converted to have the correct chromosome naming

### Get intropolis files - hg19 to hg38 liftover
if [ -f intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified.gz ] ; then
    wget https://github.com/Magdoll/images_public/raw/master/SQANTI2_support_data/intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified.gz
fi

### Get 30way PhyloP file - BigWig format
if [ -f hg38.phyloP30way.bw ] ; then
    wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP30way/hg38.phyloP30way.bw
fi
