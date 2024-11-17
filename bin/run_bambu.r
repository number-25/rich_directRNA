#!/usr/bin/env Rscript

################################################
################################################
## REQUIREMENTS                               ##
################################################
################################################

## TRANSCRIPT ISOFORM DISCOVERY AND QUANTIFICATION
    ## - ALIGNED READS IN BAM FILE FORMAT
    ## - GENOME SEQUENCE
    ## - ANNOTATION GTF FILE
    ## - THE PACKAGE BELOW NEEDS TO BE AVAILABLE TO LOAD WHEN RUNNING R

################################################
################################################
## LOAD LIBRARY                               ##
################################################
################################################
library(bambu)

################################################
################################################
## PARSE COMMAND-LINE PARAMETERS              ##
################################################
################################################
args = commandArgs(trailingOnly=TRUE)

output_tag     <- strsplit(grep('--tag*', args, value = TRUE), split = '=')[[1]][[2]]
print(output_tag)
ncore          <- strsplit(grep('--ncore*', args, value = TRUE), split = '=')[[1]][[2]]
print(ncore)
genomeseq      <- strsplit(grep('--fasta*', args, value = TRUE), split = '=')[[1]][[2]]
print(genomeseq)
genomeSequence <- Rsamtools::FaFile(genomeseq)
Rsamtools::indexFa(genomeseq)
annot_gtf      <- strsplit(grep('--annotation*', args, value = TRUE), split = '=')[[1]][[2]]
print(annot_gtf)
readlist       <- args[5:length(args)]
print(readlist)


################################################
################################################
## RUN BAMBU                                  ##
################################################
################################################
grlist <- prepareAnnotations(annot_gtf)
se     <- bambu(reads = readlist, annotations = grlist, genome = genomeSequence, ncore = ncore, verbose = TRUE)
writeBambuOutput(se, output_tag)

