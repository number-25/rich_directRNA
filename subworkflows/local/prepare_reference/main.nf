//               https://github.com/nf-core/modules/tree/master/subworkflows
//
// Uncompress and prepare reference genome files
//

include { GUNZIP as GUNZIP_FASTA } from '../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GTF } from '../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_BED } from '../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_TRANSCRIPT_FASTA } from
'../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_POLYA } from '../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_CAGE } from '../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_INTROPOLIS } from '../modules/nf-core/gunzip'

// prepare indices for reference

include { MINIMAP2_INDEX as MINIMAP2_GENOME_INDEX } from '../modules/nf-core/minimap2/index'
include { MINIMAP2_INDEX as MINIMAP2_TRANSCRIPTOME_INDEX } from '../modules/nf-core/minimap2/index'
include { SAMTOOLS_FAIDX } from '../modules/nf-core/samtools/faidx'
include { CUSTOM_GETCHROMSIZES } from
'../modules/nf-core/custom/getchromsizes'
include { SAMTOOLS_INDEX } from
'../modules/nf-core/samtools/index'


// prepare additional files

//TO-DO make these modules
include { COLLAPSE_GTF } from '../modules/local/collapsegtf
include { GTF_TO_BED } from '../modules/local/gfx2bed' // gxf2bed module
include { BIGWIG_TO_WIG } from '../modules/local/bigwigtowig'
include { WIG_TO_BED } from '../modules/local/bedops'
include { JAFFAL_PREPARE_REF } from
'../modules/local/jaffal_prepare_ref'


include { SAMTOOLS_SORT      } from '../../../modules/nf-core/samtools/sort/main'

workflow PREPARE_REFERENCE {

    take:
    genome_fasta
    transcriptome_fasta
    annotation_gtf
    cage_bed
    polyA_bed
    intropolis_bed
    phylop_bigwig



    appris_bed?
    mane_select_bed?
    mane_clinical_bed?
    // TODO nf-core: edit input (take) channels
    ch_bam // channel: [ val(meta), [ bam ] ]

    main:

    ch_versions = Channel.empty()

    // TODO nf-core: substitute modules here for the modules of your subworkflow

    SAMTOOLS_SORT ( ch_bam )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:
    genome_fasta = ch_genome_fasta
    transcriptome_fasta = ch_transcriptome_fasta
    annotation_gtf =
    annotation_bed =
    collapsed_gtf =
    cage_bed =
    polyA_bed =
    jaffal_ref =
    ch_
    // TODO nf-core: edit emitted channels
    bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

https://raw.githubusercontent.com/nf-core/rnaseq/master/subworkflows/local/prepare_genome/main.nf
