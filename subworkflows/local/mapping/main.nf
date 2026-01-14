//****************************************************************************
//* SUBWORKFLOW: MAPPING
//* Map fastq with minimap2, index a reference, convert to BAM and sort, then
//* index.
//****************************************************************************

include { MINIMAP2_ALIGN     } from '../../../modules/local/minimap2/align'
include { SAMTOOLS_SORT      } from '../../../modules/local/samtools/sort'
include { SAMTOOLS_INDEX     } from '../../../modules/local/samtools/index'

workflow MAPPING {

    take:
    ch_sample // channel: [ val(meta), [ fastqpath ] ]
    ch_genome_minimap2_idx  // channel: [minimap2 genome index ]
    ch_transcriptome_minimap2_idx // channel [minimap2 transcriptome index ]

    main:

    ch_versions = Channel.empty()

    MINIMAP2_ALIGN ( ch_sample, ch_genome_minimap2_idx )
        ch_sample_sam   = MINIMAP2_ALIGN.out.sam
        ch_versions     = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

    SAMTOOLS_SORT ( ch_sample_sam )
        ch_sample_bam   = SAMTOOLS_SORT.out.bam
        ch_versions     = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )
        ch_sample_bam_idx   = SAMTOOLS_INDEX.out.bai
        ch_versions         = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:
    bam      = ch_sample_bam           // channel: [ val(meta), [ bam ] ]
    bai      = ch_sample_bam_idx          // channel: [ val(meta), [ bai ] ]
    versions = ch_versions                     // channel: [ versions.yml ]
}
