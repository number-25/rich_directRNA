//****************************************************************************
//* SUBWORKFLOW: MAPPING
//* Map fastq with minimap2, index a reference, convert to BAM and sort, then
//* index.
//****************************************************************************

include { MINIMAP2_INDEX     } from '../../modules/local/minimap2_index'
include { MINIMAP2_ALIGN     } from '../../modules/local/minimap2_align'
include { SAMTOOLS_SORT      } from '../../modules/local/samtools/samtools_sort'
include { SAMTOOLS_INDEX     } from '../../modules/local/samtools/samtools_index'

workflow MAPPING {

    take:
    ch_sample // channel: [ val(meta), [ fastqpath ] ]
    ch_genome_fasta // channel: [ genome reference path ]
    ch_genome_minimap2_idx  // channel: [minimap2 genome index ]

    main:

    ch_versions = Channel.empty()

/*
    if (!params.minimap2_genome_idx) {
        MINIMAP2_INDEX ( ch_genome_fasta )
            ch_genome_minimap2_idx = MINIMAP2_INDEX.out.index
            ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions.first())
    }
*/
    MINIMAP2_ALIGN ( ch_sample, ch_genome_minimap2_idx )
        ch_sample_sam = MINIMAP2_ALIGN.out.sam

    SAMTOOLS_SORT ( ch_sample_sam )
        ch_sample_bam = SAMTOOLS_SORT.out.bam
        ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )
        ch_sample_bam_idx = SAMTOOLS_INDEX.out.bai

    emit:
    // do I need the index for anything else?
    //bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]
    bam      = ch_sample_bam           // channel: [ val(meta), [ bam ] ]
    //bai      = SAMTOOLS_INDEX.out.bam           // channel: [ val(meta), [ bam ] ]
    bai      = ch_sample_bam_idx          // channel: [ val(meta), [ bai ] ]
    //csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

