//****************************************************************************
//* SUBWORKFLOW: MAPPING
//* Map fastq with minimap2, index a reference, convert to BAM and sort, then
//* index.
//****************************************************************************


include { CRAMINO                } from '../../modules/local/cramino'
include { ALFRED                 } from '../../modules/local/alfred'
include { SAMTOOLS_FLAGSTAT      } from '../../modules/local/samtools/samtools_flagstat'

workflow BAM_QC {

    take:
    ch_bam // channel: [ val(meta), [ bam ] ]
    ch_genome_fasta

    main:

    ch_versions = Channel.empty()

    CRAMINO ( ch_bam )
    ch_versions = ch_versions.mix(CRAMINO.out.versions.first())

    ALFRED ( ch_bam, ch_genome_fasta )
    ch_versions = ch_versions.mix(ALFRED.out.versions.first())

    SAMTOOLS_FLAGSTAT ( ch_bam )
    //ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions.first())

    emit:
    cramino_stats      = CRAMINO.out.cramino_stats // channel: [ val(meta), [ bam ] ]
    alfred_stats      = ALFRED.out.alfred_stats          // channel: [ val(meta), [ bai ] ]
    flagstat      = SAMTOOLS_FLAGSTAT.out.flagstat          // channel: [ val(meta), [ csi ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}
