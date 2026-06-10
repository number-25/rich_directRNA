//****************************************************************************
//* SUBWORKFLOW: SQANTI - ALL
//* Run SQANTI modules
//****************************************************************************


include { SQANTI_QC     } from '../../../modules/local/cramino/cramino'
include { SQANTI_FILTER } from '../../../modules/local/alfred/qc'
include { SQANTI_RESCUE } from '../../../modules/local/alfred/qc'

workflow BAM_QC {

    take:
    skip_sqanti_qc
    skip_sqanti_filter
    skip_sqanti_rescue
    reconstruction_program

    main:

    ch_versions = channel.empty()

    // cramino
    if (skip_cramino != true) {
        CRAMINO( bam, cramino_min_length )
        ch_cramino  = CRAMINO.out.cramino_stats
        ch_versions = ch_versions.mix(CRAMINO.out.versions.first())
    } else {
        ch_cramino = null
    }

    // alfred
    if (skip_alfred != true) {
        ALFRED( bam, genome_fasta )
        ch_alfred_stats              = ALFRED.out.alfred_stats
        TRANSPOSE( ch_alfred_stats )
        ch_alfred_stats_transposed      = TRANSPOSE.out.alfred_stats_transposed
        ch_versions                     = ch_versions.mix(ALFRED.out.versions.first())
        ch_versions                     = ch_versions.mix(TRANSPOSE.out.versions)
    } else {
        ch_alfred_stats            = null
        ch_alfred_stats_transposed = null
    }

    // samtools flagstat
    if (skip_samtools_flagstat != true) {
        SAMTOOLS_FLAGSTAT( bam )
        ch_samtools_flagstat = SAMTOOLS_FLAGSTAT.out.flagstat
        ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions)
    } else {
        ch_samtools_flagstat = null
    }

    if (skip_ngs_bits != true) {
        NGS_BITS(
        mixed_bam,
        genome_fasta,
        ngs_bits_build,
        ngs_bits_skip_contamination
        )
        ch_ngs_bits = NGS_BITS.out.qcML
        ch_versions = ch_versions.mix(NGS_BITS.out.versions.first())

    } else {
        ch_ngs_bits = null
    }

    emit:
    cramino_stats           = ch_cramino // channel: [ val(meta), [ bam ] ]
    alfred_stats            = ch_alfred_stats          // channel: [ val(meta), [ bai ] ]
    alfred_stats_transpose  = ch_alfred_stats_transposed          // channel: [ val(meta), [ bai ] ]
    flagstat                = ch_samtools_flagstat           // channel: [ val(meta), [ csi ] ]
    ngs_bits_stats          = ch_ngs_bits

    versions = ch_versions                     // channel: [ versions.yml ]
}
