//****************************************************************************
//* SUBWORKFLOW: MAPPING
//* Map fastq with minimap2, index a reference, convert to BAM and sort, then
//* index.
//****************************************************************************


include { CRAMINO               } from '../../../modules/local/cramino/cramino'
include { ALFRED                } from '../../../modules/local/alfred/qc'
include { TRANSPOSE             } from '../../../modules/local/alfred/transpose'
include { SAMTOOLS_FLAGSTAT     } from '../../../modules/local/samtools/flagstat'
include { NGS_BITS              } from '../../../modules/local/ngsbits/mapping_qc'

workflow BAM_QC {

    take:
    skip_cramino
    skip_alfred
    skip_samtools_flagstat
    skip_ngs_bits
    ngs_bits_build
    ngs_bits_skip_contamination
    bam // channel: [ val(meta), [ bam ] ]
    mixed_bam // bam channel with bam.bai index included
    genome_fasta
    transcriptome_fasta
    cramino_min_length

    main:

    ch_versions = Channel.empty()

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
        if (!params.transcriptome_mapping) {
            ALFRED( bam, genome_fasta )
            ch_alfred_stats              = ALFRED.out.alfred_stats
            TRANSPOSE( ch_alfred_stats )
            ch_alfred_stats_transposed      = TRANSPOSE.out.alfred_stats_transposed
            ch_versions                     = ch_versions.mix(ALFRED.out.versions.first())
            ch_versions                     = ch_versions.mix(TRANSPOSE.out.versions)
        } else {
            ALFRED( bam, transcriptome_fasta )
            ch_alfred_stats              = ALFRED.out.alfred_stats
            TRANSPOSE( ch_alfred_stats )
            ch_alfred_stats_transposed      = TRANSPOSE.out.alfred_stats_transposed
            ch_versions                     = ch_versions.mix(ALFRED.out.versions.first())
            ch_versions                     = ch_versions.mix(TRANSPOSE.out.versions)
        }
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
        if (!params.transcriptome_mapping) {
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
