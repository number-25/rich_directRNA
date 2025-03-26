//****************************************************************************
//* SUBWORKFLOW: MAPPING
//* Map fastq with minimap2, index a reference, convert to BAM and sort, then
//* index.
//****************************************************************************


include { CRAMINO               } from '../../modules/local/cramino'
include { ALFRED                } from '../../modules/local/alfred'
include { TRANSPOSE             } from '../../modules/local/transpose'
include { SAMTOOLS_FLAGSTAT     } from '../../modules/local/samtools/samtools_flagstat'
include { NGS_BITS              } from '../../modules/local/ngsbits'

workflow BAM_QC {

    take:
    skip_cramino
    skip_alfred
    skip_samtools_flagstat
    skip_ngs_bits
    ngs_bits_build
    ngs_bits_skip_contamination
    bam // channel: [ val(meta), [ bam ] ]
    genome_fasta
    cramino_min_length

    main:

    ch_versions = Channel.empty()

    // cramino
    if (skip_cramino != true) {
        CRAMINO( bam, cramino_min_length )
        ch_versions = ch_versions.mix(CRAMINO.out.versions.first())
    }

    // alfred
    if (skip_alfred != true) {
        ALFRED( bam, genome_fasta )
        alfred_stats = ALFRED.out.alfred_stats
        TRANSPOSE( alfred_stats )
        ch_versions = ch_versions.mix(ALFRED.out.versions.first())
        ch_versions = ch_versions.mix(TRANSPOSE.out.versions)
    }

    // cramino
    if (skip_samtools_flagstat != true) {
        SAMTOOLS_FLAGSTAT( bam )
    }

    if (skip_ngs_bits != true) {
        NGS_BITS(
        bam,
        genome_fasta,
        ngs_bits_build,
        ngs_bits_skip_contamination
        )
        ch_versions = ch_versions.mix(NGS_BITS.out.versions.first())
    }
    ch_ngs_bits = null

    emit:
    cramino_stats           = CRAMINO.out.cramino_stats // channel: [ val(meta), [ bam ] ]
    alfred_stats            = ALFRED.out.alfred_stats          // channel: [ val(meta), [ bai ] ]
    alfred_stats_transpose  = TRANSPOSE.out.alfred_stats_transposed          // channel: [ val(meta), [ bai ] ]
    flagstat                = SAMTOOLS_FLAGSTAT.out.flagstat          // channel: [ val(meta), [ csi ] ]
    //ngs_bits_stats      = NGS_BITS.out.qcML
    ngs_bits_stats          = ch_ngs_bits

    versions = ch_versions                     // channel: [ versions.yml ]
}
