//****************************************************************************
//* SUBWORKFLOW: SQANTIQC
//* Run SQANTI QC on output of transcript reconstruction programs
//****************************************************************************

include { QC_SQANTI } from '../../modules/local/sqanti/sqanti_qc'

workflow SQANTI_QC {

    take:
    reconstructed_transcriptome
    annotation_gtf
    genome_fasta
    program

    main:

    ch_versions = Channel.empty()

    QC_SQANTI( reconstructed_transcriptome, annotation_gtf, genome_fasta, program)


    if (skip_cramino != true) {
        CRAMINO( bam, cramino_min_length )
        ch_versions = ch_versions.mix(CRAMINO.out.versions.first())
    }

    emit:
    cramino_stats           = CRAMINO.out.cramino_stats // channel: [ val(meta), [ bam ] ]
    alfred_stats            = ALFRED.out.alfred_stats          // channel: [ val(meta), [ bai ] ]
    alfred_stats_transpose  = TRANSPOSE.out.alfred_stats_transposed          // channel: [ val(meta), [ bai ] ]
    flagstat                = SAMTOOLS_FLAGSTAT.out.flagstat          // channel: [ val(meta), [ csi ] ]
    //ngs_bits_stats      = NGS_BITS.out.qcML
    ngs_bits_stats          = ch_ngs_bits

    versions = ch_versions                     // channel: [ versions.yml ]
}
