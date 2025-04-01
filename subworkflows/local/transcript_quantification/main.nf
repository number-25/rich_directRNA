include { TRANSIGNER    } from '../../../modules/local/transigner'
//include { OARFISH       } from '../../../modules/local/oarfish'

workflow TRANSCRIPT_QUANTIFICATION {

    take:
    ch_transcripts // channel: [ val(meta), [ bam ] ]
    ch_transcriptome_fasta

    main:

    if (!params.skip_transigner) {
        TRANSIGNER( ch_transcript, ch_transcriptome_fasta )
    }

    if (!params.skip_oarfish) {

    }

    ch_versions = Channel.empty()


    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:
    bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

