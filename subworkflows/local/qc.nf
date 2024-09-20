/* Checks the input channel and creates channel map */

include { NANOQ } from '../modules/local/check_samplesheet'

workflow INPUT_CHECK {

    take:
    samplesheet // file: /path/to/samplesheet

    main:

    CHECK_SAMPLESHEET( samplesheet )
        .samplesheet
        .splitCsv ( header:true, sep:',' )
        // .map { get_sample_info(it, params.genomes) }
        .map { it -> [ it[0], it[1], it[2], it[3] ] }
        .set { ch_sample } y
    ch_versions = Channel.empty()


    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:
    ch_sample // [ sample, replicate, sequencing_summary_file, path_to_reads ]

}


