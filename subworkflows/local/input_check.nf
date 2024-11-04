/* Checks the input channel and creates channel map */

include { CHECK_SAMPLESHEET } from '../../modules/local/check_samplesheet.nf'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet

    main:
    CHECK_SAMPLESHEET( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { get_sample_info(it) }
        .set { ch_sample }
        /*.map { row -> tuple(row.sample, row.replicate,
        file(row.sequencing_summary_path, checkIfExists: false),
        file(row.read_path, checkIfExists: false)) }
        //.view()
        .set { ch_sample }
        ch_sample.view() */

    emit:
    ch_sample // [ sample, replicate, sequencing_summary_file, path_to_reads ]
    //ch_versions = ch_versions.mix(CHECK_SAMPLESHEET.out.versions.first())
}

// Create a meta map from the samplesheet
def get_sample_info(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id           = row.sample
    meta.replicate   = (row.replicate)
    meta.sequencing_summary = row.sequencing_summary_path
    //meta.fastq = row.read_path

    // add path(s) of the fastq file to the meta map
    def fastq_meta = []
    //if (!file(row.read_path).exists()) {
    //    exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    //}
    //if (meta.single_end) {
    fastq_meta = [ meta, [ file(row.read_path) ] ]
    //} else {
    //    if (!file(row.fastq_2).exists()) {
   //         exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
   //     }
    return fastq_meta
    //return meta
}
