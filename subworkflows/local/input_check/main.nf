/* Checks the input channel and creates channel map */

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet

    main:
    ch_samplesheet = Channel.fromPath(samplesheet, checkIfExists: true)
    ch_samplesheet
        .splitCsv ( header:true, sep:',' )
        .map { get_sample_info(it) }
        .set { ch_sample }

    emit:
    ch_sample // [ sample, replicate, path_to_reads ]
}

// Create a meta map from the samplesheet
def get_sample_info(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id           = row.sample
    meta.replicate    = row.replicate

    // add path(s) of the fastq file to the meta map
    def fastq_meta = []
    //if (!file(row.read_path).exists()) {
    //    exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    //}
    //if (meta.single_end) {
    fastq_meta = [ meta, [ file(row.reads) ] ]
    //} else {
    //    if (!file(row.fastq_2).exists()) {
    //         exit 1, "ERROR: Pl since most the trains are down ease check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
    //     }
    return fastq_meta
    //return meta
}
