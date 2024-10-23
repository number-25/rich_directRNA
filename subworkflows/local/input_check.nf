/* Checks the input channel and creates channel map */

include { CHECK_SAMPLESHEET } from '../../modules/local/check_samplesheet.nf'

workflow INPUT_CHECK {

    take:
    samplesheet // file: /path/to/samplesheet

    main:

    CHECK_SAMPLESHEET( samplesheet )
        .samplesheet
        .splitCsv ( header:true, sep:',' )
        // .map { get_sample_info(it, params.genomes) }
        //nanoseq option .map { it -> [ it[0], it[1], it[2], it[3] ] }
        //rnaseq option
        .map { it -> [ it[0], it[1], it[2], it[3] ] }
        .map { create_fastq_channel(it) }
        .set { ch_sample }

    emit:
    ch_sample // [ sample, replicate, sequencing_summary_file, path_to_reads ]
    //ch_versions = ch_versions.mix(CHECK_SAMPLESHEET.out.versions.first())
}

// Create a meta map from the samplesheet
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id           = row.sample
    meta.replicate   = row.replicate
    meta.sequencing_summary = row.sequencing_summary
    meta.fastq = row.read_path

    /*
    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.read_path).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        fastq_meta = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    return fastq_meta
    */
    return meta
}
