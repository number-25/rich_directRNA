//
// Download and prepare reference files for JAFFAL
//

//include { GUNZIP as GUNZIP_JAFFAL } from '../../../modules/nf-core/gunzip'
include { CURL as CURL_JAFFAL } from '../../modules/local/curl'

workflow JAFFAL_PREPARE_REFERENCE {

    //take:

    main:

    ch_versions = Channel.empty()

    CURL_JAFFAL( 'jaffal_reference', '.gz', 'https://figshare.com/ndownloader/articles/27673314/versions/1' )
    ch_jaffal_reference = CURL_JAFFAL.out.curl
    ch_versions = ch_versions.mix(CURL_JAFFAL.out.versions)

    emit:
    jaffal_reference = ch_jaffal_reference
    versions = ch_versions                     // channel: [ versions.yml ]
}
