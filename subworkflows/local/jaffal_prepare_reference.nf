//
// Download and prepare reference files for JAFFAL
//

//include { GUNZIP as GUNZIP_JAFFAL } from '../../../modules/nf-core/gunzip'
include { CURL  as CURL_JAFFAL  } from '../../modules/local/curl'
include { UNTAR as UNTAR        } from '../../modules/nf-core/untar/main'

workflow JAFFAL_PREPARE_REFERENCE {

    //take:

    main:

    ch_versions = Channel.empty()

    CURL_JAFFAL( 'jaffal_reference', '.gz', 'https://figshare.com/ndownloader/articles/27673314/versions/1' )
    ch_jaffal_reference_tar = CURL_JAFFAL.out.curl
    UNTAR( ch_jaffal_reference_tar )
    UNTAR.out.untar
        .map { it  -> [ it[1] ]}
        .set { ch_jaffal_reference_dir }

    ch_versions = ch_versions.mix(CURL_JAFFAL.out.versions)
    ch_versions = ch_versions.mix(UNTAR.out.versions)

    emit:
    jaffal_reference_dir = ch_jaffal_reference_dir
    versions = ch_versions                     // channel: [ versions.yml ]
}
