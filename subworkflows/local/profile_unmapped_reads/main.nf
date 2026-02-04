//include { SYLPH_PREPARE_REFERENCE   } from '../../../modules/local/sylph/prepare_reference'
include { SYLPH_PROFILE             } from '../../../modules/local/sylph/profile'
include { GET_SYLPH_TAX_DB          } from '../../../modules/local/sylph/get_tax_database'
include { SYLPH_TAX                 } from '../../../modules/local/sylph/tax'

workflow PROFILE_UNMAPPED_READS {

    take:
    sylph_database
    unmapped_reads
    sylph_database_name

    main:

    ch_versions = Channel.empty()

    SYLPH_PROFILE( unmapped_reads, sylph_database )
    ch_sylph_profile = SYLPH_PROFILE.out.sylph_profile

    GET_SYLPH_TAX_DB()
    ch_sylph_tax_db = GET_SYLPH_TAX_DB.out.sylph_tax_db

    SYLPH_TAX( ch_sylph_profile, ch_sylph_tax_db, sylph_database_name )
    ch_sylph_tax = SYLPH_TAX.out.sylph_tax

    emit:

    sylph_tax = ch_sylph_tax

    versions = ch_versions                     // channel: [ versions.yml ]
}
