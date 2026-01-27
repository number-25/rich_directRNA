//include { SYLPH_PREPARE_REFERENCE   } from '../../../modules/local/sylph/prepare_reference'
include { SYLPH_PROFILE             } from '../../../modules/local/sylph/profile'
include { SYLPH_TAX                 } from '../../../modules/local/sylph/tax'

workflow PROFILE_UNMAPPED_READS {

    take:
    sylph_database
    unmapped_reads
    sylph_database_name

    main:

    ch_versions = Channel.empty()

    //ch_sylph_database = SYLPH_PREPARE_REFERENCE.out.database

    SYLPH_PROFILE( unmapped_reads, sylph_database )
    ch_sylph_profile = SYLPH_PROFILE.out.sylph_profile

    SYLPH_TAX( ch_sylph_profile, sylph_database_name )
    ch_sylph_tax = SYLPH_TAX.out.sylph_tax

    emit:

    sylph_tax = ch_sylph_tax

    versions = ch_versions                     // channel: [ versions.yml ]
}
