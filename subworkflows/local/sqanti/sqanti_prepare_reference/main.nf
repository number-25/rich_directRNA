//
// Download, uncompress and prepare reference files for SQANTI QC
//

include { GUNZIP as GUNZIP_CAGE         } from '../../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_POLYA_SITES  } from '../../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_INTROPOLIS   } from '../../../../modules/nf-core/gunzip'
include { CURL as CURL_CAGE             } from '../../../../modules/local/curl/curl'
include { CURL as CURL_POLYA_SITES      } from '../../../../modules/local/curl/curl'
include { CURL as CURL_POLYA_MOTIF      } from '../../../../modules/local/curl/curl'
include { CURL as CURL_INTROPOLIS       } from '../../../../modules/local/curl/curl'


workflow SQANTI_PREPARE_REFERENCE {

    take:
    sqanti_qc_reference
    sqanti_qc_cage
    sqanti_qc_cage_path
    sqanti_qc_polyA_sites
    sqanti_qc_polyA_sites_path
    sqanti_qc_polyA_motif
    sqanti_qc_polyA_motif_path
    sqanti_qc_intron_junctions
    sqanti_qc_intron_path

    main:

    ch_versions = Channel.empty()

    if (sqanti_qc_reference == 'human') {
        // cage data
        if (sqanti_qc_cage) {
            if (sqanti_qc_cage_path == null) { // user doesn't provide path to predownload cage data
                CURL_CAGE( 'refTSS_CAGE', 'bed.gz', 'https://figshare.com/ndownloader/files/52801133' )
                sqanti_cage_bed_gzip = CURL_CAGE.out.curl
                ch_versions = ch_versions.mix(CURL_CAGE.out.versions)
                ch_sqanti_cage_bed = GUNZIP_CAGE( [ [:], sqanti_cage_bed_gzip ] ).gunzip.map { it[1] }
                ch_versions = ch_versions.mix(GUNZIP_CAGE.out.versions)
                //ch_sqanti_cage_bed = GUNZIP_CAGE.out.
            } else {
                if (sqanti_qc_cage_path.endsWith('.gz')) {
                ch_sqanti_cage_bed = GUNZIP_CAGE( [ [:], sqanti_qc_cage_path ] ).gunzip.map { it[1] }
                ch_versions = ch_versions.mix(GUNZIP_CAGE.out.versions)
                } else {
                ch_sqanti_cage_bed = Channe.value(file(sqanti_qc_cage_path), checkIfExists:true)
                }
            }
        }
        // poly_A sites
        if (sqanti_qc_polyA_sites) {
            if (sqanti_qc_polyA_sites_path == null) { // user doesn't provide path to predownload polyA sites data
                CURL_POLYA_SITES( 'polyA_sites', 'bed.gz', 'https://figshare.com/ndownloader/files/52801130' )
                sqanti_polyA_sites_bed_gzip = CURL_POLYA_SITES.out.curl
                ch_versions = ch_versions.mix(CURL_POLYA_SITES.out.versions)
                ch_sqanti_polyA_sites_bed = GUNZIP_POLYA_SITES( [ [:], sqanti_polyA_sites_bed_gzip ] ).gunzip.map { it[1] }
                ch_versions = ch_versions.mix(GUNZIP_POLYA_SITES.out.versions)
            } else {
                if (sqanti_polyA_sites_path.endsWith('.gz')) {
                ch_sqanti_qc_polyA_sites_bed = GUNZIP_POLYA_SITES( [ [:], sqanti_qc_polyA_sites_path ] ).gunzip.map { it[1] }
                ch_versions = ch_versions.mix(GUNZIP_POLYA_SITES.out.versions)
                } else {
                ch_sqanti_qc_polyA_sites_bed = Channel.value(file(sqanti_qc_polyA_sites_path), checkIfExists:true)
                }
            }
        }

        // poly_A motifs
        if (sqanti_qc_polyA_motif) {
            if (sqanti_qc_polyA_motif_path == null) { // user doesn't provide path to predownload polyA sites data
                CURL_POLYA_MOTIF( 'polyA_motif', 'txt', 'https://figshare.com/ndownloader/files/52801139' )
                ch_versions = ch_versions.mix(CURL_POLYA_MOTIF.out.versions)
                ch_sqanti_polyA_motif = CURL_POLYA_MOTIF.out.curl
            } else {
                ch_sqanti_polyA_motif = Channel.value(file(sqanti_qc_polyA_motif_path), checkIfExists:true)
            }
        }

        // intropolis motifs
        if (sqanti_qc_intron_junctions) {
            if (sqanti_qc_intron_path == null) { // user doesn't provide path to predownload polyA sites data
                CURL_INTROPOLIS( 'intropolis', 'bed.gz', 'https://figshare.com/ndownloader/files/52801127' )
                sqanti_intron_junctions_bed_gzip = CURL_INTROPOLIS.out.curl
                ch_versions = ch_versions.mix(CURL_INTROPOLIS.out.versions)
                ch_sqanti_intron_junctions_bed = GUNZIP_INTROPOLIS( [ [:], sqanti_intron_junctions_bed_gzip ] ).gunzip.map { it[1] }
            } else {
                if (sqanti_qc_intron_path.endsWith('.gz')) {
                ch_sqanti_qc_intron_junctions_bed = GUNZIP_INTROPOLIS( [ [:], sqanti_qc_intron_path ] ).gunzip.map { it[1] }
                ch_versions = ch_versions.mix(GUNZIP_INTROPOLIS.out.versions)
                } else {
                ch_sqanti_qc_intron_bed = Channel.value(file(sqanti_qc_intron_path), checkIfExists:true)
                }
            }
        }
    }


    emit:
    sqanti_qc_cage_bed = ch_sqanti_cage_bed
    sqanti_qc_polyA_sites_bed = ch_sqanti_polyA_sites_bed
    sqanti_qc_polyA_motif = ch_sqanti_qc_polyA_motif
    sqanti_qc_intron_junctions_bed = ch_sqanti_qc_intron_junctions_bed
    versions = ch_versions                     // channel: [ versions.yml ]
}
