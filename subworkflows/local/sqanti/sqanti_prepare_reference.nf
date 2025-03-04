//
// Download, uncompress and prepare reference files for SQANTI QC
//

include { GUNZIP as GUNZIP_CAGE } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_POLYA_SITES } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_INTROPOLIS } from '../../../modules/nf-core/gunzip'

workflow PREPARE_REFERENCE {

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

    // Uncompress genome fasta file
    // Mandatory input

    if (sqanti_qc_reference == 'human') {
        // cage data
        if (sqanti_qc_cage)
            if (sqanti_qc_cage_path == null) // user doesn't provide path to predownload cage data
               CURL_CAGE('refTSS_CAGE', 'bed', 'https://reftss.riken.jp/datafiles/current/human/refTSS_v4.1_human_coordinate.hg38.bed.txt.gz')
               ch_sqanti_cage


    //} else { // mouse data - down the line
    }

    if (genome_fasta) {
        file(genome_fasta, checkIfExists: true)
        if (genome_fasta.endsWith('.gz')) {
            ch_genome_fasta = GUNZIP_GTF( [ [:], genome_fasta ] ).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
        } else {
     //which one below?
            ch_genome_fasta = Channel.value(file(genome_fasta), checkIfExists: true)
            //ch_annotation_gtf = Channel.fromPath(params.annotation_gtf, checkIfExists: true)
        }
    }

    // Genome fasta index

    if (!genome_fasta_index) {
        CUSTOM_GETCHROMSIZES( ch_genome_fasta )
        ch_genome_fasta_index = CUSTOM_GETCHROMSIZES.out.fai
        ch_versions = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions.first())
    } else {
        ch_genome_fasta_index = Channel.value(file(genome_fasta_index, checkIfExists: true))
    }

    // Genome fasta sizes

    if (!genome_fasta_sizes) {
        CUSTOM_GETCHROMSIZES( ch_genome_fasta )
        ch_genome_fasta_sizes = CUSTOM_GETCHROMSIZES.out.sizes
        ch_versions = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions.first())
    } else {
        ch_genome_fasta_sizes = Channel.value(file(genome_fasta_sizes, checkIfExists: true))
    }

    // Uncompress transcriptome fasta file
    // Mandatory input
    if (transcriptome_fasta) {
        file(transcriptome_fasta, checkIfExists: true)
        if (transcriptome_fasta.endsWith('.gz')) {
            ch_transcriptome_fasta = GUNZIP_GTF( [ [:], transcriptome_fasta ] ).gunzip.map { it[1] }
        } else {
            //which one below?
            ch_transcriptome_fasta = Channel.value(file(genome_fasta), checkIfExists: true)
            //ch_annotation_gtf = Channel.fromPath(params.annotation_gtf, checkIfExists: true)
        }
    }
    // Uncompress GTF annotation file
    // Mandatory input
    if (annotation_gtf) {
        file(annotation_gtf, checkIfExists:true)
        if (annotation_gtf.endsWith('.gz')) {
            ch_annotation_gtf = GUNZIP_GTF( [ [:], annotation_gtf ] ).gunzip.map { it[1] }
        } else {
            //which one below?
            ch_annotation_gtf = Channel.value(file(annotation_gtf), checkIfExists: true)
            //ch_annotation_gtf = Channel.fromPath(params.annotation_gtf, checkIfExists: true)
        }
    }

    // Initialise minimap2 index if provided
    // If bam input is provided, skip minimap2 genome indexing
    if (!bam_input) {
        if (genome_fasta_minimap2_index == null) {
            MINIMAP2_INDEX( ch_genome_fasta )
            ch_genome_minimap2_index = MINIMAP2_INDEX.out.bai
            ch_versions = ch.versions.mix(MINIAP2_INDEX.out.versions)
        } else {
            ch_genome_minimap2_index = Channel.value(file(genome_fasta_minimap2_index), checkIfExists: true)
        }
    }

    // SQANTI

    if (!skip_sqanti_all || !skip_sqanti_qc) {
        SQANTI_PREPARE_REFERENCE(
            sqanti_qc_reference,  // human, mouse or custom
            sqanti_qc_cage,       // boolean
            sqanti_qc_cage_path,
            sqanti_qc_polyA_sites,
            sqanti_qc_polyA_sites_path,
            sqanti_qc_polyA_motif,
            sqanti_qc_polyA_motif_path,
            sqanti_qc_intron_junctions,
            sqanti_qc_intron_path
            )
    }

    // Seed orthogonal files for SQANTI3
    if (!params.skip_sqanti_all || !skip_sqanti_qc) {
    // Uncompress cage BED interval file
    // Optional input
    // or is the pattern if (params.cage_bed)
        if (cage_bed) {
            if (cage_bed.endsWith('.gz')) {
                ch_cage_bed = GUNZIP_CAGE( [ [:], cage_bed ] ).gunzip.map { it[1] }
            } else {
                ch_cage_bed = Channel.value(file(cage_bed), checkIfExists: true)
            }
        }
    //
    // Uncompress polyA BED interval file
    //
        if (polyA_bed) {
            if (polyA_bed.endsWith('.gz')) {
                ch_polyA_bed = GUNZIP_POLYA( [ [:], polyA_bed ] ).gunzip.map { it[1] }
            } else {
                ch_polyA_bed = Channel.value(file(polyA_bed), checkIfExists: true)
            }
        }
    //
    // Initialize polyA sites file
    //
        if (polyA_sites) {
            } else {
                ch_polyA_sites = Channel.value(file(polyA_sites), checkIfExists: true)
                //ch_polyA_sites = Channel.fromPath(polyA_sites, checkIfExists = true)
            }
        }
    //
    // Uncompress intropolis BED interval file
    //
        if (intropolis_bed) {
            if (intropolis_bed.endsWith('.gz')) {
                ch_intropolis_bed = GUNZIP_INTROPOLIS ( [ [:], intropolis_bed ] ).gunzip.map { it[1] }
            } else {
                ch_intropolis_bed = Channel.value(file(intropolis_bed), checkIfExists: true)
            }
        }
    }

    //
    // Convert PhyloP bigWig to bed
    // Hold off on this - extremely memory intensive process
    //if (phylop_bigwig) {
    //    ch_phylop_wig = BIGWIG_TO_WIG( ch_phylop_bigwig ).out.phylop_wig
    //    ch_phylop_bed = WIG_TO_BED( ch_phylop_wig ).out.phylop_bed
    //    ch_versions = ch_versions.mix(WIG_TO_BED.out.versions)
    //}

//-------------------------------------------------------//
//              Prepare indices                         //

    emit:
    genome_fasta = ch_genome_fasta

    sqanti_cage_bed = ch_sqanti_cage_bed
    polyA_bed = ch_polyA_bed
    polyA_sites = ch_polyA_sites
    intropolis_bed = ch_intropolis_bed
    //phylop_bed = ch_phylop_bed
    versions = ch_versions                     // channel: [ versions.yml ]
}
