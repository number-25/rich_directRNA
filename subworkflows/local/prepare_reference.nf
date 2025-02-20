//
// Uncompress and prepare reference genome files
//

include { GUNZIP as GUNZIP_FASTA } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GTF } from '../../../modules/nf-core/gunzip'
//include { GUNZIP as GUNZIP_BED } from '../modules/nf-core/gunzip'
//include { GUNZIP as GUNZIP_TRANSCRIPT_FASTA } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_CAGE } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_POLYA } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_INTROPOLIS } from '../../../modules/nf-core/gunzip'

// prepare indices for reference
include { CUSTOM_GETCHROMSIZES } from '../../../modules/nf-core/custom/getchromsizes'

// prepare additional files
//TO-DO make these modules
//include { GXF2BED as GTF_TO_BED } from '../../../modules/local/gxf2bed' // gxf2bed module
//include { BIGWIG_TO_WIG } from '../../../modules/local/bigwigtowig'
//include { BEDOPS as WIG_TO_BED } from '../../../modules/local/bedops'
//include { SAMTOOLS_SORT } from '../../../modules/nf-core/samtools/sort/main'

workflow PREPARE_REFERENCE {

    take:
    genome_fasta                // file: /path/to/genome_fasta.fa
    genome_fasta_index          // file: /path/to/genome_fasta_index.fa.fai
    genome_fasta_minimap2_index // file: /path/to/minimap2_genome_index.fa.mai
    genome_fasta_sizes          // file: /path/to/genome_fasta.sizes
    annotation_gtf              // file: /path/to/annotation.gtf
    cage_bed                    // file: /path/to/cage.bed
    polyA_bed                   // file: /path/to/polyA.bed
    polyA_sites                 // file: /path/to/polyA_sites.txt
    intropolis_bed              // file: /path/to/intropolis.bed
    //appris_bed?
    //mane_select_bed?
    //mane_clinical_bed?
    // indices
    skip_sqanti_qc             // boolean: skip all of sqanti

    main:

    ch_versions = Channel.empty()

    // Uncompress genome fasta file
    // Mandatory input
    if (genome_fasta) {
        if (genome_fasta.endsWith('.gz')) {
            ch_genome_fasta = GUNZIP_GTF( [ [:], genome_fasta ] ).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
        } else {
            //which one below?
            ch_genome_fasta = Channel.value(file(genome_fasta), checkIfExists: true)
            //ch_annotation_gtf = Channel.fromPath(params.annotation_gtf, checkIfExists: true)
        }
    }

    // Uncompress GTF annotation file
    // Mandatory input
    if (annotation_gtf) {
        if (annotation_gtf.endsWith('.gz')) {
            ch_annotation_gtf = GUNZIP_GTF( [ [:], annotation_gtf ] ).gunzip.map { it[1] }
        } else {
            //which one below?
            ch_annotation_gtf = Channel.value(file(annotation_gtf), checkIfExists: true)
            //ch_annotation_gtf = Channel.fromPath(params.annotation_gtf, checkIfExists: true)
        }
    }


    if (params.sqanti_download_reference) {
        // download references - perhaps upload them to figshare for easiest
        // downloads?
    }


    // Seed orthogonal files for SQANTI3
    if (!params.skip_sqanti_qc) {
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
    //genome_fasta = ch_genome_fasta
    //transcriptome_fasta = ch_transcriptome_fasta
    annotation_gtf = ch_annotation_gtf
    //annotation_bed =
    cage_bed = ch_cage_bed
    polyA_bed = ch_polyA_bed
    polyA_sites = ch_polyA_sites
    intropolis_bed = ch_intropolis_bed
    //phylop_bed = ch_phylop_bed
    versions = ch_versions                     // channel: [ versions.yml ]
}
