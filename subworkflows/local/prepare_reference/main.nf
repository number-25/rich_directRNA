// SUBWORKFLOW: PREPARE_REFERENCE
// Uncompress and prepare reference files
//

include { GUNZIP as GUNZIP_FASTA                    } from '../../../modules/nf-core/gunzip'
include { CUSTOM_GETCHROMSIZES                      } from '../../../modules/nf-core/custom/getchromsizes'
include { MINIMAP2_INDEX as MINIMAP2_GENOME_INDEX   } from '../../../modules/local/minimap2/index'
include { MINIMAP2_INDEX as MINIMAP2_TXOME_INDEX    } from '../../../modules/local/minimap2/index'
include { GUNZIP as GUNZIP_TRANSCRIPTOME            } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_ANNOTATION_GTF           } from '../../../modules/nf-core/gunzip'
include { SQANTI_PREPARE_REFERENCE                  } from '../../local/sqanti/sqanti_prepare_reference'
include { JAFFAL_PREPARE_REFERENCE                  } from '../../local/jaffal_prepare_reference'
include { UNZIP                                     } from '../../../modules/local/unzip/unzip'

// prepare additional files
//TO-DO make these modules
//include { GXF2BED as GTF_TO_BED } from '../../../modules/local/gxf2bed' // gxf2bed module
//include { BIGWIG_TO_WIG } from '../../../modules/local/bigwigtowig'
//include { BEDOPS as WIG_TO_BED } from '../../../modules/local/bedops'
//include { SAMTOOLS_SORT } from '../../../modules/nf-core/samtools/sort/main'

workflow PREPARE_REFERENCE{

    take:
    genome_fasta                    // file: /path/to/genome_fasta.fa
    genome_fasta_index              // file: /path/to/genome_fasta_index.fa.fai
    genome_fasta_sizes              // file: /path/to/genome_fasta.sizes
    genome_minimap2_index           // file: /path/to/minimap2_genome_index.fa.mmi
    bam_input                       // boolean: false [default: false]
    transcriptome_fasta             // file: /path/to/genome_fasta.sizes
    transcriptome_minimap2_index    // file: /path/to/minimap2_transcriptome_index.fa.mmi
    annotation_gtf                  // file: /path/to/annotation.gtf
    //appris_bed?
    //mane_select_bed?
    //mane_clinical_bed?
    skip_jaffal                     // boolean: skip jaffal fusion gene detection [default: false]
    skip_jaffal_download            // boolean: skip jaffal fusion gene detection [default: false]
    jaffal_reference                // file: /path/to/jaffal_reference
    skip_transcript_quantification  // boolean: skip all transcript quantification [default: false]
    skip_sqanti_all                 // boolean: skip all of sqanti [default: false]
    skip_sqanti_qc                  // boolean: skip sqanti qc [default: false]
    sqanti_qc_reference             // boolean: three values options [mouse, human, custom]
    sqanti_qc_cage                  // boolean: true [default: true]
    sqanti_qc_cage_path             // file path: [default: in config]
    sqanti_qc_polyA_sites           // boolean: true [default: true]
    sqanti_qc_polyA_sites_path      // file path: [default: in config]
    sqanti_qc_polyA_motif           // boolean: true [default: true]
    sqanti_qc_polyA_motif_path      // file path: [default: in config]
    sqanti_qc_intron_junctions      // boolean: true [default: true]
    sqanti_qc_intron_path           // file path: [default: in config]

    main:

    ch_versions = Channel.empty()

    // Uncompress genome fasta file
    // Mandatory input
    if (genome_fasta) {
        file(genome_fasta, checkIfExists: true)
        if (genome_fasta.endsWith('.gz')) {
            ch_genome_fasta = GUNZIP_FASTA( [ [:], genome_fasta ] ).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
        } else {
    //which one below?
            //ch_genome_fasta = Channel.value(file(genome_fasta), checkIfExists: true)
            ch_genome_fasta = Channel.fromPath(params.genome_fasta, checkIfExists: true)
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
            ch_transcriptome_fasta = GUNZIP_TRANSCRIPTOME( [ [:], transcriptome_fasta ] ).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_TRANSCRIPTOME.out.versions)
} else {
            //ch_transcriptome_fasta = Channel.value(file(transcriptome_fasta), checkIfExists: true)
            ch_transcriptome_fasta = Channel.fromPath(params.transcriptome_fasta, checkIfExists: true)
        }
    }

    // Uncompress GTF annotation file
    // Mandatory input
    if (annotation_gtf) {
        file(annotation_gtf, checkIfExists:true)
        if (annotation_gtf.endsWith('.gz')) {
            ch_annotation_gtf = GUNZIP_ANNOTATION_GTF( [ [:], annotation_gtf ] ).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_ANNOTATION_GTF.out.versions)
        } else {
//            ch_annotation_gtf = Channel.value(file(annotation_gtf), checkIfExists: true)
            ch_annotation_gtf = Channel.fromPath(params.annotation_gtf, checkIfExists: true)
        }
    }

    // Initialise genome minimap2 index if provided
    // If bam input is provided, skip minimap2 genome indexing
    if (!bam_input) {
        //if (genome_minimap2_index == null) {
        if (!genome_minimap2_index) {
            MINIMAP2_GENOME_INDEX( ch_genome_fasta )
            ch_genome_minimap2_index = MINIMAP2_GENOME_INDEX.out.index
            ch_versions = ch_versions.mix(MINIMAP2_GENOME_INDEX.out.versions)
        } else {
            //ch_genome_minimap2_index = Channel.value(file(genome_minimap2_index), checkIfExists: true)
            ch_genome_minimap2_index = Channel.fromPath(params.genome_minimap2_index, checkIfExists: true)
        }
    } else {
        ch_genome_minimap2_index = null
    }

    // Initialise transcriptome minimap2 index if provided
    // If bam input is provided, skip minimap2 transcriptome indexing
    if (!skip_transcript_quantification) {
        //if (genome_minimap2_index == null) {
        if (!transcriptome_minimap2_index) {
            MINIMAP2_TXOME_INDEX( ch_transcriptome_fasta )
            ch_transcriptome_minimap2_index = MINIMAP2_TXOME_INDEX.out.index
            ch_versions = ch_versions.mix(MINIMAP2_TXOME_INDEX.out.versions)
        } else {
            //ch_genome_minimap2_index = Channel.value(file(genome_minimap2_index), checkIfExists: true)
            ch_transcriptome_minimap2_index = Channel.fromPath(params.transcriptome_minimap2_index, checkIfExists: true)
        }
    } else {
        ch_transcriptome_minimap2_index = null
    }

    // Prepare references for SQANTI QC
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
        ch_sqanti_qc_cage_bed               = SQANTI_PREPARE_REFERENCE.out.sqanti_qc_cage_bed
        ch_sqanti_qc_polyA_sites_bed        = SQANTI_PREPARE_REFERENCE.out.sqanti_qc_polyA_sites_bed
        ch_sqanti_qc_polyA_motif            = SQANTI_PREPARE_REFERENCE.out.sqanti_qc_polyA_motif
        ch_sqanti_qc_intron_junctions_bed   = SQANTI_PREPARE_REFERENCE.out.sqanti_qc_intron_junctions_bed

    } else {
        ch_sqanti_qc_cage_bed               = null
        ch_sqanti_qc_polyA_sites_bed        = null
        ch_sqanti_qc_polyA_motif            = null
        ch_sqanti_qc_intron_junctions_bed   = null
    }

    // Prepare reference for JAFFAL
    if (!skip_jaffal) {
        if (!skip_jaffal_download) {
            JAFFAL_PREPARE_REFERENCE()
            ch_jaffal_reference_dir = JAFFAL_PREPARE_REFERENCE.out.jaffal_reference_dir
        } else if (params.jaffal_reference.endsWith('.zip')) {
        //} else {
        // reference directory needs to be gzipped for this function to work
            ch_jaffal_reference = file(params.jaffal_reference, checkIfExists: true)
            UNZIP( ch_jaffal_ref_dir, "jaffal_reference" )
            ch_jaffal_reference_dir = UNZIP.out.unzipped_archive
            ch_versions = ch_versions.mix(UNZIP.out.versions)
        } else {
            ch_jaffal_reference_dir = file(params.jaffal_reference, checkIfExists: true)
        }
    } else {
        ch_jaffal_reference_dir = null
    }


    emit:
    genome_fasta                    = ch_genome_fasta
    genome_fasta_index              = ch_genome_fasta_index
    genome_fasta_sizes              = ch_genome_fasta_sizes
    genome_minimap2_index           = ch_genome_minimap2_index
    transcriptome_fasta             = ch_transcriptome_fasta
    transcriptome_minimap2_index    = ch_transcriptome_minimap2_index
    annotation_gtf                  = ch_annotation_gtf
    sqanti_qc_cage_bed              = ch_sqanti_qc_cage_bed
    sqanti_qc_polyA_sites_bed       = ch_sqanti_qc_polyA_sites_bed
    sqanti_qc_polyA_motif           = ch_sqanti_qc_polyA_motif
    sqanti_qc_intron_junctions_bed  = ch_sqanti_qc_intron_junctions_bed
    jaffal_reference                = ch_jaffal_reference_dir
    //phylop_bed = ch_phylop_bed
    versions                        = ch_versions // channel: [ versions.yml ]
}
