/*
----------------------------------------------------------------------------------------
    VALIDATE INPUTS
----------------------------------------------------------------------------------------
*/

// nextflow magik

//def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

//

// Check mandatory parameters (missing protocol or profile will exit the run.)
if (params.input) {
    ch_input = file(params.input) // defined in nextflow.config
} else {
    exit 1, 'Input samplesheet not specified!'
}

// Function to check if running offline
def isOffline() {
    try {
        return NXF_OFFLINE as Boolean
    }
    catch( Exception e ) {
        return false
    }
}

/*
----------------------------------------------------------------------------------------
    IMPORT LOCAL MODULES / SUBWORKFLOWS / FUNCTIONS
----------------------------------------------------------------------------------------
*/
// Check samplesheet
include { INPUT_CHECK               } from '../subworkflows/local/input_check'
// Load Unix Utils
include { GUNZIP as GUNZIP_FASTA    } from '../modules/nfcore/gunzip'
// fastq QC
include { NANOQ                     } from '../modules/local/nanoq'
// fastq mapping
include { MAPPING                   } from '../subworkflows/local/mapping'
include { SAMTOOLS_FAIDX            } from '../modules/local/samtools/samtools_faidx'
// bam QC
include { BAM_QC                    } from '../subworkflows/local/bam_qc'
// transcript reconstruction
include { BAM_TO_BED12              } from '../modules/local/flair/bam_to_bed12'
include { FLAIR_CORRECT             } from '../modules/local/flair/flair_correct'
include { FLAIR_COLLAPSE            } from '../modules/local/flair/flair_collapse'
include { BED_TO_BAM                } from '../modules/local/bedtools/bed_to_bam'
include { BAMBU                     } from '../modules/local/bambu'
//include { ISOQUANT                } from '../modules/local/isoquant'
//include { ISOQUANT_CORRECTION     } from '../modules/local/isoquant_correct'
// fusion gene detection
//include { JAFFAL             } from '../modules/local/jaffal'
// transcriptome assessment
//include { SQANTI               } from '../subworkflows/local/sqanti'
//include { SQANTI_QC            } from '../modules/local/sqanti/sqanti_qc'
//include { SQANTI_FILTER        } from '../modules/local/sqanti/sqanti_filter'
//include { SQANTI_RESCUE        } from '../modules/local/sqanti/sqanti_rescue'



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//include { FASTQC                 } from '../modules/nf-core/fastqc/main'
//include { MULTIQC                } from '../modules/nf-core/multiqc/main'
//include { paramsSummaryMap       } from 'plugin/nf-validation'
//include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
//include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_directrna_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DIRECTRNA{

    //take:
    //ch_samplesheet // channel: samplesheet read in from --input
    //main:

    ch_versions = Channel.empty()
    //ch_multiqc_files = Channel.empty()

    // INPUT_CHECK
    INPUT_CHECK ( ch_input )
        .set { ch_sample }

    // playing around with channel transformations
    //ch_sample
        //.map { it -> [ it[0], it[1] ] } // take sample, replicate, reads
        //.set { ch_fastq }

    //
    // QC of fastq files
    // Toulligqc?
    // MODULE: NANOQ
    if (!params.skip_qc) {
        NANOQ ( ch_sample )
        ch_versions = ch_versions.mix(NANOQ.out.versions.first())
    }

    //
    // MODULE: Run FastQC
    //
    //FASTQC (
    //    ch_samplesheet
    //)
    //ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    //ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // Prepare the reference files
    // SUBWORKFLOW: PREPARE_REFERENCE
    //
    if (!params.skip_prepare_reference) {
        PREPARE_REFERENCE ()
        //ch_minimap2_genome_index = PREPARE_REFERENCE.minimap2_index
        ch_transcriptome_fasta = PREPARE_REFERENCE.out.transcriptome_fasta
        ch_versions = ch_versions.mix(PREPARE_REFERENCE.out.versions.first())
        ch_cage_bed = PREPARE_REFERENCE.out.cage_bed
        ch_polyA_bed = PREPARE_REFERENCE.out.polyA_bed
        ch_intropolis_bed = PREPARE_REFERENCE.out.intropolis_bed
        ch_custom_chrom_sizes = PREPARE_REFERENCE.out.custom_chrom_sizes
        //ch_jaffal_ref = Channel.fromPath(jaffal_ref, checkIfExists = true)
    }

    // Mapping and sorting
    // SUBWORKFLOW: MAPPING
    // if the reference preparation is skipped, minimap2 will generate an index
    // for the provided reference file
    //
    if (!params.skip_mapping) {
    // the genome indexes are not created in prepare reference anymore
        if (params.genome_fasta) {
            ch_genome_fasta = Channel.fromPath(params.genome_fasta, checkIfExists: true)
            if (ch_genome_fasta.endWith('.gz')) {
                GUNZIP_FASTA( ch_genome_fasta )
                ch_genome_fasta = GUNZIP_FASTA.out.gunzip
                ch_versions = ch_versions.mix(GUNZIP.out.versions.first())
            }
            if (params.custom_genome) {
                // If custom genome is provided
                CUSTOM_GETCHROMSIZES( ch_genome_fasta )
                ch_genome_samtools_idx = CUSTOM_GETCHROMSIZES.out.fai
                ch_genome_sizes = CUSTOM_GETCHROMSIZES.out.sizes
                ch_versions = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions.first())
                //SAMTOOLS_FAIDX(ch_genome_fasta)
                //ch_genome_samtools_idx = SAMTOOLS_FAIDX.out.index
                MINIMAP2_INDEX( ch_genome_fasta )
                ch_genome_minimap2_idx = MINIMAP2_INDEX.out.index
            } else {
                ch_genome_samtools_idx = Channel.fromPath(params.genome_fasta_samtools_index, checkIfExists: true)
                ch_genome_minimap2_idx = Channel.fromPath(params.genome_fasta_minimap2_index checkIfExists: true)
                ch_genome_sizes = Channel.fromPath(params.genome_fasta_sizes checkIfExists: true)
            }
            MAPPING( ch_sample, ch_genome_minimap2_idx )
            ch_bam = MAPPING.out.bam
            ch_bam_index = MAPPING.out.bai
            ch_versions = ch_versions.mix(MAPPING.out.versions.first())
            //ch_mixed_bam = ch_bam.mix(ch_bam_indx)
            //MINIMAP2_ALIGN( ch_sample, ch_genome_fasta )
            } else {
                exit 1, 'Reference genome fasta file is not specified! please modify nextflow.config or use --genome_fasta parameter'
                    }
        } else {
            MAPPING( ch_sample, ch_genome_minimap2_idx)
            ch_bam = MAPPING.out.bam
            ch_bam_index = MAPPING.out.bai
            ch_versions = ch_versions.mix(MAPPING.out.versions.first())
            //ch_mixed_bam = ch_bam.mix(ch_bam_indx)
        }

    //
    // BAM QC
    // SUBWORKFLOW: BAM_QC
    // Execute cramino, alfred and samtools flagstat on bam output from mapping
    if (!params.skip_mapping) {
        BAM_QC( ch_bam, ch_genome_fasta )
        ch_versions = ch_versions.mix(BAM_QC.out.versions.first())
    }
    // If a raw BAM is provided and mapping is not needed
    // Will need to index it
    else {
        ch_bam = Channel.fromPath(params.bam_input, checkIfExists: true)
        BAM_QC( ch_bam, ch_genome_fasta )
        ch_versions = ch_versions.mix(BAM_QC.out.versions.first())
        }

    //
    // Transcript Reconstruction
    //
    // Source transcriptome annotation.gtf

    //ch_annotation_gtf = file(params.annotation_gtf) // check if exists
    ch_annotation_gtf = Channel.fromPath(params.annotation_gtf, checkIfExists: true) // check if exists

    if (!params.skip_flair_correct && !params.skip_flair_collapse) {
        BAM_TO_BED12( ch_bam, ch_bam_index )
        // seeing if a mixed channel with bam and bam.bai works
        //BAM_TO_BED12( ch_mixed_bam )
        ch_mapped_bed = BAM_TO_BED12.out.bed
        FLAIR_CORRECT( ch_mapped_bed, ch_genome_fasta, ch_annotation_gtf )
        ch_corrected_bed = FLAIR_CORRECT.out.flair_corrected_bed
        FLAIR_COLLAPSE( ch_corrected_bed, ch_sample, ch_annotation_gtf, ch_genome_fasta )
        ch_collapsed_bed = FLAIR_COLLAPSE.out.collapsed_isoforms_bed
        ch_collapsed_gtf = FLAIR_COLLAPSE.out_collapsed_isoforms.gtf
        ch_versions = ch_versions.mix(FLAIR_collapse.out.versions.first())
        //ch_collapsed_bed
        //   .map { it -> [ it[0], it[1] ] }
        //   .set { ch_test_bed }
        //BED_TO_BAM( ch_collapsed_bed, ch_genome_fasta_sizes )
        //ch_collapsed_bam = BED_TO_BAM.out.collapsed_bed
    } else {
        // No collapsing just correction
        if (!params.skip_flair_correct && params.skip_flair_collapse) {
            BAM_TO_BED12( ch_bam, ch_bam_index )
            //BAM_TO_BED12( ch__mixed_bam )
            ch_mapped_bed = BAM_TO_BED12.out.bed
            FLAIR_CORRECT( ch_mapped_bed, ch_genome_fasta, ch_annotation_gtf )
            ch_corrected_bed = FLAIR_CORRECT.out.flair_corrected_bed
            ch_versions = ch_versions.mix(FLAIR_CORRECT.out.versions.first())
        } else {
            BAM_TO_BED12( ch_bam, ch_bam_index )
            //BAM_TO_BED12( ch_mixed_bam )
            ch_mapped_bed = BAM_TO_BED12.out.bed
            FLAIR_COLLAPSE( ch_mapped_bed, ch_sample, ch_annotation_gtf, ch_genome_fasta )
            ch_collapsed_bed = FLAIR_COLLAPSE.out.collapsed_isoforms_bed
            ch_collapsed_gtf = FLAIR_COLLAPSE.out.collapsed_isoforms_gtf
            ch_versions = ch_versions.mix(FLAIR_COLLAPSE.out.versions.first())
            //ch_collapsed_bed
            //.map { it -> [ it[0], it[1] ] }
            //.set { ch_test_bed }
            //BED_TO_BAM( ch_collapsed_bed, ch_genome_fasta_sizes )
            //ch_collapsed_bam = BED_TO_BAM.out.collapsed_bam
        }
    }


    // Read correction tools? Which ones....
    // TC-CLEAN?
    // IsoQUANT?
    // FLAIR

    // BAMBU
    if (!params.skip_bambu) {
        BAMBU( ch_genome_fasta, ch_annotation_gtf, ch_bam )
        ch_bambu_gtf = BAMBU.out.bambu_extended_gtf
        ch_versions = ch_versions.mix(BAMBU.out.versions.first())
        // MIX genome fasta with fasta index as this will improve GFFREADs speed
        GFFREAD_GETFASTA( ch_bambu_gtf, ch_genome_fasta )
        ch_bambu_transcripts = GFFREAD_GETFASTA.out.transcripts_fa
        ch_versions = ch_versions.mix(GFFREAD_GETFASTA.out.versions.first())
        }

    // ISOQUANT
    if (params.isoquant_reconstruction && params.skip_isoquant_correction) {
        ISOQUANT( ch_bam, ch_genome_fasta, ch_annotation_gtf )
        ch_isoquant_gtf = ISOQUANT.out.isoquant_transcript_gtf
        ch_versions = ch_versions.mix(ISOQUANT.out.versions.first())
        GFFREAD_GETFASTA( ch_isoquant_gtf, ch_genome_fasta )
        ch_isoquant_transcripts = GFFREAD_GETFASTA.out.transcripts_fa
        ch_versions = ch_versions.mix(GFFREAD_GETFASTA.out.versions.first())
    } else {
        ISOQUANT_CORRECTION ( ch_bam, ch_genome_fasta, ch_annotation_gtf )
        ch_isoquant_gtf = ISOQUANT.out.isoquant_transcript_gtf
        ch_versions = ch_versions.mix(ISOQUANT_CORRECTION.out.versions.first())
        GFFREAD_GETFASTA( ch_isoquant_gtf, ch_genome_fasta )
        ch_isoquant_transcripts = GFFREAD_GETFASTA.out.transcripts_fa
        ch_versions = ch_versions.mix(GFFREAD_GETFASTA.out.versions.first())
    }

    // SQANTI?

    //
    // Fusion gene detection
    // MODULE: JAFFAL
    if (!params.skip_jaffal && !params.custom_genome) {
        JAFFAL( ch_sample, ch_jaffal_ref )
        ch_jaffal_fasta = JAFFAL.out.jaffal_fasta
        ch_jaffal_csv = JAFFAL.out.jaffal_csv
        ch_versions = ch_versions.mix(JAFFAL.out.versions.first())
        }

    //
    // Transcriptome assessment
    // SQANTI, gffcompare
    //

    //
    // Transcript quantification
    // TransSigner
    if (!params.skip_quantification && !params.skip_mapping)
        TRANSIGNER_MAP
        TRANSIGNER_
        TRANSIGNER_QUANT

    //
    // Collate statistics
    //

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .unique()
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    /*
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

    */
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
