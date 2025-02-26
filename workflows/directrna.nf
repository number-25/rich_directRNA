/*
----------------------------------------------------------------------------------------
    VALIDATE INPUTS
----------------------------------------------------------------------------------------
*/

// nextflow magik

//def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)
    // nf-schema plugins
    // Not working yet, but is promising for testing further on - it could avoid
    // using the Julia script to validate inputs - consider embarking on this
    // once a stable release it pushed/ Validate input parameters()
    // https://nextflow-io.github.io/nf-schema/latest/parameters/help_text/
    //validateParameters()

    // Print summary to stdout of supplied parameters that differ from defaults
    log.info paramsSummaryLog(workflow)

    // Create a new channel of metadata from a sample sheet passed to the pipeline through the --input parameter
    //ch_input = Channel.fromList(samplesheetToList(params.input, "assets/schema_input.json"))


//

// Check mandatory parameters (missing protocol or profile will exit the run.)
// inputs samplesheet.csv
if (params.input) {
    ch_input = file(params.input) // defined in nextflow.config
} else {
    exit 1, 'Input samplesheet not specified!'
}

// genome fasta
if (params.genome_fasta) {
    ch_genome_fasta = Channel.fromPath(params.genome_fasta, checkIfExists: true)
} else {
    exit 1, 'Reference genome fasta file is not specified! please modify nextflow.config or use --genome_fasta parameter'
}

// transcriptome
if (params.annotation_gtf) {
    ch_annotation_gtf = Channel.fromPath(params.annotation_gtf, checkIfExists: true) // check if exists
} else {
    exit 1, 'Reference transcriptome annotation file is not specified! please modify nextflow.config or use --annotation_gtf parameter'
}

if (params.transcriptome_fasta) {
    ch_transcriptome_fasta = Channel.fromPath(params.transcriptome_fasta, checkIfExists: true) // check if exists
} else {
    exit 1, 'Reference transcriptome fasta file is not specified! please modify nextflow.config or use --transcriptome_fasta parameter'
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
include { GUNZIP as GUNZIP_FASTA    } from '../modules/nf-core/gunzip'
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

// Going to be a bit of a long-think
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
//include { samplesheetToList } from 'plugin/nf-schema'
include { paramsSummaryLog } from 'plugin/nf-schema'
//include { validateParameters } from 'plugin/nf-schema'
//include { paramsSummaryMap       } from 'plugin/nf-schema'
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
    ///ch_sample
    ///.map { it -> [ it[0], it[1] ] } // take sample, replicate, reads
    ///.set { ch_fastq }

    // QC of fastq files
    /// Toulligqc?
    /// MODULE: NANOQ
    if (!params.skip_qc) {
        NANOQ ( ch_sample )
        ch_versions = ch_versions.mix(NANOQ.out.versions.first())
    }
    // MODULE: Run FastQC
    ///
    ///)
    ///ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ///ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    // Prepare the reference files
    // Only for SQANTI3 at this moment
    /// SUBWORKFLOW: PREPARE_REFERENCE
    /// Only used when a completely custom workflow is being run
    if (!params.skip_prepare_reference) {
        PREPARE_REFERENCE (
        params.genome_fasta,
        params.genome_fasta_index,
        params.genoma_fasta_sizes,
        params.annotation_gtf
        params.skip_jaffal
        params.skip_jaffal_download
        params.skip_sqanti_all
        params.skip_sqanti_qc
        params.sqanti_qc_download
        params.sqanti_qc_reference
        params.sqanti_qc_polyA_sites
        params.sqanti_qc_polyA_motif
        params.sqanti_qc_intron_junctions
        //params.cage_bed,
        //params.polyA_bed,
        //params.polyA_sites,
        //params.intropolis_bed,
        //params.skip_sqanti_qc
        )
        ch_genome_fasta = PREPARE_REFERENCE.out.genome_fasta
        ch_genome_index =
        ch_genome_sizes =
        ch_genome_minimap2_index
        if (!params.skip_sqanti_qc)

        ch_cage_bed = PREPARE_REFERENCE.out.cage_bed
        ch_polyA_bed = PREPARE_REFERENCE.out.polyA_bed
        ch_polyA_sites = PREPARE_REFERENCE.out.polyA_sites
        //ch_polyA_sites = Channel.fromPath(polyA_sites, checkIfExists = true)
        ch_intropolis_bed = PREPARE_REFERENCE.out.intropolis_bed
        ch_versions = ch_versions.mix(PREPARE_REFERENCE.out.versions)
        ch_jaffal_ref = Channel.fromPath(jaffal_ref, checkIfExists = true)
    }

    // Mapping and sorting
    // SUBWORKFLOW: MAPPING
    // if the reference preparation is skipped and a reference file isn't provided, minimap2 will generate an index
    // for the provided reference file
    //
    if (!params.skip_mapping) {
        if (params.genome_fasta.endsWith('.gz')) {
            GUNZIP_FASTA( ch_genome_fasta )
            ch_genome_fasta = GUNZIP_FASTA.out.gunzip
            ch_versions = ch_versions.mix(GUNZIP.out.versions.first())
        } else {
            ch_genome_fasta = ch_genome_fasta
        }
        // If custom genome is provided --custom_genome=true
        if (params.custom_genome) {
            CUSTOM_GETCHROMSIZES( ch_genome_fasta )
            ch_genome_index = CUSTOM_GETCHROMSIZES.out.fai
            ch_genome_sizes = CUSTOM_GETCHROMSIZES.out.sizes
            ch_versions = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions.first())
            MINIMAP2_INDEX( ch_genome_fasta )
            ch_genome_minimap2_index = MINIMAP2_INDEX.out.index
        } else {
            ch_genome_index = Channel.fromPath(params.genome_fasta_index, checkIfExists: true)
            if (!params.genome_fasta_minimap2_index) {
                ch_genome_minimap2_index = Channel.fromPath(params.genome_fasta_minimap2_index, checkIfExists: true)
            ch_genome_minimap2_index = Channel.fromPath(params.genome_fasta_minimap2_index, checkIfExists: true)
            ch_genome_sizes = Channel.fromPath(params.genome_fasta_sizes, checkIfExists: true)
            }
        MAPPING( ch_sample, ch_genome_fasta, ch_genome_minimap2_index )
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
        ch_versions = ch_versions.mix(BAM_QC.out.versions)
    // If a raw BAM is provided and mapping is not needed
    // May need to index it to work with downstream processes? Needs testing
    } else {
        //ch_bam = ch_sample
        //ch_bam_index =
        BAM_QC( ch_sample, ch_genome_fasta )
        ch_versions = ch_versions.mix(BAM_QC.out.versions)
        }

    //
    // TRANSCRIPT RECONSTRUCTION
    //

    if (!params.skip_flair_correct && !params.skip_flair_collapse) {
        BAM_TO_BED12( ch_bam, ch_bam_index )
        // seeing if a mixed channel with bam and bam.bai works
        //BAM_TO_BED12( ch_mixed_bam )
        ch_mapped_bed = BAM_TO_BED12.out.bed
        /*FLAIR_CORRECT( ch_mapped_bed, ch_genome_fasta, ch_annotation_gtf )
        ch_corrected_bed = FLAIR_CORRECT.out.flair_corrected_bed
        FLAIR_COLLAPSE( ch_corrected_bed, ch_sample, ch_annotation_gtf, ch_genome_fasta )
        ch_collapsed_bed = FLAIR_COLLAPSE.out.collapsed_isoforms_bed
        ch_collapsed_gtf = FLAIR_COLLAPSE.out_collapsed_isoforms.gtf
        ch_versions = ch_versions.mix(FLAIR_collapse.out.versions)
        */
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
        // No correcton just collapsing
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
        }
    }

    // Read correction tools? Which ones....
    // TC-CLEAN?
    // IsoQUANT?
    // FLAIR

    // BAMBU
    //if (!params.skip_bambu) {
    //    BAMBU( ch_genome_fasta, ch_annotation_gtf, ch_bam )
    //    ch_bambu_gtf = BAMBU.out.bambu_extended_gtf
    //    ch_versions = ch_versions.mix(BAMBU.out.versions.first())
        // MIX genome fasta with fasta index as this will improve GFFREADs speed
        //GFFREAD_GETFASTA( ch_bambu_gtf, ch_genome_fasta )
        //ch_bambu_transcripts = GFFREAD_GETFASTA.out.transcripts_fa
        //ch_versions = ch_versions.mix(GFFREAD_GETFASTA.out.versions.first())
        }
/*
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

    //
    // Transcript quantification
    // TransSigner
    //if (!params.skip_quantification && !params.skip_mapping)
    //    TRANSIGNER_MAP
    //    TRANSIGNER_
    //    TRANSIGNER_QUANT


    //
    // Fusion gene detection
    // MODULE: JAFFAL
    if (!params.skip_jaffal && !params.custom_genome) {
        JAFFAL( ch_sample, ch_jaffal_ref )
        ch_jaffal_fasta = JAFFAL.out.jaffal_fasta
        ch_jaffal_csv = JAFFAL.out.jaffal_csv
        ch_versions = ch_versions.mix(JAFFAL.out.versions.first())
        }
*/
    //
    // Transcriptome assessment
    // SQANTI, gffcompare
    // Different SQANTI options - create as subworkflow?

    // Not done yet
    if (!skip gff_compare) {
        GFFCOMPARE( ch_reconstructed_gtf, ch_annotation_gtf, ch_genome_fasta )
    }

    if (!skip_sqanti_qc) {
        if (!skip_sqanti


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
