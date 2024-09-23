/*
----------------------------------------------------------------------------------------
    VALIDATE INPUTS
----------------------------------------------------------------------------------------
*/

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
include { CRAMINO               } from '../modules/local/cramino'
include { ALFRED               } from '../modules/local/alfred'
include { INPUT_CHECK          } from '../subworkflows/local/input_check'



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { NANOQ                } from '../modules/nf-core/nanoq/main'
//include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
//include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_directrna_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DIRECTRNA {

    //take:
    //ch_samplesheet // channel: samplesheet read in from --input

    //main:

    ch_versions = Channel.empty()
    //ch_multiqc_files = Channel.empty()

    // INPUT_CHECK
    INPUT_CHECK ( ch_input )
        .set { ch_sample }

    ch_sample.view()

    //ch_sample
    //    .map { it -> [ it[0], it[1], it[3] } // take sample, replicate, reads
    //    .set { ch_fastq }


    //
    // MODULE: Run nanoq
    //

    //if (!params.skip_qc) {


    //NANOQ ( ch_fastq )

    //)

    //}
    //ch_versions = ch_versions.mix(NANOQ.out.versions.first())


    //
    // MODULE: Run FastQC
    //
    //FASTQC (
    //    ch_samplesheet
    //)
    //ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    //ch_versions = ch_versions.mix(FASTQC.out.versions.first())


    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
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
}
    */

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
