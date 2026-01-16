/*
----------------------------------------------------------------------------------------
    VALIDATE INPUTS
----------------------------------------------------------------------------------------
*/

// nextflow magik

    //def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)
    // Create a new channel of metadata from a sample sheet passed to the pipeline through the --input parameter
    //ch_input = Channel.fromList(samplesheetToList(params.input, "assets/schema_input.json"))


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

if (params.skip_prepare_reference){
    ch_genome_minimap2_index        = Channel.fromPath(params.genome_minimap2_index, checkIfExists: true)
    ch_transcriptome_minimap2_index = Channel.fromPath(params.transcriptome_minimap2_index, checkIfExists: true)
    ch_genome_fasta_index           = Channel.fromPath(params.genome_fasta_index, checkIfExists: true)
    ch_genome_sizes                 = Channel.fromPath(params.genome_fasta_sizes, checkIfExists: true)
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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { softwareVersionsToYAML    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
//include { samplesheetToList } from 'plugin/nf-schema'
//include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
//include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_directrna_pipeline'

/*
----------------------------------------------------------------------------------------
    IMPORT LOCAL MODULES / SUBWORKFLOWS / FUNCTIONS
----------------------------------------------------------------------------------------
*/
// Check samplesheet
include { INPUT_CHECK               } from '../subworkflows/local/input_check'
// Prepare reference files
include { PREPARE_REFERENCE         } from '../subworkflows/local/prepare_reference'
// fastq QC
include { NANOQ                     } from '../modules/local/nanoq/nanoq'
include { SEQUALI                   } from '../modules/local/sequali/sequali'
include { MULTIQC                   } from '../modules/nf-core/multiqc/main'

// fastq mapping
include { MAPPING                   } from '../subworkflows/local/mapping'
// mapping visualisation
include { BAM_TO_BEDGRAPH as BAM_TO_BEDGRAPH_FW     } from '../modules/local/bedtools/bam_to_bedgraph'
include { BAM_TO_BEDGRAPH as BAM_TO_BEDGRAPH_REV    } from '../modules/local/bedtools/bam_to_bedgraph'
include { BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG as BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_FW } from '../subworkflows/local/bedgraph_bedclip_bedgraphtobigwig'
include { BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG as BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_REV } from '../subworkflows/local/bedgraph_bedclip_bedgraphtobigwig'

// bam QC
include { BAM_QC                    } from '../subworkflows/local/bam_qc'
include { SAMTOOLS_INDEX            } from '../modules/local/samtools/index'

// transcript reconstruction
include { BAM_TO_BED12              } from '../modules/local/flair/bam_to_bed12'
include { FLAIR_CORRECT             } from '../modules/local/flair/correct'
include { FLAIR_COLLAPSE            } from '../modules/local/flair/collapse'
include { BAMBU                     } from '../modules/local/bambu/bambu'
include { ISOQUANT                  } from '../modules/local/isoquant/isoquant'
include { GTF2DB                    } from '../modules/local/isoquant/gtf2db'
include { STRINGTIE                 } from '../modules/local/stringtie/stringtie'
include { GFFREAD_GETFASTA as GFFREAD_GETFASTA_BAMBU    } from '../modules/local/gffread/gffread'
include { GFFREAD_GETFASTA as GFFREAD_GETFASTA_ISOQUANT } from '../modules/local/gffread/gffread'
include { GFFREAD_GETFASTA as GFFREAD_GETFASTA_STRINGTIE} from '../modules/local/gffread/gffread'

// fusion gene detection
include { JAFFAL                    } from '../modules/local/jaffal/bpipe'

// transcript quantification
//include { MINIMAP2_TXOME_ALIGN as MINIMAP2_FLAIR        } from '../modules/local/minimap2_txome_align
//include { MINIMAP2_TXOME_ALIGN as MINIMAP2_BAMBU        } from '../modules/local/minimap2_txome_align
//include { MINIMAP2_TXOME_ALIGN as MINIMAP2_ISOQUANT     } from '../modules/local/minimap2_txome_align
//include { MINIMAP2_TXOME_ALIGN as MINIMAP2_STRINGTIE    } from '../modules/local/minimap2_txome_align
// OARFISH
include { OARFISH as OARFISH_FLAIR                      } from '../modules/local/oarfish/raw_read'
include { OARFISH as OARFISH_BAMBU                      } from '../modules/local/oarfish/raw_read'
include { OARFISH as OARFISH_ISOQUANT                   } from '../modules/local/oarfish/raw_read'
include { OARFISH as OARFISH_STRINGTIE                  } from '../modules/local/oarfish/raw_read'
// TRANSIGNER
//include { TRANSIGNER as TRANSIGNER_FLAIR                } from '../modules/local/transigner/align'

// transcriptome assessment
// JACCARD for tools using read correction
include { BEDTOOLS_JACCARD as BEDTOOLS_JACCARD_FLAIR    } from '../modules/local/bedtools/jaccard'
include { BEDTOOLS_JACCARD as BEDTOOLS_JACCARD_ISOQUANT } from '../modules/local/bedtools/jaccard'
// GFFCOMPARE
include { GFFCOMPARE as GFFCOMPARE_FLAIR          } from '../modules/local/gffcompare/gffcompare'
include { GFFCOMPARE as GFFCOMPARE_BAMBU          } from '../modules/local/gffcompare/gffcompare'
include { GFFCOMPARE as GFFCOMPARE_ISOQUANT       } from '../modules/local/gffcompare/gffcompare'
include { GFFCOMPARE as GFFCOMPARE_STRINGTIE      } from '../modules/local/gffcompare/gffcompare'

// Going to be a bit of a long-think
//include { SQANTI_PREPARE_REFERENCE                  } from '../subworkflows/local/sqanti/'
//include { QC_SQANTI as SQANTI_QC_FLAIR              } from '../subworkflows/local/sqanti/qc'
//include { QC_SQANTI as SQANTI_QC_BAMBU              } from '../subworkflows/local/sqanti/qc'
//include { QC_SQANTI as SQANTI_QC_ISOQUANT           } from '../subworkflows/local/sqanti/qc'
//include { QC_SQANTI as SQANTI_QC_STRINGTIE          } from '../subworkflows/local/sqanti/qc'
//include { SQANTI_FILTER as SQANTI_FILTER_FLAIR      } from '../subworkflows/local/sqanti/filter'
//include { SQANTI_FILTER as SQANTI_FILTER_BAMBU      } from '../subworkflows/local/sqanti/filter'
//include { SQANTI_FILTER as SQANTI_FILTER_ISOQUANT   } from '../subworkflows/local/sqanti/filter'
//include { SQANTI_FILTER as SQANTI_FILTER_STRINGTIE  } from '../subworkflows/local/sqanti/filter'
//include { SQANTI_RESCUE as SQANTI_RESCUE_FLAIR      } from '../subworkflows/local/sqanti/rescue'
//include { SQANTI_RESCUE as SQANTI_RESCUE_BAMBU      } from '../subworkflows/local/sqanti/rescue'
//include { SQANTI_RESCUE as SQANTI_RESCUE_ISOQUANT   } from '../subworkflows/local/sqanti/rescue'
//include { SQANTI_RESCUE as SQANTI_RESCUE_STRINGTIE  } from '../subworkflows/local/sqanti/rescue'

// transcript reconstruction subworkflow?
// include { TRANSCRIPT_RECONSTRUCTION

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DIRECTRNA{

    //take:
    //main:

    ch_versions = Channel.empty()
    ch_versions.view()
    ch_multiqc_files = Channel.empty()
    //def multiqc_report      = []

    // INPUT_CHECK
    INPUT_CHECK ( ch_input )
        .set { ch_sample }

    ch_sequencing_type = channel.value(params.sequencing_type)

    // QC of fastq files
    /// MODULES: NANOQ & SEQUALI
    if (!params.skip_qc) {
        if (!params.bam_input) {
            if (!params.skip_nanoq) {
                NANOQ( ch_sample )
                ch_nanoq_stats = NANOQ.out.stats.collect{it[1]}.flatten()
                ch_multiqc_files = ch_multiqc_files.mix(ch_nanoq_stats.ifEmpty([]))
                ch_versions = ch_versions.mix(NANOQ.out.versions.first())
            }
            if (!params.skip_sequali) {
                SEQUALI( ch_sample )
                ch_sequali_stats = SEQUALI.out.sequali_json.collect{it[1]}.flatten()
                ch_multiqc_files = ch_multiqc_files.mix(ch_sequali_stats.ifEmpty([]))
                ch_versions = ch_versions.mix(SEQUALI.out.versions.first())
            }
        }
    }

    // Prepare the reference files
    /// SUBWORKFLOW: PREPARE_REFERENCE
    if (!params.skip_prepare_reference) {
        PREPARE_REFERENCE (
            params.genome_fasta,
            params.genome_fasta_index,
            params.genome_fasta_sizes,
            params.genome_minimap2_index,
            params.bam_input,
            params.transcriptome_fasta,
            params.transcriptome_minimap2_index,
            params.annotation_gtf,
            params.skip_jaffal,             // boolean [default: false]
            params.skip_jaffal_download,    // boolean [default: false]
            params.jaffal_reference,        // path
            params.skip_transcript_quantification, // boolean [default: false]
            params.skip_sqanti_all,         // boolean [default: false]
            params.skip_sqanti_qc,          // boolean [defeault: false]
            params.sqanti_qc_reference,     // value: human, mouse or custom
            params.sqanti_qc_cage,          // boolean [default: true]
            params.sqanti_qc_cage_path,          // boolean [default: path]
            params.sqanti_qc_polyA_sites,   // boolean [default: true]
            params.sqanti_qc_polyA_sites_path,   // boolean [default: path]
            params.sqanti_qc_polyA_motif,   // boolean [default: true]
            params.sqanti_qc_polyA_motif_path,   // boolean [default: path]
            params.sqanti_qc_intron_junctions, // boolean [default: true]
            params.sqanti_qc_intron_path // boolean [default: path]
        )
        // initialize genome + transcriptome references
        ch_genome_fasta                 = PREPARE_REFERENCE.out.genome_fasta
        ch_genome_fasta_index           = PREPARE_REFERENCE.out.genome_fasta_index
        ch_genome_sizes                 = PREPARE_REFERENCE.out.genome_fasta_sizes
        ch_genome_minimap2_index        = PREPARE_REFERENCE.out.genome_minimap2_index
        ch_transcriptome_fasta          = PREPARE_REFERENCE.out.transcriptome_fasta
        ch_transcriptome_minimap2_index = PREPARE_REFERENCE.out.transcriptome_minimap2_index
        ch_annotation_gtf               = PREPARE_REFERENCE.out.annotation_gtf
        ch_jaffal_reference_dir         = PREPARE_REFERENCE.out.jaffal_reference
        // Combine genome fasta with genome fasta index into single channel -
        // some software expect both files in a single path/channel
        ch_genome_fasta_with_index = ch_genome_fasta.combine(ch_genome_fasta_index)

        // initialize sqanti qc references
        if (!params.skip_sqanti_qc) {
            if (params.sqanti_qc_cage) {
                ch_sqanti_qc_cage_bed = PREPARE_REFERENCE.out.sqanti_qc_cage_bed
            }
            if (params.sqanti_qc_polyA_sites) {
                ch_sqanti_qc_polyA_sites_bed = PREPARE_REFERENCE.out.sqanti_qc_polyA_sites_bed
            }
            if (params.sqanti_qc_polyA_motif) {
                ch_sqanti_qc_polyA_motif = PREPARE_REFERENCE.out.sqanti_qc_polyA_motif
            }
            if (params.sqanti_qc_intron_junctions) {
                ch_sqanti_qc_intron_junctions_bed = PREPARE_REFERENCE.out.sqanti_qc_intron_junctions_bed
            }
        }
    } else {
        ch_genome_fasta_with_index = ch_genome_fasta.combine(ch_genome_fasta_index)
    }

    // Mapping and sorting
    // SUBWORKFLOW: MAPPING
    // TODO mapping to transcriptome
    if (!params.bam_input) {
        MAPPING( ch_sample, ch_genome_minimap2_index,
        ch_transcriptome_minimap2_index )
        ch_bam = MAPPING.out.bam
        ch_bam_index = MAPPING.out.bai
        ch_bam_index_path = MAPPING.out.bai.flatten().last()
        ch_mixed_bam = ch_bam.combine(ch_bam_index_path)
        ch_versions = ch_versions.mix(MAPPING.out.versions.first())
    } else {
        ch_bam = ch_sample
        SAMTOOLS_INDEX( ch_bam )
        ch_bam_index = SAMTOOLS_INDEX.out.bai
        ch_mixed_bam = ch_bam.combine(ch_bam_index)
    }

    // BAM TO BIGWIG for visualisation
    // uses SUBWORKFLOW: BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG
    if (!params.skip_bam_to_bigwig) {
        BAM_TO_BEDGRAPH_FW( ch_bam, ch_genome_sizes, '+' )
        BAM_TO_BEDGRAPH_REV( ch_bam, ch_genome_sizes, '-' )
        ch_bedgraph_fw = BAM_TO_BEDGRAPH_FW.out.bedgraph
        ch_bedgraph_rev = BAM_TO_BEDGRAPH_REV.out.bedgraph
        BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_FW( ch_bedgraph_fw, ch_genome_sizes, '+' )
        BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_REV( ch_bedgraph_rev, ch_genome_sizes, '-' )
    }

    // BAM QC
    // SUBWORKFLOW: BAM_QC
    // Execute cramino, alfred, samtools flagstat, ngs-bits on bam output from mapping
    if (!params.skip_bam_qc) {
        ch_skip_cramino = params.skip_cramino
        ch_skip_alfred = params.skip_alfred
        ch_skip_samtools_flagstat = params.skip_samtools_flagstat
        ch_cramino_min_length = params.cramino_min_length
        ch_skip_ngs_bits = params.skip_ngs_bits
        ch_ngs_bits_build = params.ngs_bits_build
        //TO-DO NEED TO ADD THIS CONTAMINATION TO A CUSTOM CONFIG
        ch_ngs_bits_skip_contamination = params.ngs_bits_skip_contamination
        BAM_QC(
            ch_skip_cramino,
            ch_skip_alfred,
            ch_skip_samtools_flagstat,
            ch_skip_ngs_bits,
            ch_ngs_bits_build,
            ch_ngs_bits_skip_contamination,
            ch_bam,
            ch_mixed_bam,
            ch_genome_fasta,
            ch_cramino_min_length
            )
        if (!params.skip_ngs_bits){
            ch_ngs_bits_stats = BAM_QC.out.ngs_bits_stats.collect{it[1]}.flatten()
            ch_multiqc_files = ch_multiqc_files.mix(ch_ngs_bits_stats.ifEmpty([]))
        }
        if (!params.skip_samtools_flagstat){
            ch_samtools_flagstat_stats = BAM_QC.out.flagstat.collect{it[1]}.flatten()
            ch_multiqc_files = ch_multiqc_files.mix(ch_samtools_flagstat_stats.ifEmpty([]))
        }
        ch_versions = ch_versions.mix(BAM_QC.out.versions)
    }


    // Stand alone read correction tools? Which ones....

    // TRANSCRIPT RECONSTRUCTION
    //
    // FLAIR
    if (!params.skip_flair) {
        if (!params.skip_flair_correct) {
            BAM_TO_BED12( ch_bam, ch_bam_index )
            ch_mapped_bed = BAM_TO_BED12.out.bed
            FLAIR_CORRECT( ch_mapped_bed, ch_genome_fasta, ch_annotation_gtf )
            ch_flair_corrected_bed = FLAIR_CORRECT.out.flair_corrected_bed
            BEDTOOLS_JACCARD_FLAIR( ch_flair_corrected_bed, ch_mapped_bed, 'flair' )
        }
        if (!params.skip_flair_collapse) {
            if (!params.skip_flair_correct) {
                FLAIR_COLLAPSE( ch_sample, ch_flair_corrected_bed, ch_annotation_gtf, ch_genome_fasta )
                ch_flair_collapsed_bed = FLAIR_COLLAPSE.out.collapsed_isoforms_bed
                ch_flair_collapsed_gtf = FLAIR_COLLAPSE.out.collapsed_isoforms_gtf
                ch_flair_collapsed_fa = FLAIR_COLLAPSE.out.collapsed_isoforms_fa
            } else {
                BAM_TO_BED12( ch_bam, ch_bam_index )
                ch_mapped_bed = BAM_TO_BED12.out.bed
                FLAIR_COLLAPSE( ch_sample, ch_mapped_bed, ch_annotation_gtf, ch_genome_fasta )
                ch_flair_collapsed_bed = FLAIR_COLLAPSE.out.collapsed_isoforms_bed
                ch_flair_collapsed_gtf = FLAIR_COLLAPSE.out.collapsed_isoforms_gtf
                ch_flair_collapsed_fa = FLAIR_COLLAPSE.out.collapsed_isoforms_fa
            }
        }
    }

    // BAMBU
    if (!params.skip_bambu) {
        BAMBU( ch_genome_fasta, ch_annotation_gtf, ch_bam )
        ch_bambu_supported_gtf = BAMBU.out.bambu_supported_gtf
        ch_versions = ch_versions.mix(BAMBU.out.versions.first())
        GFFREAD_GETFASTA_BAMBU( ch_bambu_supported_gtf, ch_genome_fasta_with_index, 'bambu' )
        ch_bambu_transcripts = GFFREAD_GETFASTA_BAMBU.out.transcripts_fa
        ch_versions = ch_versions.mix(GFFREAD_GETFASTA_BAMBU.out.versions.first())
        }

    // ISOQUANT
    if (!params.skip_isoquant) {
        GTF2DB( ch_annotation_gtf )
        ch_isoquant_database = GTF2DB.out.isoquant_database
        ISOQUANT( ch_mixed_bam, ch_isoquant_database, ch_genome_fasta)
        ch_isoquant_gtf = ISOQUANT.out.isoquant_transcript_gtf
        ch_versions = ch_versions.mix(ISOQUANT.out.versions.first())
        GFFREAD_GETFASTA_ISOQUANT( ch_isoquant_gtf, ch_genome_fasta_with_index, 'isoquant' )
        ch_isoquant_transcripts = GFFREAD_GETFASTA_ISOQUANT.out.transcripts_fa
        ch_versions = ch_versions.mix(GFFREAD_GETFASTA_ISOQUANT.out.versions.first())
    }

    // STRINGTIE
    if (!params.skip_stringtie) {
        STRINGTIE( ch_bam, ch_annotation_gtf )
        ch_stringtie_gtf = STRINGTIE.out.stringtie_gtf
        ch_versions = ch_versions.mix(STRINGTIE.out.versions.first())
        GFFREAD_GETFASTA_STRINGTIE( ch_stringtie_gtf, ch_genome_fasta_with_index, 'stringtie' )
        ch_stringtie_transcripts = GFFREAD_GETFASTA_STRINGTIE.out.transcripts_fa
        ch_versions = ch_versions.mix(GFFREAD_GETFASTA_STRINGTIE.out.versions.first())
    }

    // Fusion gene detection
    // MODULE: JAFFAL
    if (!params.skip_jaffal && !params.custom_genome) {
        JAFFAL( ch_sample, ch_jaffal_reference_dir )
        ch_jaffal_fasta = JAFFAL.out.jaffal_fasta
        ch_jaffal_csv = JAFFAL.out.jaffal_results
        ch_versions = ch_versions.mix(JAFFAL.out.versions.first())
        }

    //
    // Transcriptome assessment
    // gffcompare
    if (!params.skip_gffcompare) {
        if (!params.skip_flair) {
            GFFCOMPARE_FLAIR( ch_genome_fasta_with_index, ch_flair_collapsed_gtf, ch_annotation_gtf, 'flair' )
            ch_flair_gffcompare_stats = GFFCOMPARE_FLAIR.out.gffcompare_stats.collect{it[1]}.flatten()
            ch_multiqc_files = ch_multiqc_files.mix(ch_flair_gffcompare_stats.ifEmpty([]))
        }
        if (!params.skip_bambu) {
            GFFCOMPARE_BAMBU( ch_genome_fasta_with_index, ch_bambu_supported_gtf, ch_annotation_gtf, 'bambu' )
            ch_bambu_gffcompare_stats = GFFCOMPARE_BAMBU.out.gffcompare_stats.collect{it[1]}.flatten()
            ch_multiqc_files = ch_multiqc_files.mix(ch_bambu_gffcompare_stats.ifEmpty([]))
        }
        if (!params.skip_isoquant) {
            GFFCOMPARE_ISOQUANT(ch_genome_fasta_with_index, ch_isoquant_gtf, ch_annotation_gtf, 'isoquant' )
            ch_isoquant_gffcompare_stats = GFFCOMPARE_ISOQUANT.out.gffcompare_stats.collect{it[1]}.flatten()
            ch_multiqc_files = ch_multiqc_files.mix(ch_isoquant_gffcompare_stats.ifEmpty([]))
        }
        if (!params.skip_stringtie) {
            GFFCOMPARE_STRINGTIE( ch_genome_fasta_with_index, ch_stringtie_gtf, ch_annotation_gtf, 'stringtie' )
            ch_stringtie_gffcompare_stats = GFFCOMPARE_STRINGTIE.out.gffcompare_stats.collect{it[1]}.flatten()
            ch_multiqc_files = ch_multiqc_files.mix(ch_stringtie_gffcompare_stats.ifEmpty([]))
        }
    }

/*
    // TODO
    // if (!skip_sqanti_all) {
        if (!skip_sqanti_qc) {
            if (run_flair){
                SQANTI_QC_FLAIR( ch_flair_collapsed_gtf, ch_annotation_gtf, ch_genome_fasta_with_index, 'flair' )
            }
            if (run_bambu){
                SQANTI_QC_BAMBU( ch_bambu_gtf, ch_annotation_gtf, ch_genome_fasta_with_index, 'bambu' )
            }
            if (run_isoquant){
                SQANTI_QC_ISOQUANT( ch_isoquant_gtf, ch_annotation_gtf, ch_genome_fasta_with_index, 'isoquant' )
            }
            if (run_stringtie){
                SQANTI_QC_STRINGTIE( ch_stringtie_gtf, ch_annotation_gtf, ch_genome_fasta_with_index, 'stringtie' )
            }
        }
*/

    //
    // Transcript quantification
    // Oarfish
   // if (!params.skip_quantification && !params.skip_mapping && params.!skip_oarfish)
    if (!params.skip_transcript_quantification) {
        if (!params.skip_oarfish) {
            if (!params.skip_flair) {
                OARFISH_FLAIR( ch_flair_collapsed_fa, ch_transcriptome_minimap2_index, ch_sequencing_type, 'flair' )
            }
            if (!params.skip_bambu) {
                OARFISH_BAMBU( ch_bambu_transcripts, ch_transcriptome_minimap2_index, ch_sequencing_type, 'bambu' )
            }
            if (!params.skip_isoquant) {
                OARFISH_ISOQUANT( ch_isoquant_transcripts, ch_transcriptome_minimap2_index, ch_sequencing_type, 'isoquant' )
            }
            if (!params.skip_stringtie) {
                OARFISH_STRINGTIE( ch_stringtie_transcripts, ch_transcriptome_minimap2_index, ch_sequencing_type, 'stringtie' )
            }
        }
    }

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

    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    /* summary_params      = paramsSummaryMap(
    //    workflow, parameters_schema: "nextflow_schema.json")
    ///ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    //ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
    //    file(params.multiqc_methods_description, checkIfExists: true) :
    //    file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    //ch_methods_description                = Channel.value(
    //    methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
    ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )
*/

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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
