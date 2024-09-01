//               https://github.com/nf-core/modules/tree/master/subworkflows
//
// Uncompress and prepare reference genome files
//

include { GUNZIP as GUNZIP_FASTA } from '../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GTF } from '../modules/nf-core/gunzip'
//include { GUNZIP as GUNZIP_BED } from '../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_TRANSCRIPT_FASTA } from
include { GUNZIP as GUNZIP_CAGE } from '../modules/nf-core/gunzip'
'../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_POLYA } from '../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_INTROPOLIS } from '../modules/nf-core/gunzip'

// prepare indices for reference

include { MINIMAP2_INDEX as MINIMAP2_GENOME_INDEX } from '../modules/nf-core/minimap2/index'
include { MINIMAP2_INDEX as MINIMAP2_TRANSCRIPTOME_INDEX } from '../modules/nf-core/minimap2/index'
include { SAMTOOLS_FAIDX } from '../modules/nf-core/samtools/faidx'
include { CUSTOM_GETCHROMSIZES } from
'../modules/nf-core/custom/getchromsizes'
include { SAMTOOLS_INDEX } from
'../modules/nf-core/samtools/index'


// prepare additional files

//TO-DO make these modules
include { GTF_TO_BED } from '../modules/local/gfx2bed' // gxf2bed module
include { BIGWIG_TO_WIG } from '../modules/local/bigwigtowig'
include { WIG_TO_BED } from '../modules/local/bedops'
include { JAFFAL_PREPARE_REF } from
'../modules/local/jaffal_prepare_ref'


include { SAMTOOLS_SORT      } from '../../../modules/nf-core/samtools/sort/main'

workflow PREPARE_REFERENCE {

    take:
    genome_fasta
    transcriptome_fasta
    annotation_gtf
    cage_bed
    polyA_bed
    intropolis_bed
    phylop_bigwig

    minimap2_index
    minimap2_transcriptome_index
    samtools_genome_index
    custom_chrom_sizes

    appris_bed?
    mane_select_bed?
    mane_clinical_bed?
    // TODO nf-core: edit input (take) channels

    main:

    ch_versions = Channel.empty()

    // Uncompress genome fasta file if required
    //
    if (genome_fasta.endsWith('.gz')) {
        ch_genome_fasta    = GUNZIP_FASTA ( [ [:], genome_fasta ] ).gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_genome_fasta = Channel.value(file(genome_fasta))
    }

    // Uncompress transcriptome fasta file if required
    //
    if (transcriptome_fasta.endsWith('.gz')) {
        ch_transcriptome_fasta    = GUNZIP_TRANSCRIPT_FASTA ( [ [:], transcriptome_fasta ] ).gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_TRANSCRIPT_FASTA.out.versions)
    } else {
        ch_transcriptome_fasta =
        Channel.value(file(transcriptome_fasta))
    }

    //
    // Uncompress GTF annotation file
    //
    if (annotation_gtf) {
        if (gtf.endsWith('.gz')) {
            ch_annotation_gtf      = GUNZIP_GTF ( [ [:],
            annotation_gtf ] ).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
        } else {
            ch_annotation_gtf = Channel.value(file(annotation_gtf))
        }
            }

    //
    // Uncompress cage BED interval file
    //
    if (cage_bed) {
        if (cage_bed.endsWith('.gz')) {
            ch_cage_bed = GUNZIP_CAGE ( [ [:], cage_bed ] ).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_CAGE.out.versions)
        } else {
            ch_cage_bed = Channel.value(file(cage_bed))
        }
    }

    //
    // Uncompress polyA BED interval file
    //
    if (polyA_bed) {
        if (polyA_bed.endsWith('.gz')) {
            ch_polyA_bed = GUNZIP_POLYA ( [ [:], polyA_bed ] ).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_POLYA.out.versions)
        } else {
            ch_polyA_bed = Channel.value(file(polyA_bed))
        }
    }

    //
    // Uncompress intropolis BED interval file
    //
    if (intropolis_bed) {
        if (intropolis_bed.endsWith('.gz')) {
            ch_intropolis_bed = GUNZIP_INTROPOLIS ( [ [:], polyA_bed ] ).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_INTROPOLIS.out.versions)
        } else {
            ch_intropolis_bed = Channel.value(file(intropolis_bed))
        }
    }

    //
    // Uncompress intropolis BED interval file
    //
    if (intropolis_bed) {
        if (intropolis_bed.endsWith('.gz')) {
            ch_intropolis_bed = GUNZIP_INTROPOLIS ( [ [:], polyA_bed ] ).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_INTROPOLIS.out.versions)
        } else {
            ch_intropolis_bed = Channel.value(file(intropolis_bed))
        }
    }

//-------------------------------------------------------//
//              Prepare indices                         //

    //
    // Create genome index for minimap2
    //

    ch_minimap2_index = MINIMAP2_GENOME_INDEX( ch_genome_fasta).index
    ch_versions = MINIMAP2_GENOME_INDEX( ch_genome_fasta).versions

    //
    // Create transcriptome index for minimap2
    //
    ch_minimap2_transcriptome_index = MINIMAP2_TRANSCRIPTOME_INDEX(
    ch_transcriptome_fasta.out.index)
    ch_versions = ch_versions.mix(MINIMAP2_TRANSCRIPTOME_INDEX.out.versions)

    //
    // Create samtools genome index and create custom chrom sizes
    //
    ch_samtools_genome_index =  CUSTOM_GETCHROMSIZES(ch_genome_fasta).fai
    ch_custom_chrom_sizes = CUSTOM_GETCHROMSIZES(ch_genome_fasta).sizes
    ch_versions = CUSTOM_GETCHROMSIZES(ch_genome_fasta).versions    //


    SAMTOOLS_SORT ( ch_bam )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:
    genome_fasta = ch_genome_fasta
    transcriptome_fasta = ch_transcriptome_fasta
    annotation_gtf = ch_annotation_gtf
    annotation_bed =
    cage_bed = ch_cage_bed
    polyA_bed = ch_polyA_bed
    intropolis_bed = ch_intropolis_bed

    minimap2_index = ch_minimap2_index
    minimap2_transcriptome_index = ch_minimap2_transcriptome_index
    samtools_genome_index = ch_samtools_genome_index
    custom_chrom_sizes = ch_custom_chrom_sizes
    jaffal_ref =


    // TODO nf-core: edit emitted channels
    bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

https://raw.githubusercontent.com/nf-core/rnaseq/master/subworkflows/local/prepare_genome/main.nf
