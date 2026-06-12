process ISOQUANT {
    tag "$meta.id"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/isoquant:3.7.0--hdfd78af_0' :
        'biocontainers/isoquant:3.7.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bam_bai)
    path annotation_gtf
    path genome_fasta_index

    output:
    tuple val(meta), path("*/*.corrected_reads.bed.gz")           , optional: true
    tuple val(meta), path("*/*.discovered_gene_counts.tsv")       , optional: true
    tuple val(meta), path("*/*.discovered_gene_tpm.tsv")          , optional: true
    tuple val(meta), path("*/*.discovered_transcript_counts.tsv") , optional: true
    tuple val(meta), path("*/*.discovered_transcript_tpm.tsv")    , optional: true
    tuple val(meta), path("*/*.extended_annotation.gtf")          , emit: isoquant_new_reference_transcriptome_gtf , optional: true
    tuple val(meta), path("*/*.gene_counts.tsv")
    tuple val(meta), path("*/*.gene_tpm.tsv")                     , optional: true
    tuple val(meta), path("*/*.read_assignments.tsv.gz")          , optional: true
    tuple val(meta), path("*/*.transcript_counts.tsv")
    tuple val(meta), path("*/*.transcript_model_reads.tsv.gz")    , emit: isoquant_transcript_models
    tuple val(meta), path("*/*.transcript_models.gtf")            , emit: isoquant_transcript_gtf
    tuple val(meta), path("*/*.transcript_tpm.tsv")               , optional: true
    //tuple val(meta)   , path("*/*.transcript_model_tpm.tsv")         , optional: true
    //tuple val(meta)   , path("*/*.transcript_model_counts.tsv")      , optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}_${meta.replicate}_isoquant"
    def input_bam   = task.ext.input_bam ?: "--bam $bam"
    def ref_genome  = task.ext.ref_genome ?: "--reference $genome_fasta_index"
    def ref_gtf     = task.ext.ref_gtf ?: "--genedb $annotation_gtf"
    //def output = task.ext.output ?: "--output isoquant_${meta.id}_${meta.replicate}"
    """
    export HOME=\$(pwd)

    isoquant.py \\
        $args \\
        $input_bam \\
        $ref_genome \\
        $ref_gtf \\
        -o . \\
        --prefix $prefix \\
        --threads $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isoquant: \$(isoquant.py -v | sed 's#IsoQuant ##')
    END_VERSIONS
    """

    stub:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}_${meta.replicate}_isoquant"
    def input_bam  = task.ext.input_bam ?: "--bam $bam"
    def ref_genome = task.ext.ref_genome ?: "--reference $genome_fasta_index"
    def ref_gtf    = task.ext.ref_gtf ?: "--genedb $annotation_gtf"
    //def output   = task.ext.output ?: "--output isoquant_${meta.id}_${meta.replicate}"
    //touch ${prefix}.bam

    """
    touch ${prefix}.read_assignments.tsv.gz
    touch ${prefix}.corrected_reads.bed.gz
    touch ${prefix}.gene_counts.tsv
    touch ${prefix}.transcript_counts.tsv
    touch ${prefix}.gene_tpm.tsv
    touch ${prefix}.transcript_tpm.tsv
    touch ${prefix}.transcript_models.gtf
    touch ${prefix}.extended_annotation.gtf
    touch ${prefix}.transcript_model_reads.tsv.gz
    touch ${prefix}.transcript_model_tpm.tsv
    touch ${prefix}.transcript_model_counts.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isoquant: \$(isoquant.py -v | sed 's#IsoQuant ##')
    END_VERSIONS
    """
}
