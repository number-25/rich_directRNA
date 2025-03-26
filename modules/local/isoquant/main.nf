process ISOQUANT {
    tag "$meta.id"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/isoquant:3.6.2--hdfd78af_0':
        'biocontainers/isoquant:3.6.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam)
    path annotation_gtf
    path genome_fasta

    output:
    tuple val(meta), path("*.read_assignments.tsv.gz")
    tuple val(meta), path("*.corrected_reads.bed.gz")
    tuple val(meta), path("*.transcript_tpm.tsv")
    tuple val(meta), path("*.transcript_counts.tsv")
    tuple val(meta), path("*.gene_tpm.tsv")
    tuple val(meta), path("*.gene_counts.tsv")
    tuple val(meta), path("*.transcript_models.gtf"), emit: isoquant_transcript_gtf
    tuple val(meta), path("*.transcript_model_reads.tsv.gz"), emit: isoquant_transcript_models
    tuple val(meta), path("*.transcript_model_tpm.tsv")
    tuple val(meta), path("*.transcript_model_counts.tsv")
    tuple val(meta), path("*.extended_annotation.gtf"), emit: isoquant_new_reference_transcriptome_gtf, optional: true
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.replicate}_isoquant"
    //def dRNA_preset = task.ext.dRNA_preset ?: "-d nanopore"
    //def strand_preset = task.ext.dRNA_preset ?: "--stranded none"
    def input_bam = task.ext.input_bam ?: "--bam $bam"
    def ref_genome = task.ext.ref_genome ?: "--reference $genome_fasta"
    def ref_gtf = task.ext.ref_genome ?: "--genedb $annotation_gtf"
    //def complete = task.ext.kmer ?: "--complete_genedb"
    //def output = task.ext.kmer ?: "--output $annotation_gtf"
    """
    isoquant.py \\
        $args \\
        $input_bam \\
        $ref_genome \\
        $ref_gtf \\
        --prefix $prefix \\
        --threads $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isoquant: \$(isoquant --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.replicate}_isoquant"
    def input_bam = task.ext.input_bam ?: "--bam $bam"
    def ref_genome = task.ext.ref_genome ?: "--reference $genome_fasta"
    def ref_gtf = task.ext.ref_genome ?: "--genedb $annotation_gtf"

    """
    touch ${prefix}.bam
    touch ${prefix}.read_assignments.tsv.gz
    touch ${prefix}.corrected_reads.bed.gz
    touch ${prefix}.transcript_tpm.tsv
    touch ${prefix}.transcript_counts.tsv
    touch ${prefix}.gene_tpm.tsv
    touch ${prefix}.gene_counts.tsv
    touch ${prefix}.transcript_models.gtf
    touch ${prefix}.transcript_model_reads.tsv.gz
    touch ${prefix}.transcript_model_tpm.tsv
    touch ${prefix}.transcript_model_counts.tsv
    touch ${prefix}.extended_annotation.gtf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isoquant: \$(isoquant --version)
    END_VERSIONS
    """
}
