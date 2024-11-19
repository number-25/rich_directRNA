process ALFRED {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/alfred:0.3.2--hf9970c3_0':
        'biocontainers/alfred:0.3.2--hf9970c3_0' }"

    input:
    // Where applicable please provide/convert compressed files as input/output
    //               e.g. "*.fastq.gz" and NOT "*.fastq", "*.bam" and NOT "*.sam" etc.
    tuple val(meta), path(bam)
    path(genome_fasta)

    output:
    tuple val(meta), path("*.tsv.gz"), emit: alfred_stats
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.replicate}_alfred"

    """
    alfred qc \\
        -r ${genome_fasta} \\
        -s \\
        -o ${prefix}.tsv.gz \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        alfred: \$(alfred --version |& sed '1!d ; s/Alfred version: v//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.replicate}_alfred"
    """
    touch ${prefix}.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        alfred: \$(alfred --version |& sed '1!d ; s/Alfred version: v//')
    END_VERSIONS
    """
}
