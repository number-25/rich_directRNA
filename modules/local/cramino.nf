process CRAMINO {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cramino:0.15.0--h2e7e638_0':
        'quay.io/biocontainers/cramino:0.15.0--h2e7e638_0' }"

    input:
    // TODO nf-core: Where applicable please provide/convert compressed files as input/output
    //               e.g. "*.fastq.gz" and NOT "*.fastq", "*.bam" and NOT "*.sam" etc.
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.stats"), emit: cramino_stats
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--spliced'
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.replicate}_cramino"
    """
    cramino \\
        -t $task.cpus \\
        $args \\
        $bam \\
        > ${prefix}.stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cramino: \$(cramino --version |& sed 's/cramino //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: '--spliced'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cramino: \$(cramino --version |& sed 's/cramino //')
    END_VERSIONS
    """
}
