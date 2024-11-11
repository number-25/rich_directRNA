process FLAIRCOLLAPSE {
    tag "$meta.id"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/flair:2.0.0--pyhdfd78af_1':
        'biocontainers/flair:2.0.0--pyhdfd78af_1' }"

    input:
    // nf-core: Where applicable please provide/convert compressed files as input/output
    tuple val(meta), path(mapped_bed)
    path(annotation_gtf)
    path(genome_fasta)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    flair \\
        correct \\
        -q ${mapped_bed} \\
        -f ${annotation_bed} \\
        -g ${genome_fasta} \\
        --threads $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flaircollapse: \$(flair --version |& sed 's/flair //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flaircollapse: \$(flair --version |& sed 's/flair //')
    END_VERSIONS
    """
}
