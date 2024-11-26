process GFFREAD {
    tag "$meta.id"
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffutils:0.13--pyh7cba7a3_0':
        'quay.io/sangerpathogens/gffutils:0.13' }"

    input:
    tuple val(meta), path(gtf)
    path(genome_fasta)

    output:
    tuple val(meta), path("*.fa"), emit: transcripts_fa
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.replicate}_transcripts"
    """
    gffread \\
    -w ${prefix}.fa \\
    -g ${genome_fasta} \\
    ${gtf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version |& sed '1!d ; s/ //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.replicate}_transcripts"
    """
    touch ${prefix}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version |& sed '1!d ; s/ //')
    END_VERSIONS
    """
}
