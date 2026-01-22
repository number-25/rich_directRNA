process SYLPH_PROFILE {
    tag "SYLPH_PROFILE"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sylph:0.9.0--ha6fb395_0' :
        'biocontainers/sylph:0.9.0--ha6fb395_0' }"

    input:
    tuple val(meta), path(sample)
    path sylph_database

    output:
    tuple val(meta), path("*.tsv")  , emit: sylph_profile
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.replicate}"

    """
    sylph profile \\
        $sylph_database \\
        $sample \\
        $task.cpus \\
        -o ${prefix}_sylph.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph_profile: \$(sylph --version | cut -d" " -f2')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.replicate}"

    """
    touch ${prefix}_sylph.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph_profile: \$(sylph --version | cut -d" " -f2')
    END_VERSIONS
    """
}
